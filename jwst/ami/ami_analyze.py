#  Module for applying the LG-PLUS algorithm to an AMI exposure
import logging
import numpy as np

from .find_affine2d_parameters import find_rotation
from . import instrument_data
from . import nrm_core
from .utils import img_median_replace
from stdatamodels.jwst.datamodels import dqflags

from astropy import units as u

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def apply_LG_plus(input_model, 
                oversample, rotation,
                psf_offset, rotsearch_parameters, 
                src, bandpass, usebp, firstfew, chooseholes, affine2d
                # **kwargs?
                ):
    """
    Short Summary
    -------------
    Applies the image plane algorithm to an AMI image

    Parameters
    ----------
    input_model : data model object
        AMI science image to be analyzed

    # filter_model : filter model object
    #     filter throughput reference data

    oversample : integer
        Oversampling factor

    rotation : float (degrees)
        Initial guess at rotation of science image relative to model

    Returns
    -------
    output_model : Fringe model object
        Fringe analysis data

    """
    # Create copy of input_model to avoid overwriting input
    input_copy = input_model.copy()

    # If the input image is 2D, expand all relevant extensions to be 3D
    # Incl. those not currently used?
    if len(input_model.data.shape) == 2:
        input_copy.data = np.expand_dims(input_copy.data, axis=0)
        input_copy.dq = np.expand_dims(input_copy.dq, axis=0)
        # input_copy.err = np.expand_dims(input_copy.err, axis=0)
        # input_copy.var_poisson = np.expand_dims(input_copy.var_poisson, axis=0)
        # input_copy.var_rnoise = np.expand_dims(input_copy.var_rnoise, axis=0)
        # input_copy.var_flat = np.expand_dims(input_copy.var_flat, axis=0)



    # If the input data were taken in full-frame mode, extract a region
    # equivalent to the SUB80 subarray mode to make execution time acceptable.
    if input_model.meta.subarray.name.upper() == 'FULL':
        log.info("Extracting 80x80 subarray from full-frame data")
        xstart = 1045
        ystart = 1
        xsize = 80
        ysize = 80
        xstop = xstart + xsize - 1
        ystop = ystart + ysize - 1
        input_copy.data = input_copy.data[:, ystart - 1:ystop, xstart - 1:xstop]
        input_copy.dq = input_copy.dq[:, ystart - 1:ystop, xstart - 1:xstop]
        input_copy.err = input_copy.err[:, ystart - 1:ystop, xstart - 1:xstop]

    # Replace NaN's and DO_NOT_USE pixels in the input image
    # with median of surrounding pixel values in a 3x3 box
    # box_size = 3
    # input_copy = img_median_replace(input_copy, box_size)

    # Instead, run bp_fix code here??
    # TODO: run_bp_fix

    data = input_copy.data
    dim = data.shape[-1] # 80 px 

    # Create mask of DO_NOT_USE and JUMP_DET pixels
    bpdata = np.array(input_copy.dq)
    DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
    JUMP_DET = dqflags.pixel["JUMP_DET"]
    dq_dnu = bpdata & DO_NOT_USE == DO_NOT_USE
    dq_jump = bpdata & JUMP_DET == JUMP_DET
    dqmask = dq_dnu | dq_jump

    # now run bpfix, updating the data and dq array in input_copy


    # Initialize transformation parameters:
    #   mx, my: dimensionless magnifications
    #   sx, sy: dimensionless shears
    #   x0, y0: offsets in pupil space
    mx = 1.0
    my = 1.0
    sx = 0.0
    sy = 0.0
    xo = 0.0
    yo = 0.0

    psf_offset_ff = None

    # there should be an offset search?


    # get filter, pixel scale from input_model,
    # make bandpass array for find_rotation, instrument_data calls
    filt = input_copy.meta.instrument.filter
    pscale_deg = utils.degrees_per_pixel(input_copy)
    PIXELSCALE_r = np.deg2rad(pscale_deg)

    # # spec type needs to be a user input?

    holeshape = 'hex'
    # # throughput ref file is too coarsely sampled, use webbpsf data instead
    # # get throughput here instead of in instrument_data
    if bandpass is not None:
        self.log.info('User-defined bandpass provided: OVERWRITING ALL NIRISS-SPECIFIC FILTER/BANDPASS VARIABLES')
        # bandpass can be user-defined synphot object or appropriate array
        if isinstance(bandpass, synphot.spectrum.SpectralElement):
            self.log.info('User-defined synphot spectrum provided')
            wl, wt = bandpass._get_arrays(bandpass.waveset)
            bandpass = np.array((wt,wl)).T
        else:
            self.log.info('User-defined bandpass array provided')
            bandpass = np.array(bandpass)

    else:
        # get the filter and source spectrum
        self.log.info(f'Getting WebbPSF throughput data for {filt}.')
        filt_spec = utils.get_filt_spec(filt)
        self.log.info(f'Getting source spectrum for spectral type {src}.')
        src_spec = utils.get_src_spec(src) # always going to be A0V currently
        nspecbin = 19 # how many wavelngth bins used across bandpass -- affects runtime
        # **NOTE**: As of WebbPSF version 1.0.0 filter is trimmed to where throughput is 10% of peak
        # For consistency with WebbPSF simultions, use trim=0.1
        bandpass = utils.combine_src_filt(filt_spec, 
                                      src_spec, 
                                      trim=0.01, 
                                      nlambda=nspecbin,
                                      verbose=False, 
                                      plot=False) 
            

    rotsearch_d = np.append(np.arange(rotsearch_parameters[0], rotsearch_parameters[1], rotsearch_parameters[2]),
                            rotsearch_parameters[1])

    log.info(f'Initial values to use for rotation search {rotsearch_d}')
    if affine2d is None:
        # affine2d object, can be overridden by user input affine?
        affine2d = find_rotation(data[:, :, :], psf_offset, rotsearch_d,
                                 mx, my, sx, sy, xo, yo,
                                 PIXELSCALE_r, dim, bandpass, oversample, holeshape)

    niriss = instrument_data.NIRISS(filt, 
                                    bandpass=bandpass,
                                    affine2d=affine2d,
                                    src=src,
                                    firstfew=firstfew,
                                    usebp=usebp,
                                    chooseholes=chooseholes)
    # more args to pass to instrument_data: src, usebp, firstfew, chooseholes

    ff_t = nrm_core.FringeFitter(niriss, 
                                psf_offset_ff=psf_offset_ff,
                                oversample=oversample)

    #output_model = ff_t.fit_fringes_all(input_copy)
    # FOR NOW: DEBUGGING
    n_resid_arr, cp_arr, va_arr = ff_t.fit_fringes_all(input_copy)


    # Copy header keywords from input to output
    # output_model.update(input_model, only="PRIMARY")

    # return output_model
