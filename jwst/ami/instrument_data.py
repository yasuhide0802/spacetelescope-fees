#
#  Module for defining data format, wavelength info, an mask geometry for these
#   instrument: NIRISS AMI
#

import logging
import numpy as np
from scipy.integrate import simpson as simps

from .mask_definitions import NRM_mask_definitions
from . import utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

um = 1.0e-6


class NIRISS:
    def __init__(self, filt,
                 objname="obj",
                 src="A0V",
                 chooseholes=None,
                 affine2d=None,
                 bandpass=None,
                 usebp=True,
                 firstfew=None,
                 rotsearch_parameters=None,
                 oversample=None,
                 psf_offset=None,
                 **kwargs):
        """
        Short Summary
        ------------
        Initialize NIRISS class
        Either user has webbpsf and filter file can be read, or this will use a
          tophat and give a warning

        Parameters
        ----------
        filt: string
            filter name

        objname: string
            name of object

        src: string
            spectral type

        chooseholes: list
            None, or e.g. ['B2', 'B4', 'B5', 'B6'] for a four-hole mask

        affine2d: Affine2d object
            Affine2d object

        bandpass: list
            None or [(wt,wlen),(wt,wlen),...].  Monochromatic would be e.g. [(1.0, 4.3e-6)]
            Explicit bandpass arg will replace *all* niriss filter-specific variables with
            the given bandpass, so you can simulate 21cm psfs through something called "F430M"!
        """
        self.chooseholes = chooseholes
        self.objname = objname
        self.filt = filt
        self.throughput = bandpass
        self.firstfew = firstfew
        if firstfew is not None:
            log.info(f'Analyzing only the first {firstfew:d} integrations')


        self.lam_c, self.lam_w = utils.get_cw_beta(self.throughput)
        ####### 


        self.wls = [self.throughput, ]
        # Wavelength info for NIRISS bands F277W, F380M, F430M, or F480M
        self.wavextension = ([self.lam_c,], [self.lam_w,])
        self.nwav = 1

        # only one NRM on JWST:
        self.telname = "JWST"
        self.instrument = "NIRISS"
        self.arrname = "jwst_g7s6c"
        self.holeshape = "hex"
        self.mask = NRM_mask_definitions(maskname=self.arrname, chooseholes=chooseholes, holeshape=self.holeshape)
        pscale_deg = utils.degrees_per_pixel(input_model)
        self.pscale_rad = np.deg2rad(pscale_deg)
        dim = input_model.data.shape[-1] # 80 pixels

        # save affine deformation of pupil object or create a no-deformation object.
        # We apply this when sampling the PSF, not to the pupil geometry.
        # This will set a default Ideal or a measured rotation, for example,
        # and include pixel scale changes due to pupil distortion.
        # Separating detector tilt pixel scale effects from pupil distortion effects is
        # yet to be determined... see comments in Affine class definition.
        if affine2d is None:
            self.affine2d = utils.Affine2d(mx=1.0, my=1.0,
                               sx=0.0, sy=0.0,
                               xo=0.0, yo=0.0, name="Ideal")
        else:
            self.affine2d = affine2d

        # finding centroid from phase slope only considered cv_phase data
        # when cv_abs data exceeds this cvsupport_threshold.
        # Absolute value of cv data normalized to unity maximum
        # for the threshold application.
        # Data reduction gurus: tweak the threshold value with experience...
        # Gurus: tweak cvsupport with use...
        self.cvsupport_threshold = {"F277W": 0.02, "F380M": 0.02, "F430M": 0.02, "F480M": 0.02}
        self.threshold = self.cvsupport_threshold[filt]


    def set_pscale(self, pscalex_deg=None, pscaley_deg=None):
        """
        Short Summary
        ------------
        Override pixel scale in header

        Parameters
        ----------
        pscalex_deg: float, degrees
            pixel scale in x-direction

        pscaley_deg: float, degrees
            pixel scale in y-direction

        Returns
        -------
        None

        """
        if pscalex_deg is not None:
            self.pscalex_deg = pscalex_deg
        if pscaley_deg is not None:
            self.pscaley_deg = pscaley_deg
        self.pscale_mas = 0.5 * (pscalex_deg + pscaley_deg) * (60 * 60 * 1000)
        self.pscale_rad = utils.mas2rad(self.pscale_mas)

    def read_data_model(self, input_model):
        """
        Short Summary
        -------------
        Retrieve info from input data model and store in NIRISS class

        Parameters
        ----------
        input_model: instance Data Model
            DM object for input

        Returns
        -------
        Data and parameters from input data model
        """
        # The info4oif_dict will get pickled to disk when we write txt files of results.
        # That way we don't drag in objects like instrument_data into code that reads text results
        # and writes oifits files - a simple built-in dictionary is the only object used in this transfer.
 

        # # Whatever we did set is averaged for isotropic pixel scale here
        # self.pscale_mas = 0.5 * (pscalex_deg + pscaley_deg) * (60 * 60 * 1000)
        # self.pscale_rad = utils.mas2rad(self.pscale_mas)

        # the following is done by updatewithheaderinfo in implaneia, 
        # don't need to pickle a dict here though.
        # all instrumentdata attributes will be available when oifits files written out?
        scidata = input_model.data.copy()
        bpdata = input_model.dq.copy()
        # make dq mask again? for now
        DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]
        JUMP_DET = dqflags.pixel["JUMP_DET"]
        dq_dnu = bpdata & DO_NOT_USE == DO_NOT_USE
        dq_jump = bpdata & JUMP_DET == JUMP_DET
        dqmask = dq_dnu | dq_jump

        pscale_deg = utils.degrees_per_pixel(input_model)
        self.pscale_rad = np.deg2rad(pscale_deg)

        self.pav3 = input_model.meta.pointing.pa_v3
        self.vparity = input_model.meta.wcsinfo.vparity
        self.v3iyang = input_model.meta.wcsinfo.v3yangle
        self.parangh = input_model.meta.wcsinfo.roll_ref
        self.ra = input_model.meta.target.ra
        self.dec = input_model.meta.target.dec
        self.crpix1 = input_model.meta.wcsinfo.crpix1
        self.crpix2 = input_model.meta.wcsinfo.crpix2
        self.pupil = input_model.meta.instrument.pupil

        datestr = input_model.meta.visit.start_time.replace(' ','T')
        self.date = datestr # is this the right start time?
        self.year = datestr[:4]
        self.month = datestr[5:7]
        self.day = datestr[8:10]
        effinttm = input_model.meta.exposure.effective_exposure_time
        nints = input_model.meta.exposure.nints
        # if 2d input, model has already been expanded to 3d, so check 0th dimension
        if input_model.data.shape[0] == 1: 
            self.itime = effinttm*nints # CHECK THIS
        else:
            self.itime = effinttm
            if self.firstfew is not None:
                if scidata.shape[0] > self.firstfew:
                    scidata = scidata[:self.firstfew, :, :]
                    dqmask = dqmask[:self.firstfew, :, :]
            self.nwav = nints # update from 1 (number of slices)
            [self.wls.append(self.wls[0]) for f in range(self.nwav-1)]

        # Get integer-pixel position of target from siaf & header info
        siaf = pysiaf.Siaf('NIRISS')
        # select AMI aperture by name
        xoffset, yoffset = input_model.meta.dither.x_offset, input_model.meta.dither.y_offset
        apername = input_model.meta.aperture.name
        nis_ami = siaf[apername]
        xtarg_detpx, ytarg_detpx = nis_ami.idl_to_sci(xoffset, yoffset) # decimal pixel position in subarray, 1 indexed?
        self.peak0, self.peak1 = int(np.floor(xtarg_detpx)), int(np.floor(ytarg_detpx))

        # do all the stuff with rotating centers, save as attributes instead of info4oif_dict
        ctrs_sky = self.mast2sky()
        oifctrs = np.zeros(self.mask.ctrs.shape)
        oifctrs[:,0] = ctrs_sky[:,1].copy() * -1
        oifctrs[:,1] = ctrs_sky[:,0].copy() * -1
        self.ctrs_eqt = oifctrs
        self.ctrs_inst = self.mask.ctrs
        self.hdia = self.mask.hdia
        self.nslices = self.nwav

        # Trim refpix from all slices
        scidata = scidata[:,4:, :]
        dqmask = dqmask[:,4:, :]

        self.rootfn = input_model.meta.filename.replace('.fits','')

        

        # all info needed to write out oifits should be stored in NIRISS object attributes
        return scidata, dqmask

    def reset_nwav(self, nwav):
        """
        Reset self.nwav

        Parameters
        ----------
        nwav: integer
            length of axis3 for 3D input

        Returns
        -------
        """
        self.nwav = nwav

#####
    def mast2sky(self):
        """
        Rotate hole center coordinates:
            Clockwise by the V3 position angle - V3I_YANG from north in degrees if VPARITY = -1
            Counterclockwise by the V3 position angle - V3I_YANG from north in degrees if VPARITY = 1
        Hole center coords are in the V2, V3 plane in meters.
        Return rotated coordinates to be put in info4oif_dict.
        implane2oifits.ObservablesFromText uses these to calculate baselines.
        """
        pa = self.pav3
        mask_ctrs = copy.deepcopy(self.mask.ctrs)
        # rotate by an extra 90 degrees (RAC 9/21)
        # these coords are just used to orient output in OIFITS files
        # NOT used for the fringe fitting itself
        mask_ctrs = utils.rotate2dccw(mask_ctrs,np.pi/2.)
        vpar = self.vparity # Relative sense of rotation between Ideal xy and V2V3
        v3iyang = self.v3i_yang
        rot_ang = pa - v3iyang # subject to change!

        if pa != 0.0:
            # Using rotate2sccw, which rotates **vectors** CCW in a fixed coordinate system,
            # so to rotate coord system CW instead of the vector, reverse sign of rotation angle.  Double-check comment
            if vpar == -1:
                # rotate clockwise  <rotate coords clockwise?>
                ctrs_rot = utils.rotate2dccw(mask_ctrs, np.deg2rad(-rot_ang))
                log.debug(f'Rotating mask hole centers clockwise by {rot_ang:.3f} degrees')
            else:
                # counterclockwise  <rotate coords counterclockwise?>
                ctrs_rot = utils.rotate2dccw(mask_ctrs, np.deg2rad(rot_ang))
                log.debug(f'Rotating mask hole centers counterclockwise by {rot_ang:.3f} degrees')
        else:
            ctrs_rot = mask_ctrs
        return ctrs_rot


#####