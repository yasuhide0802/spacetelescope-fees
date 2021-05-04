#! /usr/bin/env python

import numpy as np

from jwst.lib import reffile_utils
from jwst.lib import pipe_utils

from ..stpipe import Step
from .. import datamodels

from stcal.ramp_fitting.ramp_fit import ramp_fit
from stcal.ramp_fitting.ramp_fit import BUFSIZE

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["RampFitStep"]


def get_ref_subs(model, readnoise_model, gain_model, nframes):
    """
    Get readnoise array for calculation of variance of noiseless ramps, and
    the gain array in case optimal weighting is to be done. The returned
    readnoise has been multiplied by the gain.

    Parameters
    ----------
    model : data model
        input data model, assumed to be of type RampModel

    readnoise_model : instance of data Model
        readnoise for all pixels

    gain_model : instance of gain Model
        gain for all pixels

    nframes : int
        number of frames averaged per group; from the NFRAMES keyword. Does
        not contain the groupgap.

    Returns
    -------
    readnoise_2d : float, 2D array
        readnoise subarray

    gain_2d : float, 2D array
        gain subarray
    """
    if reffile_utils.ref_matches_sci(model, gain_model):
        gain_2d = gain_model.data
    else:
        log.info('Extracting gain subarray to match science data')
        gain_2d = reffile_utils.get_subarray_data(model, gain_model)

    if reffile_utils.ref_matches_sci(model, readnoise_model):
        readnoise_2d = readnoise_model.data.copy()
    else:
        log.info('Extracting readnoise subarray to match science data')
        readnoise_2d = reffile_utils.get_subarray_data(model, readnoise_model)

    return readnoise_2d, gain_2d


def compute_int_times(input_model):
    """
    input_model: RampModel
        Compute integration time based on time series observetion and int_times.
    """
    int_times = None
    if pipe_utils.is_tso(input_model) and hasattr(input_model, 'int_times'):
        int_times = input_model.int_times

    return int_times


def create_output_image_model(input_model, image_info):
    """
    Creates output ImageModel using keywords from the inut RampModel and the
    computed arrays from the ramp fitting.

    Parameter
    ---------
    input_model: RampModel
        The input RampModel to be fitted.  It contains all needed keywords.

    image_info: tuple
        This contains the fitted arrays for the output image.
        (data, dq, var_poisson, var_rnoise, err)

    Return
    ------
    image_model: ImageModel
        An image model created from the input ramp model keywords and the
        computed ramp fitting arrays.
    """
    nints, ngroups, nrows, ncols = input_model.data.shape
    imshape = (nrows, ncols)

    # Create image model
    image_model = datamodels.ImageModel(
        data=np.zeros(imshape, dtype=np.float32),
        dq=np.zeros(imshape, dtype=np.float32),
        var_poisson=np.zeros(imshape, dtype=np.float32),
        var_rnoise=np.zeros(imshape, dtype=np.float32),
        var_err=np.zeros(imshape, dtype=np.float32))

    # Update all keys from input
    image_model.update(input_model)

    # Set image arrays
    data, dq, var_poisson, var_rnoise, err = image_info
    image_model.data = data
    image_model.dq = dq
    image_model.var_poisson = var_poisson
    image_model.var_rnoise = var_rnoise
    image_model.err = err

    return image_model


def create_optional_results_model(opt_info):
    """
    Creates the optional RampFitOutputModel from the computed arrays from ramp
    fitting

    Parameter
    ---------
    opt_info: tuple
        A tuple containing the optional results arrays from ramp fitting.

    Return
    ------
    opt_model: RampFitOutputModel
        The optional RampFitOuptput model computed during ramp fitting.
    """
    slope, sigslope, var_poisson, var_rnoise, yint = opt_info[:5]
    sigyint, pedestal, weights, crmag = opt_info[5:]
    opt_model = datamodels.RampFitOutputModel(
        slope=slope,
        sigslope=sigslope,
        var_poisson=var_poisson,
        var_rnoise=var_rnoise,
        yint=yint,
        sigyint=sigyint,
        pedestal=pedestal,
        weights=weights,
        crmag=crmag)

    return opt_model


def create_output_integration_model(input_model, int_info):
    """
    Creates the integration CubeModel from the computed arrays from ramp
    fittingl

    Parameter
    ---------
    input_model: RampModel
        The input RampModel to be fitted.  It contains all needed keywords.

    int_info: tuple
        A tuple containing the integration arrays from ramp fitting.

    Return
    ------
    int_model: CubeModel
        The integration Cubemodel computed during ramp fitting.
    """
    nints, ngroups, nrows, ncols = input_model.data.shape

    cubeshape = (nints, nrows, ncols)
    int_model = datamodels.CubeModel(
        data=np.zeros(cubeshape, dtype=np.float32),
        err=np.zeros(cubeshape, dtype=np.float32),
        dq=np.zeros(cubeshape, dtype=np.float32),
        var_poisson=np.zeros(cubeshape, dtype=np.float32),
        var_rnoise=np.zeros(cubeshape, dtype=np.float32))
    int_model.int_times = None

    # Get keywords from input model
    int_model.update(input_model)


    data, dq, var_poisson, var_rnoise, int_times, err = int_info
    int_model.data = data
    int_model.err = err
    int_model.dq = dq
    int_model.var_poisson = var_poisson
    int_model.var_rnoise = var_rnoise
    int_model.int_times = int_times

    return int_model


class RampFitStep (Step):

    """
    This step fits a straight line to the value of counts vs. time to
    determine the mean count rate for each pixel.
    """

    spec = """
        int_name = string(default='')
        save_opt = boolean(default=False) # Save optional output
        opt_name = string(default='')
        maximum_cores = option('none', 'quarter', 'half', 'all', default='none') # max number of processes to create
    """

    # Prior to 04/26/17, the following were also in the spec above:
    #      algorithm = option('OLS', 'GLS', default='OLS') # 'OLS' or 'GLS'
    #      weighting = option('unweighted', 'optimal', default='unweighted') \
    #      # 'unweighted' or 'optimal'
    # As of 04/26/17, the only allowed algorithm is 'ols', and the
    #      only allowed weighting is 'optimal'.

    algorithm = 'ols'      # Only algorithm allowed for Build 7.1
#    algorithm = 'gls'       # 032520

    weighting = 'optimal'  # Only weighting allowed for Build 7.1

    reference_file_types = ['readnoise', 'gain']

    def process(self, input):

        with datamodels.RampModel(input) as input_model:
            max_cores = self.maximum_cores
            readnoise_filename = self.get_reference_file(input_model, 'readnoise')
            gain_filename = self.get_reference_file(input_model, 'gain')

            log.info('Using READNOISE reference file: %s', readnoise_filename)
            readnoise_model = datamodels.ReadnoiseModel(readnoise_filename)
            log.info('Using GAIN reference file: %s', gain_filename)
            gain_model = datamodels.GainModel(gain_filename)

            # Try to retrieve the gain factor from the gain reference file.
            # If found, store it in the science model meta data, so that it's
            # available later in the gain_scale step, which avoids having to
            # load the gain ref file again in that step.
            if gain_model.meta.exposure.gain_factor is not None:
                input_model.meta.exposure.gain_factor = \
                    gain_model.meta.exposure.gain_factor

            log.info('Using algorithm = %s' % self.algorithm)
            log.info('Using weighting = %s' % self.weighting)

            buffsize = BUFSIZE
            if self.algorithm == "GLS":
                buffsize //= 10

            # TODO do subarray stuff here
            nframes = input_model.meta.exposure.nframes
            gain_2d, readnoise_2d = get_ref_subs(
                    input_model, readnoise_model, gain_model, nframes)

            # Save old value in case needed here because the function can return NoneType
            old_int_times = input_model.int_times
            input_model.int_times = compute_int_times(input_model)

            # The out_model and int_model are JWST data models, but need to be
            # converted to simple arrays and the models created here, not in
            # the ramp fitting code.
            image_info, int_info, opt_info, gls_opt_model = ramp_fit(
                input_model, buffsize, self.save_opt, 
                readnoise_2d, gain_2d, 
                self.algorithm, self.weighting, max_cores
            )

            readnoise_model.close()
            gain_model.close()

        # Save the GLS optional fit product, if it exists (NOT USED RIGHT NOW)
        # When GLS is implemented, this will not be a data model
        if gls_opt_model is not None:
            self.save_model(
                gls_opt_model, 'fitoptgls', output_file=self.opt_name
            )

        # Save the OLS optional fit product, if it exists
        if opt_info is not None:
            opt_model = create_optional_results_model(opt_info)
            self.save_model(opt_model, 'fitopt', output_file=self.opt_name)

        if image_info is not None:
            input_model.int_times = old_int_times
            out_model = create_output_image_model(input_model, image_info)
            out_model.meta.bunit_data = 'DN/s'
            out_model.meta.bunit_err = 'DN/s'
            out_model.meta.cal_step.ramp_fit = 'COMPLETE'
            if (input_model.meta.exposure.type in ['NRS_IFU', 'MIR_MRS']) or (
                input_model.meta.exposure.type in ['NRS_AUTOWAVE', 'NRS_LAMP'] and
                    input_model.meta.instrument.lamp_mode == 'IFU'):

                out_model = datamodels.IFUImageModel(out_model)


        if int_info is not None:
            int_model = create_output_integration_model(input_model, int_info)
            int_model.meta.bunit_data = 'DN/s'
            int_model.meta.bunit_err = 'DN/s'
            int_model.meta.cal_step.ramp_fit = 'COMPLETE'

        return out_model, int_model
