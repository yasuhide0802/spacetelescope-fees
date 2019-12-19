import time
import logging

import numpy as np
from ..datamodels import dqflags
from ..lib import reffile_utils
from . import twopoint_difference as twopt
from . import yintercept as yint
import multiprocessing

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def detect_jumps (input_model, gain_model, readnoise_model,
                  rejection_threshold, do_yint, signal_threshold, max_cores,
                  max_jump_to_flag_neighbors, min_jump_to_flag_neighbors,
                  flag_4_neighbors):
    """
    This is the high-level controlling routine for the jump detection process.
    It loads and sets the various input data and parameters needed by each of
    the individual detection methods and then calls the detection methods in
    turn.

    Note that the detection methods are currently setup on the assumption
    that the input science and error data arrays will be in units of
    electrons, hence this routine scales those input arrays by the detector
    gain. The methods assume that the read noise values will be in units
    of DN.

    The gain is applied to the science data and error arrays using the
    appropriate instrument- and detector-dependent values for each pixel of an
    image.  Also, a 2-dimensional read noise array with appropriate values for
    each pixel is passed to the detection methods.
    """
    if max_cores is None:
        numslices = 1
    else:
        num_cores = multiprocessing.cpu_count()
        log.info("Found %d possible cores to use for jump detection " % num_cores)
        if max_cores == 'quarter':
            numslices = num_cores // 4 or 1
        elif max_cores == 'half':
            numslices = num_cores // 2 or 1
        elif max_cores == 'all':
            numslices = num_cores
        else:
            numslices = 1

    # Load the data arrays that we need from the input model
    output_model = input_model.copy()
    data = input_model.data
    err  = input_model.err
    gdq  = input_model.groupdq
    pdq  = input_model.pixeldq


    # Get 2D gain and read noise values from their respective models
    if reffile_utils.ref_matches_sci(input_model, gain_model):
        gain_2d = gain_model.data
    else:
        log.info('Extracting gain subarray to match science data')
        gain_2d = reffile_utils.get_subarray_data(input_model, gain_model)

    if reffile_utils.ref_matches_sci(input_model, readnoise_model):
        readnoise_2d = readnoise_model.data
    else:
        log.info('Extracting readnoise subarray to match science data')
        readnoise_2d = reffile_utils.get_subarray_data(input_model, readnoise_model)

    # Flag the pixeldq where the gain is <=0 or NaN so they will be ignored
    wh_g = np.where( gain_2d <= 0.)
    if len(wh_g[0] > 0):
        pdq[wh_g] = np.bitwise_or( pdq[wh_g], dqflags.pixel['NO_GAIN_VALUE'] )
        pdq[wh_g] = np.bitwise_or( pdq[wh_g], dqflags.pixel['DO_NOT_USE'] )

    wh_g = np.where( np.isnan( gain_2d ))
    if len(wh_g[0] > 0):
        pdq[wh_g] = np.bitwise_or( pdq[wh_g], dqflags.pixel['NO_GAIN_VALUE'] )
        pdq[wh_g] = np.bitwise_or( pdq[wh_g], dqflags.pixel['DO_NOT_USE'] )

    # Apply gain to the SCI, ERR, and readnoise arrays so they're in units
    # of electrons

    data *= gain_2d
    err  *= gain_2d
    readnoise_2d *= gain_2d

    # Apply the 2-point difference method as a first pass
    log.info('Executing two-point difference method')
    start = time.time()
    nrows = data.shape[-2]
    ncols = data.shape[-1]
    num_groups = data.shape[1]
    num_ints = data.shape[0]
    median_slopes = np.zeros((num_ints, nrows, ncols), dtype=np.float32)
    all_ratios = np.zeros((num_ints, nrows, ncols, num_groups-1), dtype=np.float32)
    yincrement = int(nrows / numslices)
    print("yincrement ", yincrement)
    print("num slices ",numslices)
    slices = []
    # Slice up data, gdq, readnoise_2d into slices
    # Each element of slices is a tuple of
    # (data, gdq, readnoise_2d, rejection_threshold, nframes)
    for i in range(numslices - 1):
        slices.insert(i, (data[:, :, i * yincrement:(i + 1) * yincrement, :],
                          gdq[:, :, i * yincrement:(i + 1) * yincrement, :],
                          readnoise_2d[i * yincrement:(i + 1) * yincrement, :],
                          rejection_threshold, num_groups))
    # last slice get the rest
    slices.insert(numslices - 1, (data[:, :, (numslices - 1) * yincrement:nrows, :],
                                 gdq[:, :, (numslices - 1) * yincrement:nrows, :],
                                 readnoise_2d[(numslices - 1) * yincrement:nrows, :],
                                 rejection_threshold, num_groups))
    if numslices == 1:
        median_slopes, gdq, all_ratios = twopt.find_crs(data, gdq, readnoise_2d, rejection_threshold, num_groups)
        elapsed = time.time() - start
    else:
        log.info("Creating %d processes for jump detection " % numslices)
        pool = multiprocessing.Pool(processes=numslices)
        real_result = pool.starmap(twopt.find_crs, slices)
        k = 0
        # Reconstruct median_slopes and gdq from the slice results
        for resultslice in real_result:
            print("k ",k,"yincrement ",yincrement)
            if len(real_result) == k + 1:  # last result
                median_slopes[:, k * yincrement:nrows, :] = resultslice[0]
                gdq[:, :, k * yincrement:nrows, :] = resultslice[1]
                all_ratios[:, k * yincrement:nrows, :, :] = resultslice[2]
            else:
                print("k ", k, "yincrement ", yincrement,resultslice[0].shape, resultslice[1].shape, resultslice[2].shape)
                median_slopes[:, k * yincrement:(k + 1) * yincrement, :] = resultslice[0]
                gdq[:, :, k * yincrement:(k + 1) * yincrement, :] = resultslice[1]
                all_ratios[:,  k * yincrement:(k + 1) * yincrement, :, :] = resultslice[2]
            k += 1
        pool.terminate()
        pool.close()
        elapsed = time.time() - start
    log.info('Elapsed time of twopoint difference = %g sec' % elapsed)
    if flag_4_neighbors:
        cr_int, cr_group, cr_row, cr_col = np.where(np.bitwise_and(gdq, dqflags.group['JUMP_DET']))
        number_pixels_with_cr = len(cr_int)
        for j in range(number_pixels_with_cr):
            if all_ratios[cr_int[j],  cr_row[j], cr_col[j], cr_group[j]-1] < max_jump_to_flag_neighbors and \
                all_ratios[cr_int[j], cr_row[j], cr_col[j], cr_group[j]-1] > min_jump_to_flag_neighbors:
                    if cr_row[j] != 0:
                        gdq[cr_int[j], cr_group[j], cr_row[j]-1, cr_col[j]] = np.bitwise_or(
                            gdq[cr_int[j], cr_group[j], cr_row[j]-1, cr_col[j]],
                            dqflags.group['JUMP_DET'])
                    if cr_row[j] != nrows-1:
                        gdq[cr_int[j], cr_group[j], cr_row[j] + 1, cr_col[j]] = np.bitwise_or(
                            gdq[cr_int[j], cr_group[j], cr_row[j] + 1, cr_col[j]],
                            dqflags.group['JUMP_DET'])
                    if cr_col[j] != 0:
                        gdq[cr_int[j], cr_group[j], cr_row[j], cr_col[j]- 1] = np.bitwise_or(
                            gdq[cr_int[j], cr_group[j], cr_row[j], cr_col[j]-1],
                            dqflags.group['JUMP_DET'])
                    if cr_col[j] != ncols -1:
                        gdq[cr_int[j], cr_group[j], cr_row[j], cr_col[j] + 1] = np.bitwise_or(
                            gdq[cr_int[j], cr_group[j], cr_row[j], cr_col[j] + 1],
                            dqflags.group['JUMP_DET'])

    elapsed = time.time() - start
    log.info('Total elapsed time = %g sec' % elapsed)

    # Apply the y-intercept method as a second pass, if requested
    if do_yint:

        # Set up the ramp time array for the y-intercept method
        group_time = output_model.meta.exposure.group_time
        times = np.array([(k+1)*group_time for k in range(num_groups)])
        median_slopes /= group_time

        # Now apply the y-intercept method
        log.info('Executing yintercept method')
        start = time.time()
        yint.find_crs(data, err, gdq, times, readnoise_2d,
                        rejection_threshold, signal_threshold, median_slopes)
        elapsed = time.time() - start
        log.debug('Elapsed time = %g sec' % elapsed)

    # Update the DQ arrays of the output model with the jump detection results
    output_model.groupdq = gdq
    output_model.pixeldq = pdq

    return output_model
