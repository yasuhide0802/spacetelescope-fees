#  Module for undersampling correction
#
import logging
import numpy as np

from stdatamodels.jwst.datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def undersampling_correction(input_model, signal_threshold):
    """
    Correct for undersampling.

    Parameters
    ----------
    input_model : `~jwst.datamodels.RampModel`
        The input science data to be corrected

    signal_threshold : float
        Science value above which a group will be flagged as UNDERSAMP

    Returns
    -------
    output_model : `~jwst.datamodels.RampModel`
        Data model with undersampling_correction applied; add UNDERSAMP flag
        to groups exceeding signal_threshold
    """
    data = input_model.data
    gdq = input_model.groupdq

    # Create the output model as a copy of the input
    output_model = input_model.copy()

    log.info('Using signal_threshold: %.2f', signal_threshold)

    gdq_new = flag_pixels(data, gdq, signal_threshold)

    # Save the flags in the output GROUPDQ array
    output_model.groupdq = gdq_new

    return output_model


def flag_pixels(data, gdq, signal_threshold):
    """
    Flag first group in each ramp that exceeds signal_threshold as UNDERSAMP and DO_NOT_USE,
    skipping groups already flagged as DO_NOT_USE; then flag all subsequent groups in the ramp.

    Parameters
    ----------
    data : float, 4D array
        science array

    gdq : int, 4D array
        group dq array

    signal_threshold : float
        Science value above which a group will be flagged as UNDERSAMP and DO_NOT_USE

    Returns
    -------
    gdq : int, 4D array
        updated group dq array
    """
    n_ints, n_grps, n_rows, n_cols = gdq.shape
    num_pix = n_cols * n_rows

    lowest_exc_1d = np.zeros(num_pix) + n_grps

    ncols = data.shape[3]
    nrows = data.shape[2]

    # Loop over all groups to locate exceedances
    for ii_int in range(n_ints):
        for ii_grp in range(n_grps):
            data_1d = data[ii_int, ii_grp, :, :].reshape(num_pix)  # vectorize slice
            gdq_1d = gdq[ii_int, ii_grp, :, :].reshape(num_pix)

            wh_not_dnu = np.logical_not(gdq_1d & dqflags.group['DO_NOT_USE'])

            # In the current group for all ramps, locate pixels that :
            #  a) exceed the signal_threshold, and
            #  b) have not been previously flagged as an exceedance, and
            #  c) were not flagged in an earlier step as DO_NOT_USE
            wh_exc_1d = np.where((data_1d > signal_threshold) &
                                 (lowest_exc_1d == n_grps) & wh_not_dnu)

            # ... and mark those pixels, as current group is their first exceedance
            if len(wh_exc_1d[0] > 0):  # For ramps previously unflagged ...
                lowest_exc_1d[wh_exc_1d] = ii_grp

    # Loop over exceedances and flag as UNDERSAMP and DO_NOT_USE, then flag
    #   that pixel's 4 nearest neighbors
    lowest_exc_2d = lowest_exc_1d.reshape((n_rows, n_cols))

    for ii_int in range(n_ints):
        for ii_grp in range(n_grps):
            wh_set_flag = np.where(lowest_exc_2d == ii_grp)
            # set arrays of coordinates of each exceedance
            yy = wh_set_flag[0]
            xx = wh_set_flag[1]

            if len(wh_set_flag) > 0 and len(xx) > 0:  # there are exceedances
                gdq[ii_int, ii_grp:, yy, xx] = \
                    np.bitwise_or(gdq[ii_int, ii_grp:, yy, xx],
                                  dqflags.group['UNDERSAMP']
                                  | dqflags.group['DO_NOT_USE'])

                # Set the same flags for the 4 nearest neighbors
                wh_west = np.where(xx > 0)  # pixels whose western neighbors to be flagged
                if len(wh_west[0]) > 0:
                    gdq[ii_int, ii_grp:, yy, xx[wh_west] - 1] = \
                        np.bitwise_or(gdq[ii_int, ii_grp:, yy, xx[wh_west] - 1],
                                      dqflags.group['UNDERSAMP'] |
                                      dqflags.group['DO_NOT_USE'])

                wh_east = np.where(xx < ncols - 1)  # pixels whose eastern neighbors to be flagged
                if len(wh_east[0]) > 0:
                    gdq[ii_int, ii_grp:, yy, xx[wh_east] + 1] = \
                        np.bitwise_or(gdq[ii_int, ii_grp:, yy, xx[wh_east] + 1],
                                      dqflags.group['UNDERSAMP'] |
                                      dqflags.group['DO_NOT_USE'])

                wh_north = np.where(yy > 0)  # pixels whose northern neighbors to be flagged
                if len(wh_north[0]) > 0:
                    gdq[ii_int, ii_grp:, yy[wh_north] - 1, xx] = \
                        np.bitwise_or(gdq[ii_int, ii_grp:, yy[wh_north] - 1, xx],
                                      dqflags.group['UNDERSAMP'] |
                                      dqflags.group['DO_NOT_USE'])

                wh_south = np.where(yy < nrows - 1)  # pixels whose southern neighbors to be flagged
                if len(wh_south[0]) > 0:
                    gdq[ii_int, ii_grp:, yy[wh_south] + 1, xx] = \
                        np.bitwise_or(gdq[ii_int, ii_grp:, yy[wh_south] + 1, xx],
                                      dqflags.group['UNDERSAMP'] |
                                      dqflags.group['DO_NOT_USE'])

    return gdq
