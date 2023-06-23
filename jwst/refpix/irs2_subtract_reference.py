import logging

import numpy as np
from scipy.ndimage import convolve1d
from scipy.stats import norm

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def correct_model(input_model, irs2_model,
                  scipix_n_default=16, refpix_r_default=4, pad=8):
    """Correct an input NIRSpec IRS2 datamodel using reference pixels.

    Parameters
    ----------
    input_model: ramp model
        The input science data model.

    irs2_model: IRS2 model
        The reference file model for IRS2 correction.

    scipix_n_default: int
        Number of regular samples before stepping out to collect
        reference samples.

    refpix_r_default: int
        Number of reference samples before stepping back in to collect
        regular samples.

    pad: int
        The effective number of pixels sampled during the pause at the end
        of each row (new-row overhead).  The padding is needed to preserve
        the phase of temporally periodic signals.

    Returns
    -------
    output_model: ramp model
        The science data with reference output and reference pixels
        subtracted.
    """
    #
    """Readout parameters
    scipix_n   16         Number of regular samples before stepping out
                          to collect reference samples
    refpix_r    4         Number of reference samples before stepping back
                          in to collect regular samples
    NFOH        1 row     New frame overhead (714?)
    NROH        8 pixels  New row overhead (`pad`)
    JOH         1 pixel   Jump overhead for stepping out to or in from
                          reference pixels
    TPIX       10 microseconds  Pixel dwell time
    tframe     14.5889 s  Frame readout time

    The image and reference data will be rearranged into a 1-D array
    containing the values in time order, i.e. the element number * TPIX
    is the relative time at which a pixel was read out.  This array will
    have elements corresponding to the gaps (overheads) when no pixel was
    being read.  This 1-D array has length 1,458,176, which is equal to
    712 * 2048:

    ((scipix_n + refpix_r + 2) * (512 // scipix_n) + NROH) * 2048

    The total frame readout time is:
    (((scipix_n + refpix_r + 2) * (512 // scipix_n) + NROH) * 2048 + NFOH)
        * TPIX
    This agrees with the above value of tframe (14.5889 s) if NFOH = 714.
    """

    # Only copy in SCI data array for now; that's all we need. The rest
    # of the input model will be copied to output at the end of the step.
    data = input_model.data.copy()
    input_model.meta.cal_step.refpix = 'not specified yet'

    # Load the reference file data.
    # The reference file data are complex, but they're stored as float, with
    # alternating real and imaginary parts.  We therefore check for twice
    # as many rows as we actually want, and we'll divide that number by two
    # when allocating the arrays alpha and beta.
    nrows = len(irs2_model.irs2_table.field("alpha_0"))
    expected_nrows = 2 * 712 * 2048
    if nrows != expected_nrows:
        log.error("Number of rows in reference file = {},"
                  " but it should be {}.".format(nrows, expected_nrows))
        input_model.meta.cal_step.refpix = 'SKIPPED'
        return input_model
    alpha = np.ones((4, nrows // 2), dtype=np.complex64)
    beta = np.zeros((4, nrows // 2), dtype=np.complex64)

    alpha[0, :] = float_to_complex(irs2_model.irs2_table.field("alpha_0"))
    alpha[1, :] = float_to_complex(irs2_model.irs2_table.field("alpha_1"))
    alpha[2, :] = float_to_complex(irs2_model.irs2_table.field("alpha_2"))
    alpha[3, :] = float_to_complex(irs2_model.irs2_table.field("alpha_3"))

    beta[0, :] = float_to_complex(irs2_model.irs2_table.field("beta_0"))
    beta[1, :] = float_to_complex(irs2_model.irs2_table.field("beta_1"))
    beta[2, :] = float_to_complex(irs2_model.irs2_table.field("beta_2"))
    beta[3, :] = float_to_complex(irs2_model.irs2_table.field("beta_3"))

    scipix_n = input_model.meta.exposure.nrs_normal
    if scipix_n is None:
        log.warning("Keyword NRS_NORM not found; using default value %d" %
                    scipix_n_default)
        scipix_n = scipix_n_default

    refpix_r = input_model.meta.exposure.nrs_reference
    if refpix_r is None:
        log.warning("Keyword NRS_REF not found; using default value %d" %
                    refpix_r_default)
        refpix_r = refpix_r_default

    # Convert from sky (DMS) orientation to detector orientation.
    detector = input_model.meta.instrument.detector
    if detector == "NRS1":
        data = np.swapaxes(data, 2, 3)
    elif detector == "NRS2":
        data = np.swapaxes(data, 2, 3)[:, :, ::-1, ::-1]
    else:
        log.warning("Detector '%s'; not changing orientation (sky vs detector)"
                    % detector)

    n_int = data.shape[0]               # number of integrations in file
    ny = data.shape[-2]                 # 2048
    nx = data.shape[-1]                 # 3200

    # Create a mask that indicates the locations of normal vs interspersed
    # reference pixels. True flags normal pixels, False is reference pixels.
    irs2_mask = make_irs2_mask(nx, ny, scipix_n, refpix_r)

    '''
    from matplotlib import pyplot as plt
    import copy
    fig, axs = plt.subplots(1, 1, figsize=(12, 10))
    fig.suptitle('Original data')
    middle_of_img = nx // 2
    orig_data = copy.deepcopy(data)
    vminmax = [min(data[0][0][middle_of_img]), max(data[0][-1][middle_of_img])]
    im = axs.imshow(data[0][-1], aspect="auto", origin='lower', vmin=vminmax[0], vmax=vminmax[1])
    plt.colorbar(im, ax=axs)
    axs.tick_params(axis='both', which='both', bottom=True, top=True, right=True, direction='in', labelbottom=True)
    axs.minorticks_on()
    plt.show(block=False)
    '''

    # If the IRS2 reference file includes data quality info, use that to
    # set bad reference pixel values to zero.
    if hasattr(irs2_model, 'dq_table') and len(irs2_model.dq_table) > 0:
        output = irs2_model.dq_table.field("output")
        odd_even = irs2_model.dq_table.field("odd_even")
        mask = irs2_model.dq_table.field("mask")
        # Set interleaved reference pixel values to zero if they are flagged
        # as bad in the DQ extension of the CRDS reference file.
        clobber_ref(data, output, odd_even, mask)
        '''
        from matplotlib import pyplot as plt
        fig, axs = plt.subplots(1, 1, figsize=(12, 10))
        fig.suptitle('masked data')
        middle_of_img = nx // 2
        vminmax = [min(data[0][0][middle_of_img]), max(data[0][-1][middle_of_img])]
        im = axs.imshow(data[0][-1], aspect="auto", origin='lower', vmin=vminmax[0], vmax=vminmax[1])
        plt.colorbar(im, ax=axs)
        axs.tick_params(axis='both', which='both', bottom=True, top=True, right=True, direction='in', labelbottom=True)
        axs.minorticks_on()
        plt.show(block=False)
        fig, axs = plt.subplots(1, 1, figsize=(12, 10))
        fig.suptitle('difference')
        diff = orig_data[0][-1] - data[0][-1]
        vminmax = [min(diff[middle_of_img]), max(diff[middle_of_img])]
        im = axs.imshow(diff, aspect="auto", origin='lower', vmin=vminmax[0], vmax=vminmax[1])
        plt.colorbar(im, ax=axs)
        axs.tick_params(axis='both', which='both', bottom=True, top=True, right=True, direction='in', labelbottom=True)
        axs.minorticks_on()
        plt.show(block=True)
        '''

    else:
        log.warning("DQ extension not found in reference file")

    # Create new array to save result of IRS2 processing: slice out the reference
    # pixels from the original data array using the mask of the science pixels
    corr_data = data[..., irs2_mask==True]

    # Compute and apply the correction to one integration at a time
    for integ in range(n_int):
        log.info(f'Working on integration {integ+1}')

        data0 = data[integ, :, :, :]
        data0 = subtract_reference(data0, alpha, beta, irs2_mask, scipix_n, refpix_r, pad)
        corr_data[integ, :, :, :] = data0

    # Convert corrected data back to sky orientation
    output_model = input_model.copy()
    #temp_data = data[:, :, :, nx - ny:]
    if detector == "NRS1":
        output_model.data = np.swapaxes(corr_data, 2, 3)
    elif detector == "NRS2":
        output_model.data = np.swapaxes(corr_data[:, :, ::-1, ::-1], 2, 3)
    else:                       # don't change orientation
        output_model.data = corr_data

    # Strip interleaved ref pixels from the PIXELDQ, GROUPDQ, and ERR extensions.
    strip_ref_pixels(output_model, irs2_mask)

    return output_model


def float_to_complex(data):
    """Convert real and imaginary parts to complex"""
    nelem = len(data)

    return data[0:-1:2] + 1j * data[1:nelem:2]


def make_irs2_mask(nx, ny, scipix_n, refpix_r):

    # Number of (scipix_n + refpix_r) per output, assuming four amplifier
    # outputs and one reference output.
    irs2_nx = max((ny, nx))

    # Length of the reference output section.
    refout = irs2_nx // 5
    part = refout - (scipix_n // 2 + refpix_r)
    k = part // (scipix_n + refpix_r)
    # `part` consists of k * (scipix_n + refpix_r) + stuff_at_end
    stuff_at_end = part - k * (scipix_n + refpix_r)

    # Create the mask that flags normal pixels as True.
    irs2_mask = np.ones(irs2_nx, dtype=bool)
    irs2_mask[0:refout] = False

    # Check whether the interspersed reference pixels are in the same
    # locations regardless of readout direction.
    if stuff_at_end == scipix_n // 2:
        # Yes, they are in the same locations.
        for i in range(refout + scipix_n // 2, irs2_nx + 1,
                       scipix_n + refpix_r):
            irs2_mask[i:i + refpix_r] = False
    else:
        # Set the flags for each readout direction separately.
        nelem = refout                  # number of elements per output
        temp = np.ones(nelem, dtype=bool)
        for i in range(scipix_n // 2, nelem + 1,
                       scipix_n + refpix_r):
            temp[i:i + refpix_r] = False
        j = refout
        irs2_mask[j:j + nelem] = temp.copy()
        j += nelem
        irs2_mask[j:j + nelem] = temp[::-1].copy()
        j += nelem
        irs2_mask[j:j + nelem] = temp.copy()
        j += nelem
        irs2_mask[j:j + nelem] = temp[::-1].copy()

    return irs2_mask


def clobber_ref(data, output, odd_even, mask, scipix_n=16, refpix_r=4):
    """Set some interleaved reference pixel values to zero.

    Long Description
    ----------------
    This is an explanation of the arithmetic for computing `ref` in the loop
    over the list of bit numbers that is returned by `decode_mask`.
    Reads of reference pixels are interleaved with reads of science data.  The
    pattern of science pixels (S) and reference pixels (r) looks like this:

    SSSSSSSSrrrrSSSSSSSSSSSSSSSSrrrrSSSSSSSSSSSSSSSSrrrr ... rrrrSSSSSSSS

    Within each amplifier output, a row starts and ends with 8 (scipix_n / 2)
    science pixels, and the row contains 32 blocks of 4 reference pixels.
    There are 20 (scipix_n + refpix_r) pixels from the start of one block of
    reference pixels to the start of the next.  `k` is an integer between
    0 and 31, inclusive, an index to identify the block of reference pixels
    that we need to modify (we'll set two of the pixels to zero).  `odd_even`
    is either 1 or 2, indicating that we should set either the first or the
    second pair of reference pixels to 0.

    The same set of interleaved reference pixels will be set to 0 regardless
    of integration number, group number, or image line number.

    Parameters
    ----------
    data : 4-D ndarray
        The data array in detector orientation.  This includes both the
        science and interleaved reference pixel values.  `data` will be
        modified in-place to set some of the reference pixel values to zero.
        The science data values will not be modified.

    output : 1-D ndarray, int16
        An array of amplifier output numbers, 1, 2, 3, or 4, read from the
        OUTPUT column in the DQ extension of the CRDS reference file.

    odd_even : 1-D ndarray, int16
        An array of integer values, which may be either 1 or 2, read from the
        ODD_EVEN column in the DQ extension of the CRDS reference file.

    mask : 1-D ndarray, uint32
        The MASK column read from the CRDS reference file.

    scipix_n : int
        Number of regular (science) samples before stepping out to collect
        reference samples.

    refpix_r : int
        Number of reference samples before stepping back in to collect
        regular samples.
    """

    nx = data.shape[-1]                 # 3200
    amplifier = nx // 5                 # 640
    nrows = len(output)
    total_bad_pixels = []
    '''
    for row in range(nrows):
        # `offset` is the offset in pixels from the beginning of the row
        # to the start of the current amp output.  `offset` starts with
        # 640 in order to skip over the reference output.
        offset = output[row] * (nx // 5)                # nx // 5 is 640
        # The readout direction alternates from one amp output to the next.
        if output[row] // 2 * 2 == output[row]:
            odd_even_row = 3 - odd_even[row]            # 1 --> 2;  2 --> 1
        else:
            odd_even_row = odd_even[row]
        bits = decode_mask(output[row], mask[row])
        total_bad_pixels.append(len(bits))
        log.debug("output {}  odd_even {}  mask {}  DQ bits {}"
                  .format(output[row], odd_even[row], mask[row], bits))
        for k in bits:
            ref = (offset + scipix_n // 2 + k * (scipix_n + refpix_r) +
                   2 * (odd_even_row - 1))
            log.debug("bad interleaved reference at pixels {} {}"
                      .format(ref, ref + 1))
            data[..., ref:ref + 2] = 0.
    '''
    # loop through the reference table rows
    for row in range(nrows):
        bits = decode_mask(output[row], mask[row])
        total_bad_pixels.append(len(bits))
        counter_even, counter_odd = 0, 0
        offset = int(output[row] * amplifier)
        # only loop through the x-axis ref pixel groups when there are bad pixels
        if bits:
            for ri in range(scipix_n//2, amplifier, scipix_n+refpix_r):
                ri = ri + offset
                # set the corresponding columns of reference pixels to 0.0
                if ri % 2 != 0:
                    counter_odd += 1
                    counter = counter_odd
                else:
                    counter_even += 1
                    counter = counter_even
                for k in bits:
                    if counter == k+1:
                        data[:, :, :, ri: ri+2] = 0.
                        log.debug("bad interleaved reference at pixels {} {}"
                                  .format(ri, ri+2))
    log.debug("total bad reference pixels = {}".format(sum(total_bad_pixels)))


def decode_mask(output, mask):
    """Interpret the MASK column of the DQ table.

    As per the ESA CDP3 document:
    "There is also a DQ extension that holds a binary table with three
    columns (OUTPUT, ODD_EVEN, and MASK) and eight rows. In the current
    IRS2 implementation, one jumps 32 times to odd and 32 times to even
    reference pixels, which are then read twice consecutively. Therefore,
    the masks are 32 bit unsigned integers that encode bad interleaved
    reference pixels/columns from left to right (increasing column index)
    in the native detector frame. When a bit is set, the corresponding
    reference data should not be used for the correction."

    Parameters
    ----------
    output : int
        An amplifier output number, 1, 2, 3, or 4.

    mask : uint32
        A mask value.

    Returns
    -------
    bits : list
        A list of the indices of bits set in the `mask` value.
    """

    # The bit number corresponds to a count of groups of reads of the
    # interleaved reference pixels. The 32-bit unsigned integer encoding
    # has increasing index, following the amplifier readout direction.

    flags = np.array([2**n for n in range(32)], dtype=np.uint32)
    temp = np.bitwise_and(flags, mask)
    bits = np.where(temp > 0)[0]
    bits = list(bits)
    if output // 2 * 2 == output:
        # account for the readout orientation of even amplifiers
        bits = [31 - bit for bit in bits]
    # order indeces increasing from left to right always
    bits.sort()

    return bits


def shift(arr, num):
    """ Function to reproduce the IDL shift function. It preallocates an empty
    array and assigns the slice.

    Parameters
    ----------
    arr: numpy array
        The 1-D array that we wish to shift

    num: integer
        The number of places we want to shift the array, positive shifts to the right
        and negative shifts to the left. Both shifts are circular, e.g.
        a = np.arange(5)
        shift(a, 3) returns array([2., 3., 4., 0., 1.])
        shift(a, -3) returns array([3., 4., 0., 1., 2.])

    Returns
    -------
    result: np.array
        Shifted 1-D array of same length

    """
    result = arr.copy()
    if num > 0:
        result[num:] = arr[:-num]
        result[:num] = arr[-num:]
    else:
        result[num:] = arr[:-num]
        result[:num] = arr[-num:]
    return result


def ols_line(x, y):
    """Fit a straight line using ordinary least squares."""

    xf = x.ravel()
    yf = y.ravel()
    if len(xf) < 1 or len(yf) < 1:
        return 0., 0.

    groups = float(len(xf))
    sum_x2 = (xf**2).sum()
    sum_xy = (xf * yf).sum()
    sum_xf = xf.sum()
    sum_yf = yf.sum()

    slope = (groups * sum_xy - sum_xf * sum_yf) / (groups * sum_x2 - sum_xf**2)
    intercept = (sum_yf - slope * sum_xf) / groups

    return intercept, slope


def remove_slopes(data0, ngroups, ny, row):

    # Fitting and removal of slopes per frame to remove issues at frame boundaries.
    # IDL:  time = findgen(row, s[2])

    time_arr = np.arange(ny * row, dtype=np.float32).reshape((ny, row))
    row4plus4 = np.array([0, 1, 2, 3, 2044, 2045, 2046, 2047])
    time_arr -= np.mean(time_arr).round(decimals=0)
    print('np.shape(time_arr) = ', np.shape(time_arr), time_arr)

    # For ab_3, it is OK to use the same index order as the IDL code.
    ab_3 = np.zeros((2, ngroups, 5), dtype=np.float32)   # for top+bottom ref pixel rows
    for i in range(5):
        for k in range(ngroups):
            # mask is 2-D, since both row4plus4 and : have more than one element.
            mask = np.where(data0[i, k, row4plus4, :] != 0.)
            (intercept, slope) = ols_line(time_arr[row4plus4, :][mask],
                                          data0[i, k, row4plus4, :][mask])
            ab_3[0, k, i] = intercept
            ab_3[1, k, i] = slope

    for i in range(5):
        for k in range(ngroups):
            # weight is 0 where data0 is 0, else 1.
            weight = np.where(data0[i, k, :, :] != 0., 1., 0.)
            data0[i, k, :, :] -= (time_arr * ab_3[1, k, i] + ab_3[0, k, i]) * weight


def replace_bad_pixels(data0, ngroups, ny, row):

    # Use cosine weighted interpolation to replace 0.0 values and bad
    # pixels and gaps. (initial guess)

    # s[1] = nx  s[2] = ny  s[3] = ngroups
    w_ind = np.arange(1, 32, dtype=np.float32) / 32.
    w = np.sin(w_ind * np.pi)
    print('np.shape(w) = ', np.shape(w), w)
    for kk in range(5):
        for jj in range(ngroups):
            dat = data0[kk, jj, :, :].reshape(row * ny)
            mask = np.where(dat != 0., 1., 0.)
            numerator = convolve1d(dat, w, mode='wrap')
            denominator = convolve1d(mask, w, mode='wrap')
            div_zero = denominator == 0.          # check for divide by zero
            numerator = np.where(div_zero, 0., numerator)
            denominator = np.where(div_zero, 1., denominator)
            dat = numerator / denominator
            dat = dat.reshape(ny, row)
            mask = mask.reshape(ny, row)
            data0[kk, jj] += dat * (1. - mask)
    print('dat[:3, :3] = ', dat[:3, :3])
    print('data0[1, 20, :3, :3] = ', data0[1, 20, :3, :3])


def subtract_reference(data0, alpha, beta, irs2_mask, scipix_n, refpix_r, pad):
    """Subtract reference output and pixels for the current integration. This routine
    is based off of the NASA Goddard IDL routines part_a.pro, part_b.pro, and part_c.pro
    obtained from https://jwst.nasa.gov/content/forScientists/publications.html -
    Improved Reference Sampling and Subtraction: A Technique for Reducing the Read
    Noise of Near-infrared Detector Systems, Journal Article, by B.J.
    Rauscher et al., April 2017

    Parameters
    ----------
    data0: ramp data
        The science data for the current integration.  The shape is
        expected to be (ngroups, ny, 3200), where ngroups is the number of
        groups, and ny is the pixel height of the image.  The width 3200
        of the image includes the "normal" pixel data, plus the embedded
        reference pixels, and the reference output.

    alpha: ndarray
        This is a 2-D array of values read from the reference file.  The
        first axis is the sector number (but only for the normal pixel
        data and reference pixels, not the reference output).  The second
        axis has length 2048 * 712, corresponding to the time-ordered
        arrangement of the data.  For each sector, the correction is
        applied as follows:  data * alpha[i] + reference_output * beta[i].

    beta: ndarray
        Data read from the reference file.  See `alpha` for details.

    irs2_mask: Boolean, 1-D array
        True means the element corresponds to a normal pixel in the raw,
        IRS2-format data.  False corresponds either to a reference output
        pixel or to one of the interspersed reference pixel values.

    scipix_n: int
        Number of regular samples before stepping out to collect
        reference samples.

    refpix_r: int
        Number of reference samples before stepping back in to collect
        regular samples.

    pad: int
        The effective number of pixels sampled during the pause at the end
        of each row (new-row overhead).

    Returns
    -------
    data0: ramp data
        The science data for the current integration, with reference output
        and embedded reference pixels subtracted and also removed, leaving
        only the normal pixel data (including the reference pixels on each
        edge).  The shape is expected to be (ngroups, ny, nx), where
        nx = ny = 2048.
    """

    # The gain is not dealt with at this point, it will be taken care of
    # at the gain_scale step of the pipeline.

    # This is the effective number of pixels sampled during the
    # pause at the end of each row. The padding is needed to preserve phase
    # of temporally periodic signals.
    # See expression in equation 1 in IRS2_Handoff.pdf.
    # row = 712, if scipix_n = 16, refpix_r = 4, pad = 8.
    row = (scipix_n + refpix_r + 2) * 512 // scipix_n + pad

    # IDL definitions:  s = size(data0)
    # If data0 is the data for one integration, then in IDL:
    # s[0] would be 3
    # s[1] = shape[2] = nx, the length of the X axis
    # s[2] = shape[1] = ny, the length of the Y axis
    # s[3] = shape[0] = ngroups, the number of groups (or frames)
    shape = data0.shape
    ngroups = shape[0]
    ny = shape[1]
    nx = shape[2]

    nn = np.arange(ngroups, dtype=float)
    # subtract mean(nn) to minimize correlation of slope and offset
    nn -= np.mean(nn)

    # hnorm is an array of column indices of normal pixels.
    # len(hnorm) = 512; len(href) = 128
    # len(hnorm1) = 512; len(href1) = 128
    ind_n = np.arange(512, dtype=np.intp)
    hnorm = ind_n + refpix_r * ((ind_n + scipix_n // 2) // scipix_n)

    # href is an array of column indices of reference pixels.
    ind_ref = np.arange(512 // scipix_n * refpix_r, dtype=np.intp)
    href = ind_ref + scipix_n * (ind_ref // refpix_r) + scipix_n // 2

    hnorm1 = ind_n + (refpix_r + 2) * ((ind_n + scipix_n // 2) // scipix_n)
    href1 = ind_ref + (scipix_n + 2) * (ind_ref // refpix_r) + scipix_n // 2 + 1

    # remove linear trends per pixel
    data0 = data0.reshape(ngroups, ny * nx)
    print('np.shape(data0) = ', np.shape(data0))
    # the equivalent of IDL nn##data0 in python is matrix multiplication np.matmul
    sxy = np.matmul(nn, data0)
    print('np.shape(sxy) = ', np.shape(sxy))
    sx = sum(nn)
    sxx = sum(nn**2)
    sy = np.matmul(nn*0 + 1, data0)
    print('np.shape(sy) = ', np.shape(sy))
    b = (sy * sxx - sxy * sx) / (ngroups * sxx - sx**2)
    print('np.shape(b) = ', np.shape(b))
    # subtract the offset only
    for i in range(ngroups):
        data0[i] -= b
    sxy, sy, b = 0, 0, 0
    mtxprod = np.matmul( nn*0 + 1, data0**2 )
    sig_data = np.sqrt( mtxprod / (ngroups - 1))
    print('np.shape(sig_data) = ', np.shape(sig_data))
    sig_data = sig_data.reshape((ny, nx))
    data0 = data0.reshape((ngroups, ny, nx))
    print('np.shape(data0) = ', np.shape(data0))
    print('np.shape(sig_data) = ', np.shape(sig_data))

    nx5 = nx // 5
    print('nx5 = ', nx5)
    mask = np.zeros((ny, nx), dtype=bool)
    print('np.shape(mask) = ', np.shape(mask))
    # set reference output to no masking
    mask[:, :nx5] = True
    print('np.shape(mask[:, :nx5]) = ', np.shape(mask[:, :nx5]))
    nmask = mask.copy()
    for i in range(1, 5):
        # the IDL indices used are 4:2043 because it includes 2043
        # python slicing does not include 2044 so it ends at 2043
        nmask[4: 2044, hnorm + i*nx5 ] = True
    print('np.shape(nmask[4: 2044, hnorm + 1*nx5 ]) = ', np.shape(nmask[4: 2044, hnorm + 1*nx5 ]))
    nmask[:, nx5: nx5+4] = False
    print('np.shape(nmask[:, nx5: nx5+4]) = ', np.shape(nmask[:, nx5: nx5+4]))
    nmask[:, nx-4:] = False
    print('np.shape(nmask[:, nx-4:]) = ', np.shape(nmask[:, nx-4:]))
    for amp in range(1, 5):   # skip the first since it requires no mask
        print('amp = ', amp, amp * nx5, (amp+1) * nx5)
        msk = np.zeros((ny, nx5), dtype=bool)
        nmask_slice = nmask[:, amp * nx5: (amp+1) * nx5]
        hhnorm = np.where(nmask_slice == True)
        print('hhnorm = ', np.shape(hhnorm))
        hhref = np.where(nmask_slice == False)
        print('hhref = ', np.shape(hhref))
        # there is no equivalent to IDL HISTOGAUSS in python so do separately
        # IDL shows a plot, this is skipped in the python version
        dat1 = sig_data[:, amp * nx5: (amp+1) * nx5]
        dat2 = sig_data[:, :nx5]
        dat = dat1 - dat2
        print('np.shape(dat1), np.shape(dat2), np.shape(dat) = ', np.shape(dat1), np.shape(dat2), np.shape(dat))
        datnorm = dat[hhnorm]
        print('np.shape(datnorm) = ', np.shape(datnorm))
        print('np.mean(datnorm), np.std(datnorm) = ', np.mean(datnorm), np.std(datnorm))
        
        # this is what is stored in IDL histogauss variable bb
        #   bb = coefficients of the Gaussian fit: Height, mean, sigma
        #      bb[0]= the height of the Gaussian
        #      bb[1]= the mean
        #      bb[2]= the standard deviation
        #      bb[3]= the half-width of the 95% conf. interval of the standard
        #            mean
        #      bb[4]= 1/(N-1)*total( (y-mean)/sigma)^2 ) = a measure of normality
        (mu, sigma) = norm.fit(datnorm)   # best fit of data
        print('hhnorm gaussian fit, mu, sigma =', mu, sigma)
        msk[hhnorm] = np.where(np.abs(datnorm - mu) < 4 * sigma, True, False)
        print('np.shape(msk[hhnorm] ) = ', np.shape(msk[hhnorm]), msk[hhnorm] )
        datref = dat[hhref]
        print('np.shape(datref) = ', np.shape(datref))
        (mu, sigma) = norm.fit(datref)   # best fit of data
        print('hhref gaussian fit, mu, sigma = ', mu, sigma)
        msk[hhref] = np.where(np.abs(datref - mu) < 4 * sigma, True, False)
        print('np.shape(msk[hhref] ) = ', np.shape(msk[hhref]), msk[hhref] )
        mask[:, amp * nx5: (amp+1) * nx5] = msk

    # expand the mask to cover neighbor pixels in x and y
    print('np.shape(mask) = ', np.shape(mask))
    mask[1] *= shift(mask[1], -1) * shift(mask[1], 1)
    mask[0] *= shift(mask[0], -1) * shift(mask[0], 1)
    # reset the reference output to no masking
    mask[:, :nx5] = True

    # apply outlier mask and reshape for later use
    print('np.shape(data0) = ', np.shape(data0))
    print('applying mask')
    for ng in range(ngroups):
        data0[ng] *= mask
    print('mask applied! ')
    mask0 = mask.copy()
    print('np.shape(mask0) = ', np.shape(mask0))
    mask0 = mask0.reshape(ny, 5, nx5)
    print('np.shape(mask0) = ', np.shape(mask0))
    # IDL mask0 order at this point is       (nx5, 5, ny)  -> indices 0, 1, 2
    # numpy outlimask order at this point is ( ny, 5, nx5) -> indices 0, 1, 2
    # IDL transpose mask0 to order    (nx5, ny, 5)  -> indices 0, 2, 1
    # this corresponds to numpy order ( 5, ny, nx5)  -> indices 1, 0, 2
    mask0 = np.transpose(mask0, (1, 0, 2))
    # revert order to match data for readout direction
    mask0[0, :, :] = mask0[0, :, ::-1]
    mask0[2, :, :] = mask0[2, :, ::-1]
    mask0[4, :, :] = mask0[4, :, ::-1]
    print('transposed and reversed np.shape(mask0) = ', np.shape(mask0))

    # IDL:  data0 = reform(data0, s[1]/5, 5, s[2], s[3], /over)
    #                             nx/5,   5, ny,   ngroups    (IDL)
    print('np.shape(data0) = ', np.shape(data0))
    data0 = data0.reshape((ngroups, ny, 5, nx // 5))
    print('reshaped np.shape(data0) = ', np.shape(data0))

    # current order:  nx/5, 5, ny, ngroups    (IDL)
    # current order:  ngroups, ny, 5, nx/5    (numpy)
    #                 0        1   2  3       current numpy indices
    # transpose to:   nx/5, ny, ngroups, 5    (IDL)
    # transpose to:   5, ngroups, ny, nx/5    (numpy)
    #                 2  0        1   3       transpose order for numpy
    # Therefore:      0 1 2 3  -->  2 0 1 3   transpose order for numpy
    # Here is another way to look at it:
    # IDL:    0 1 2 3  -->  0 2 3 1
    #         3 2 1 0       1 3 2 0 (IDL indices, but reversed to numpy order)
    # numpy:  0 1 2 3  -->  2 0 1 3
    # IDL:  data0 = transpose(data0, [0,2,3,1])
    data0 = np.transpose(data0, (2, 0, 1, 3))

    # Flip the direction of the X axis for every other output, so the readout
    # direction in data0 will be the same for every output.
    data0[0, :, :, :] = data0[0, :, :, ::-1]
    data0[2, :, :, :] = data0[2, :, :, ::-1]
    data0[4, :, :, :] = data0[4, :, :, ::-1]
    print('transposed and reversed np.shape(data0) = ', np.shape(data0))

    # convert to time sequences of normal pixels and reference pixels.
    # IDL (lines 122 - 126 in part_A.pro):
    # d0 = fltarr(s[1] / 5 + pad + 2 * (512 / scipix_n), s[2], s[3], 5)
    # Note:  nx // 5 + pad + 2 * (512 // scipix_n) = 640 + 8 +64 = 712.
    # hnorm1[-1] = 703, and hnorm[-1] = 639, so 703 - 639 = 64.
    # 8 is the pad value.

    d0 = np.zeros((5, ngroups, ny, row))   # (5, ngroups, 2048, 712)
    print('np.shape(d0) = ', np.shape(d0))
    # IDL:  d0[hnorm1,*,*,*] = data0[hnorm,*,*,*]
    # IDL:  d0[href1,*,*,*] = data0[href,*,*,*]
    # IDL:  data0 = temporary(d0)
    d0[:, :, :, hnorm1] = data0[:, :, :, hnorm]
    d0[:, :, :, href1] = data0[:, :, :, href]
    del data0
    data0 = d0.copy()
    del d0

    # Fitting and removal of slopes per frame to remove issues at frame boundaries
    remove_slopes(data0, ngroups, ny, row)
    print('removed slopes per frame')

    # Use cosine weighted interpolation to replace 0.0 values and bad
    # pixels and gaps. (initial guess)
    replace_bad_pixels(data0, ngroups, ny, row)
    print('finished cosine weighted interpolation')

    # Fill in bad pixels, gaps, and reference data locations in the normal
    # data, using Fourier filtering/interpolation
    print('entering fill_bad_regions')
    fourier_filter_replace(data0, ngroups, ny, nx, row, scipix_n, refpix_r, pad,
        hnorm, hnorm1, mask0)
    print('done with fill_bad_regions')


def fourier_filter_replace(data0, ngroups, ny, nx, row, scipix_n, refpix_r, pad,
        hnorm, hnorm1, mask0):

    # Use Fourier filter/interpolation to replace
    # (a) bad pixel, gaps, and reference data in the time-ordered normal data
    # (b) gaps and normal data in the time-ordered reference data
    # This "improves" upon the cosine interpolation performed above.

    # Parameters for the filter to be used.
    # length of apodization cosine filter
    elen = 110000 // (scipix_n + refpix_r + 2)   # elen = 5000

    # max unfiltered frequency
    blen = (512 + 512 // scipix_n * (refpix_r + 2) + pad) // \
           (scipix_n + refpix_r + 2) * ny // 2 - elen // 2   # blen = 30268

    # Construct the filter [1, cos, 0, cos, 1].
    temp_a1 = (np.cos(np.arange(elen, dtype=np.float32) *
                      np.pi / float(elen)) + 1.) / 2.
    temp_a2 = np.concatenate((np.ones(blen, dtype=np.float32),
                              temp_a1,
                              np.zeros(row * ny // 2 - 2 * blen - 2 * elen,
                                       dtype=np.float32),
                              temp_a1[::-1].copy(),
                              np.ones(blen, dtype=np.float32)))
    aa = np.matmul(np.transpose(temp_a2.copy()), np.ones(ngroups, dtype=np.float32))
    afilter = temp_a2.copy()
    print('np.shape(aa) = ', np.shape(aa), aa[:3])
    print('np.shape(afilter) = ', np.shape(afilter), afilter[:3])
    exit()





    # Setup various lists of indices that will be used in subsequent
    # sections for keeping/shuffling reference pixels in various arrays
    #
    # The comments are for scipix_n = 16, refpix_r = 4
    n0 = 512 // scipix_n
    n1 = scipix_n + refpix_r + 2
    ht = np.arange(n0 * n1, dtype=np.int32).reshape((n0, n1))   # (32, 22)
    ht[:, 0:(scipix_n - refpix_r) // 2 + 1] = -1
    ht[:, scipix_n // 2 + 1 + 3 * refpix_r // 2:] = -1
    hs = ht.copy()
    # ht is like href1, but extended over gaps and first and last norm pix
    mask = (ht >= 0)
    ht = ht[mask]               # 1-D, length = 2 * refpix_r * 512 / scipix_n

    # IDL:  hs[scipix_n/2 + 1-refpix_r/2:scipix_n/2 + refpix_r + refpix_r/2,*] =
    #       hs[reform([transpose(reform(indgen(refpix_r),refpix_r/2,2)),
    #           transpose(reform(indgen(refpix_r),refpix_r/2,2))],refpix_r * 2)
    #           + scipix_n/2 + 1,*]  ; WIRED for R=2^(int)

    indr = np.arange(refpix_r, dtype=np.intp).reshape((2, refpix_r // 2))
    # indr_t =
    # [[0 2]
    #  [1 3]]
    indr_t = indr.transpose()

    # Before flattening, two_indr_t =
    # [[0 2 0 2]
    #  [1 3 1 3]]
    # After flattening, two_indr_t = [0 2 0 2 1 3 1 3].
    two_indr_t = np.concatenate((indr_t, indr_t), axis=1).flatten()
    two_indr_t += (scipix_n // 2 + 1)     # [9 11 9 11 10 12 10 12]
    hs[:, scipix_n // 2 + 1 - refpix_r // 2:
       scipix_n // 2 + 1 + refpix_r // 2 + refpix_r] = hs[:, two_indr_t]
    mask = (hs >= 0)
    hs = hs[mask]    # hs is now 1-D

    if refpix_r % 4 == 2:
        len_hs = len(hs)
        temp_hs = hs.reshape(len_hs // 2, 2)
        temp_hs = temp_hs[:, ::-1]
        hs = temp_hs.flatten()

    # Construct the reference data: this is done in a big loop over the
    # four "sectors" of data in the image, corresponding to the amp regions.
    # Data from each sector is operated on independently and ultimately
    # the corrections are subtracted from each sector independently.
    shape_d = data0.shape
    for k in range(1, 5):
        log.debug(f'processing sector {k}')
        print(f'processing sector {k}')

        # At this point in the processing data0 has shape (5, ngroups, 2048, 712),
        # assuming normal IRS2 readout settings. r0k contains a subset of the
        # data from 1 sector of data0, with shape (ngroups, 2048, 256)
        r0k = np.zeros((shape_d[1], shape_d[2], shape_d[3]), dtype=np.float32)
        temp = data0[k, :, :, hs].copy()
        temp = np.transpose(temp, (1, 2, 0))
        r0k[:, :, ht] = temp
        del temp

        # data0 has shape (5, ngroups, ny, row).  See the section above where
        # d0 was created, then copied (moved) to data0.
        # sd[1] = shape_d[3]   row (712)
        # sd[2] = shape_d[2]   ny (2048)
        # sd[3] = shape_d[1]   ngroups
        # sd[4] = shape_d[0]   5
        # s is used below, so for convenience, here are the values again:
        # s[1] = shape[2] = nx
        # s[2] = shape[1] = ny
        # s[3] = shape[0] = ngroups

        # IDL and numpy differ in where they apply the normalization for the
        # FFT.  This really shouldn't matter.
        normalization = float(shape_d[2] * shape_d[3])

        if beta is not None:
            # IDL:  refout0 = reform(data0[*,*,*,0], sd[1] * sd[2], sd[3])
            refout0 = data0[0, :, :, :].reshape((shape_d[1], shape_d[2] * shape_d[3]))

            # IDL:  refout0 = fft(refout0, dim=1, /over)
            # Divide by the length of the axis to be consistent with IDL.
            refout0 = np.fft.fft(refout0, axis=1) / normalization

        # IDL:  r0 = reform(r0, sd[1] * sd[2], sd[3], 5, /over)
        r0k = r0k.reshape((shape_d[1], shape_d[2] * shape_d[3]))
        r0k = r0k.astype(np.complex64)
        r0k_fft = np.fft.fft(r0k, axis=1) / normalization

        # Note that where the IDL code uses alpha, we use beta, and vice versa.
        # IDL:  for k=0,3 do oBridge[k]->Execute,
        #           "for i=0, s3-1 do r0[*,i] *= alpha"
        r0k_fft *= beta[k - 1]

        # IDL:  for k=0,3 do oBridge[k]->Execute,
        #           "for i=0, s3-1 do r0[*,i] += beta * refout0[*,i]"
        if beta is not None:
            r0k_fft += (alpha[k - 1] * refout0)
        del refout0

        # IDL:  for k=0,3 do oBridge[k]->Execute,
        #           "r0 = fft(r0, 1, dim=1, /overwrite)", /nowait
        r0k = np.fft.ifft(r0k_fft, axis=1) * normalization
        del r0k_fft

        # sd[1] = shape_d[3]   row (712)
        # sd[2] = shape_d[2]   ny (2048)
        # sd[3] = shape_d[1]   ngroups
        # sd[4] = shape_d[0]   5
        # IDL:  r0 = reform(r0, sd[1], sd[2], sd[3], 5, /over)
        r0k = r0k.reshape(shape_d[1], shape_d[2], shape_d[3])
        r0k = r0k.real
        r0k = r0k[:, :, hnorm1]

        # Subtract the correction from the data in this sector
        data0[k, :, :, hnorm1] -= np.transpose(r0k, (2, 0, 1))
        del r0k

    # End of loop over 4 sectors

    # Original data0 array has shape (5, ngroups, 2048, 712). Now that
    # correction has been applied, remove the interleaved reference pixels.
    # This leaves data0 with shape (5, ngroups, 2048, 512).
    data0 = data0[:, :, :, hnorm1]
    print('after sector loop, np.shape(data0), np.shape(hnorm1) = ', np.shape(data0), np.shape(hnorm1))

    # Unflip the data in the sectors that have opposite readout direction
    data0[2, :, :, :] = data0[2, :, :, ::-1]
    data0[4, :, :, :] = data0[4, :, :, ::-1]

    # IDL:  data0 = transpose(data0, [0,3,1,2])  0, 1, 2, 3 --> 0, 3, 1, 2
    # current order:  512, ny, ngroups, 5     (IDL)
    # current order:  5, ngroups, ny, 512     (numpy)
    #                 0  1        2   3       current numpy indices
    # transpose to:   512, 5, ny, ngroups     (IDL)
    # transpose to:   ngroups, ny, 5, 512     (numpy)
    #                 1        2   0  3       transpose order for numpy
    # Therefore:      0 1 2 3  -->  1 2 0 3   transpose order for numpy
    # After transposing, data0 will have shape (ngroups, 2048, 5, 512).
    data0 = np.transpose(data0, (1, 2, 0, 3))

    # Reshape data0 back to its normal (ngroups, 2048, 2048), which has
    # the interleaved reference pixels stripped out.
    # IDL:  data0 = reform(data0[*, 1:*, *, *], s[2], s[2], s[3], /over)
    # Note:  ny x ny, not ny x nx.
    data0 = data0[:, :, 1:, :].reshape((ngroups, ny, ny))

    # b_offset is the average over the ramp that we subtracted near the
    # beginning; add it back in.
    # Shape of b_offset is (2048, 3200), but data0 is (ngroups, 2048, 2048),
    # so a mask is applied to b_offset to remove the reference pix locations.
    data0 += b_offset[..., irs2_mask]

    return data0
