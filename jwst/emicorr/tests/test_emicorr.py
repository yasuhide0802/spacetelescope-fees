"""

Unit tests for EMI correction

"""


import numpy as np
from jwst.emicorr import emicorr
from stdatamodels.jwst.datamodels import Level1bModel

default_subarray_cases = emicorr.default_subarray_cases

def mk_emi_mdl(data, subarray, readpatt, detector):
    # create input_model
    input_model = Level1bModel(data=data)
    input_model.meta.instrument.detector = detector
    input_model.meta.exposure.type = 'MIR_4QPM'
    input_model.meta.subarray.name = subarray
    input_model.meta.exposure.readpatt = readpatt
    input_model.meta.subarray.xsize = 288
    input_model.meta.subarray.xstart = 1
    input_model.meta.exposure.nsamples = 1
    return input_model


def test_get_subarcase():
    subarray, readpatt, detector = 'MASK1550', 'FASTR1', 'MIRIMAGE'
    subarray_info_r = emicorr.get_subarcase(default_subarray_cases, subarray, readpatt, detector)
    subname_r, rowclocks_r, frameclocks_r, freqs2correct_r = subarray_info_r

    # set up a fake configuration
    subarray_info_f = emicorr.get_subarcase(default_subarray_cases, 'FULL', 'kkkkk', detector)
    subname_f, rowclocks_f, frameclocks_f, freqs2correct_f = subarray_info_f

    # test if we get the right configuration
    compare_real = ['MASK1550', 82, 23968, ["Hz390", "Hz10"]]
    subname_real, rowclocks_real, frameclocks_real, freqs2correct_real = compare_real

    assert subname_real == subname_r
    assert rowclocks_real == rowclocks_r
    assert frameclocks_real == frameclocks_r
    assert freqs2correct_real[0] == freqs2correct_r[0]


def test_apply_emicorr():
    data = np.ones((1, 5, 20, 20))
    input_model = mk_emi_mdl(data, 'MASK1550', 'FASTR1', 'MIRIMAGE')
    emicorr_model = None
    save_onthefly_reffile = None
    output_model = emicorr.apply_emicorr(input_model, emicorr_model, save_onthefly_reffile)

    compare_model = Level1bModel(data=data)

    # correction should equal input for data of all 1.0
    assert compare_model.data.all() == output_model.data.all()


def test_sloper():
    data = np.ones((5, 5, 5))
    outarray, intercept = emicorr.sloper(data)

    compare_arr, compare_intercept = np.zeros((5, 5)), np.ones((5, 5))

    assert compare_arr.all() == outarray.all()
    assert compare_intercept.all() == intercept.all()


def test_minmed():
    data = np.ones((5, 5, 5))
    compare_arr = data.copy()

    medimg = emicorr.minmed(data)

    assert compare_arr.all() == medimg.all()


def test_iter_stat_sig_clip():
    data = np.ones((5, 5, 5, 5))
    data[0, 0, 0, 0] = 0.55
    data[1, 1, 1, 1] = 1.11
    data[2, 2, 2, 2] = 2.22
    data[3, 3, 3, 3] = 3.33
    data[4, 4, 4, 4] = 4.44
    u = np.where(data > 1)
    dmean, dmedian, dsigma, dmask = emicorr.iter_stat_sig_clip(data[u])

    compare_mean, compare_median = 2.7750000000000004, 0.0
    compare_sigma, compare_mask = 1.326703756361177, np.ones(np.size(data), dtype='b')+1

    assert compare_mean == dmean
    assert compare_median == dmedian
    assert compare_sigma == dsigma
    assert compare_mask.all() == dmask.all()


def test_rebin():
    data = np.ones(10)
    data[1] = 0.55
    data[5] = 1.55
    data[9] = 2.0
    rebinned_data = emicorr.rebin(data, [7])

    compare_arr = np.array([1., 0.55, 1. , 1., 1.55, 1., 1.])

    assert compare_arr.all() == rebinned_data.all()
