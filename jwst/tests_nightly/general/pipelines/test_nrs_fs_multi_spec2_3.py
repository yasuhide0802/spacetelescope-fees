import os
from astropy.io import fits as pf
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

BIGDATA = os.environ['TEST_BIGDATA']

def test_nrs_fs_multi_spec2_3():
    """

    Regression test of calwebb_spec2 pipeline performed on NIRSpec fixed-slit data
    using the ALLSLITS subarray and detector NRS2.

    """

    Spec2Pipeline.call(BIGDATA+'/pipelines/jw84600002001_02101_00001_nrs2_rate.fits',
                       config_file='calwebb_spec2.cfg')

    na = 'jw84600002001_02101_00001_nrs2_cal.fits'
    nb = BIGDATA+'/pipelines/jw84600002001_02101_00001_nrs2_cal_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci',1],h['err',1],h['dq',1],h['relsens',1],h['wavelength',1],
                                    h['pathloss_pointsource',1],h['wavelength_pointsource',1],
                                    h['pathloss_uniformsource',1],h['wavelength_uniformsource',1]])
    newhref = pf.HDUList([href['primary'],href['sci',1],href['err',1],href['dq',1],href['relsens',1],href['wavelength',1],
                                          href['pathloss_pointsource',1],href['wavelength_pointsource',1],
                                          href['pathloss_uniformsource',1],href['wavelength_uniformsource',1]])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb) 

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    na = 'jw84600002001_02101_00001_nrs2_s2d.fits'
    nb = BIGDATA+'/pipelines/jw84600002001_02101_00001_nrs2_s2d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['sci',1],h['wht',1],h['con',1],h['relsens',1]])
    newhref = pf.HDUList([href['primary'],href['sci',1],href['wht',1],href['con',1],href['relsens',1]])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb) 

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

    na = 'jw84600002001_02101_00001_nrs2_x1d.fits'
    nb = BIGDATA+'/pipelines/jw84600002001_02101_00001_nrs2_x1d_ref.fits'
    h = pf.open(na)
    href = pf.open(nb)
    newh = pf.HDUList([h['primary'],h['extract1d',1]])
    newhref = pf.HDUList([href['primary'],href['extract1d',1]])
    result = pf.diff.FITSDiff(newh, newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001)

    print (' Fitsdiff comparison between product file - a:', na)
    print (' ... and the reference file - b:', nb) 

    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

