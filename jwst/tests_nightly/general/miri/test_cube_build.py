import os
from astropy.io import fits as pf
from jwst.cube_build.cube_build_step import CubeBuildStep
from jwst import datamodels

BIGDATA = os.environ['TEST_BIGDATA']

def test_cubebuild_miri():
    """

    Regression test of cube_build performed on MIRI MRS data.

    """
    try:
        os.remove("cubebuild1_output.fits")
    except:
        pass

    input_model = datamodels.IFUImageModel(BIGDATA+'/miri/test_cube_build/jw10001001001_01101_00001_mirifushort_cal.fits')
    cube_model = CubeBuildStep.call(input_model, output_type='multi')
    cube_model.save('cubebuild1_output.fits')

    h = pf.open('jw10001001001_01101_00001_mirifushort_s3d.fits')
    href = pf.open(BIGDATA+'/miri/test_cube_build/jw10001001001_01101_00001_mirifushort_s3d_ref.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['dq'],h['wmap']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['dq'],href['wmap']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.000001
    )
    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        print(result.report())
        raise AssertionError(e)

