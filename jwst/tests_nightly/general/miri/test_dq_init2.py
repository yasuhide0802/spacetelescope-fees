import os
import pytest
from astropy.io import fits as pf
from jwst.dq_init.dq_init_step import DQInitStep

from ..helpers import add_suffix

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]

def test_dq_init_miri2():
    """

    Regression test of dq_init step performed on uncalibrated MIRI data.

    """
    output_file_base, output_file = add_suffix('dqinit2_output.fits', 'dq_init')

    try:
        os.remove(output_file)
    except:
        pass

    DQInitStep.call(_bigdata+'/miri/test_dq_init/jw80600012001_02101_00003_mirimage_uncal.fits',
                    output_file=output_file_base
    )

    h = pf.open(output_file)
    href = pf.open(_bigdata+'/miri/test_dq_init/jw80600012001_02101_00003_mirimage_dqinit.fits')
    newh = pf.HDUList([h['primary'],h['sci'],h['err'],h['pixeldq'],h['groupdq']])
    newhref = pf.HDUList([href['primary'],href['sci'],href['err'],href['pixeldq'],href['groupdq']])
    result = pf.diff.FITSDiff(newh,
                              newhref,
                              ignore_keywords = ['DATE','CAL_VER','CAL_VCS','CRDS_VER','CRDS_CTX'],
                              rtol = 0.00001
    )
    result.report()
    try:
        assert result.identical == True
    except AssertionError as e:
        raise AssertionError(e)
