"""Test suite for ensuring correct FGS pointing"""
import os.path
import logging
from pathlib import Path
import pytest

from astropy.time import Time
from numpy import array
from numpy import isclose

from jwst.datamodels import Level1bModel
from jwst.lib import engdb_mast
from jwst.lib import engdb_tools
from jwst.lib import set_telescope_pointing as stp
from jwst.lib.tests.engdb_mock import EngDB_Mocker

# Set logging for the module to be tested.
logger = logging.getLogger(stp.__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.addHandler(handler)

# Get database paths.
DATA_PATH = Path(__file__).parent / 'data'
db_1029 = DATA_PATH / 'engdb_jw01029'
siaf_db = DATA_PATH / 'siaf.db'

# Time frame
OBSSTART = '2022-05-23T00:36:08.000'
OBSEND = '2022-05-23T00:36:22.480'

# Define minimal model meta structure
WCS_META = {
    'meta': {
        'aperture': {
            'name': 'FGS1_FULL',
        },
        'exposure': {
            'type': 'FGS_ACQ1',
        },
        'guidestar': {
            'gs_ra': 45.1234,
            'gs_dec': -45.1234,
        },
        'instrument': {
            'detector': 'GUIDER1',
        },
        'observation': {
            'date_beg': OBSSTART,
            'date_end': OBSEND,
            'date': '2017-01-01',
        },
    }
}


@pytest.mark.parametrize('multi_fixture', ['engdb_jw01029', 'mast'], indirect=True)
@pytest.mark.parametrize('exp_type, expected',
                         [('fgs_acq1', stp.GuideStarPosition(position=(17.7885990143, 42.6522407532), corner=(1193, 343), size=(128, 128))),
                          ('fgs_acq2', stp.GuideStarPosition(position=(17.8298149109, 42.65200042725), corner=(1268, 386), size=(32, 32)))])
def test_gs_position_acq(multi_fixture, exp_type, expected):
    """Test the mnemonics reading"""
    engdb = multi_fixture

    # Perform operation
    mnemonics = stp.get_mnemonics(Time(OBSSTART).mjd, Time(OBSEND).mjd, 60.,
                                  stp.FGS_ACQ_MNEMONICS, engdb_url=engdb.base_url)
    gs_position = stp.gs_position_acq(stp.FGS_ACQ_MNEMONICS, mnemonics, exp_type)

    # Test
    assert gs_position == expected


@pytest.mark.parametrize('multi_fixture', ['engdb_jw01029', 'mast'], indirect=True)
def test_fgs_pointing(multi_fixture):
    engdb = multi_fixture

    # Setup model
    model = make_level1b()

    # Update wcs
    stp.update_wcs(model, engdb_url=engdb.base_url)

    # Test results
    assert isclose(model.meta.wcsinfo.pc1_1, -0.9997617158628777, atol=1e-15)
    assert isclose(model.meta.wcsinfo.pc1_2, -0.02166140686177685, atol=1e-15)
    assert isclose(model.meta.wcsinfo.pc2_1, -0.02166140686177685, atol=1e-15)
    assert isclose(model.meta.wcsinfo.pc2_2, 0.9997653641994049, atol=1e-15)

    assert isclose(model.meta.wcsinfo.crpix1, 36.23964226769749, atol=1e-15)
    assert isclose(model.meta.wcsinfo.crpix2, 68.3690778810028, atol=1e-15)
    assert isclose(model.meta.wcsinfo.crval1, 45.1234, atol=1e-15)
    assert isclose(model.meta.wcsinfo.crval2, -45.1234, atol=1e-15)


# ---------
# Utilities
# ---------
@pytest.fixture
def engdb_jw01029():
    """Setup the test engineering database"""
    with EngDB_Mocker(db_path=db_1029):
        engdb = engdb_tools.ENGDB_Service(base_url='http://localhost')
        yield engdb


@pytest.fixture
def mast():
    """Use the Mast database."""
    try:
        engdb = engdb_mast.EngdbMast(base_url=engdb_mast.MAST_BASE_URL)
    except RuntimeError as exception:
        pytest.skip(f'Live MAST Engineering Service not available: {exception}')
    yield engdb


@pytest.fixture
def multi_fixture(request):
    """Allow a test to use multiple fixtures"""
    return request.getfixturevalue(request.param)


def make_level1b():
    data = array([1.])
    data.shape = (1, 1, 1, 1)
    model = Level1bModel(data)
    model.update(WCS_META)
    return model
