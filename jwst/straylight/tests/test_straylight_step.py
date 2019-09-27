"""
Unit tests for straylight step configuration
"""

from jwst.datamodels import IFUImageModel
from jwst.straylight import StraylightStep
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
import numpy as np
import pytest


@pytest.fixture(scope='module')
def miri_mrs():
    """ Set up MIRI MRS Short data """

    image = IFUImageModel((20, 20))
    image.data = np.random.random((20, 20))
    image.meta.instrument.name = 'MIRI'
    image.meta.instrument.detector = 'MIRIFUSHORT'
    image.meta.exposure.type = 'MIR_MRS'
    image.meta.instrument.channel = '12'
    image.meta.instrument.band = 'SHORT'
    image.meta.filename = 'test_miri.fits'
    return image


@pytest.fixture(scope='module')
def miri_mrs_long():
    """ Set up MIRI MRS Long data """

    image = IFUImageModel((30, 30))
    image.data = np.random.random((30, 30))
    image.meta.instrument.name = 'MIRI'
    image.meta.instrument.detector = 'MIRIFULONG'
    image.meta.exposure.type = 'MIR_MRS'
    image.meta.instrument.channel = '34'
    image.meta.instrument.band = 'LONG'
    image.meta.filename = 'test_miri_long.fits'
    return image


def test_call_straylight2(_jail, miri_mrs_long):
    """ test step is skipped for MRS IFULONG data """

    collect_pipeline_cfgs('./config')
    # Test the ModShepard power is in the correct range
    step = StraylightStep.from_config_file('config/straylight.cfg')
    step.override_straylight = 'dummy.asdf'
    step.method = 'ModShepard'
    result = step.run(miri_mrs_long)
    assert result.meta.cal_step.straylight == 'SKIPPED'
