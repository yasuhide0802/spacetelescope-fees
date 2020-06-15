import os

import pytest
from astropy.io.fits.diff import FITSDiff

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

"""
nrs1_group_subarray.fits                the input (uncal) file
nrs1_group_subarray_group_scale.fits    output from group_scale
nrs1_group_subarray_rate.fits           output
"""


@pytest.fixture(scope="module")
def run_pipeline(jail, rtdata_module):
    """Run calwebb_detector1 pipeline on NIRSpec subarray data."""
    rtdata = rtdata_module
    rtdata.get_data("nirspec/fs/nrs1_group_subarray.fits")

    collect_pipeline_cfgs('config')
    args = ["config/calwebb_detector1.cfg", rtdata.input,
            "--steps.group_scale.save_results=True"]
    Step.from_cmdline(args)

    return rtdata


@pytest.mark.xfail(reason='known crds bestref error due to rmap')
@pytest.mark.bigdata
@pytest.mark.parametrize("output", [
    'nrs1_group_subarray_group_scale.fits',
    'nrs1_group_subarray_rate.fits',],
    ids=['group_scale', 'rate'])
def test_nirspec_detector1(run_pipeline, fitsdiff_default_kwargs, output):
    """
    Regression test of calwebb_detector1 pipeline performed on NIRSpec data.
    """
    rtdata = run_pipeline
    rtdata.output = output
    rtdata.get_truth(os.path.join("truth/test_nirspec_subarray", output))

    diff = FITSDiff(rtdata.output, rtdata.truth, **fitsdiff_default_kwargs)
    assert diff.identical, diff.report()
