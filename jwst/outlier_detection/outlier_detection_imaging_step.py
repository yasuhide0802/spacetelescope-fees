"""
Submodule for performing outlier detection on imaging data.
"""

import copy
import logging
import os

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.resample import resample
from jwst.resample.resample_utils import build_driz_weight
from jwst.stpipe.utilities import record_step_status
from jwst.stpipe import Step

from .utils import create_median, flag_crs_in_models, flag_crs_in_models_with_resampling, OutlierDetectionStepBase
from ._fileio import remove_file, save_median

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["OutlierDetectionImagingStep"]


class OutlierDetectionImagingStep(Step, OutlierDetectionStepBase):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.

    Input images can be listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.

    Parameters
    -----------
    input_data : asn file or ~jwst.datamodels.ModelContainer
        Single filename association table, or a datamodels.ModelContainer.

    """

    class_alias = "outlier_detection_imaging"

    spec = """
        weight_type = option('ivm','exptime',default='ivm')
        pixfrac = float(default=1.0)
        kernel = string(default='square') # drizzle kernel
        fillval = string(default='INDEF')
        maskpt = float(default=0.7)
        snr = string(default='5.0 4.0')
        scale = string(default='1.2 0.7')
        backg = float(default=0.0)
        save_intermediate_results = boolean(default=False)
        resample_data = boolean(default=True)
        good_bits = string(default="~DO_NOT_USE")  # DQ flags to allow
        allowed_memory = float(default=None)  # Fraction of memory to use for the combined image
        in_memory = boolean(default=False)
    """

    def process(self, input_models):
        """Perform outlier detection processing on input data."""

        # determine the asn_id (if not set by the pipeline)
        asn_id = self._get_asn_id(input_models)
        self.log.info(f"Outlier Detection asn_id: {asn_id}")

        snr1, snr2 = [float(v) for v in self.snr.split()]
        scale1, scale2 = [float(v) for v in self.scale.split()]

        if not isinstance(input_models, ModelContainer):
            input_models = ModelContainer(input_models, save_open=self.in_memory)

        if len(input_models) < 2:
            log.warning(f"Input only contains {len(input_models)} exposures")
            log.warning("Outlier detection will be skipped")
            record_step_status(input_models, "outlier_detection", False)
            return input_models

        if self.resample_data:
            # Start by creating resampled/mosaic images for
            # each group of exposures
            output_path = self.make_output_path(basepath=input_models[0].meta.filename,
                            suffix='')
            output_path = os.path.dirname(output_path)
            resamp = resample.ResampleData(
                input_models,
                output=output_path,
                single=True,
                blendheaders=False,
                wht_type=self.weight_type,
                pixfrac=self.pixfrac,
                kernel=self.kernel,
                fillval=self.fillval,
                good_bits=self.good_bits,
                in_memory=self.in_memory,
                asn_id=asn_id,
                allowed_memory=self.allowed_memory,
            )
            median_wcs = resamp.output_wcs
            drizzled_models = resamp.do_drizzle(input_models)
        else:
            # for non-dithered data, the resampled image is just the original image
            drizzled_models = input_models
            for i in range(len(input_models)):
                drizzled_models[i].wht = build_driz_weight(
                    input_models[i],
                    weight_type=self.weight_type,
                    good_bits=self.good_bits)
            # copy for when saving median and input is a filename?
            median_wcs = copy.deepcopy(input_models[0].meta.wcs)

        # Perform median combination on set of drizzled mosaics
        median_data = create_median(drizzled_models, self.maskpt)

        if self.save_intermediate_results:
            # make a median model
            with datamodels.open(drizzled_models[0]) as dm0:
                median_model = datamodels.ImageModel(median_data)
                median_model.update(dm0)
                median_model.meta.wcs = median_wcs

            save_median(median_model, self.make_output_path, asn_id)
            del median_model
        else:
            # since we're not saving intermediate results if the drizzled models
            # were written to disk, remove them
            if not self.in_memory:
                for fn in drizzled_models._models:
                    remove_file(fn)

        # Perform outlier detection using statistical comparisons between
        # each original input image and its blotted version of the median image
        if self.resample_data:
            flag_crs_in_models_with_resampling(
                input_models,
                median_data,
                median_wcs,
                snr1,
                snr2,
                scale1,
                scale2,
                self.backg,
            )
        else:
            flag_crs_in_models(input_models, median_data, snr1)
        return input_models
