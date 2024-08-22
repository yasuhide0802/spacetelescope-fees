"""
Submodule for performing outlier detection on spectra.
"""
import copy

from stdatamodels.jwst import datamodels

from jwst.datamodels import ModelContainer
from jwst.stpipe.utilities import record_step_status
from jwst.stpipe import Step

from ..resample import resample_spec, resample_utils
from .utils import create_median, flag_crs_in_models, flag_crs_in_models_with_resampling, OutlierDetectionStepBase
from ._fileio import remove_file

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["OutlierDetectionSpecStep"]


class OutlierDetectionSpecStep(Step, OutlierDetectionStepBase):
    """Flag outlier bad pixels and cosmic rays in DQ array of each input image.

    Input images can be listed in an input association file or already opened
    with a ModelContainer.  DQ arrays are modified in place.

    Parameters
    -----------
    input_data : ~jwst.datamodels.ModelContainer or str
        ModelContainer or association file that can initialize a ModelContainer.
    """

    class_alias = "outlier_detection_spec"

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

        if self.resample_data is True:
            # Start by creating resampled/mosaic images for
            #  each group of exposures
            resamp = resample_spec.ResampleSpecData(
                input_models,
                single=True,
                blendheaders=False,
                wht_type=self.weight_type,
                pixfrac=self.pixfrac,
                kernel=self.kernel,
                fillval=self.fillval,
                good_bits=self.good_bits,
                in_memory=self.in_memory,
                asn_id=asn_id,
            )
            median_wcs = resamp.output_wcs
            drizzled_models = resamp.do_drizzle(input_models)
            if self.save_intermediate_results:
                for model in drizzled_models:
                    model.meta.filename = self.make_output_path(
                        basepath=model.meta.filename,
                        suffix="_outlier_s2d.fits",
                    )
                    log.info("Writing out resampled spectra...")
                    model.save(model.meta.filename)
        else:
            drizzled_models = input_models
            for i in range(len(input_models)):
                drizzled_models[i].wht = resample_utils.build_driz_weight(
                    input_models[i],
                    weight_type=self.weight_type,
                    good_bits=self.good_bits)
            # copy for when saving median and input is a filename?
            median_wcs = copy.deepcopy(input_models[0].meta.wcs)

        # Perform median combination on set of drizzled mosaics
        # create_median should be called as a method from parent class
        median_data = create_median(drizzled_models, self.maskpt)

        if self.save_intermediate_results:
            # Initialize intermediate products used in the outlier detection
            median_model = datamodels.ImageModel(median_data)
            median_model.meta = drizzled_models[0].meta
            median_model.meta.filename = self.make_output_path(
                basepath=input_models[0].meta.filename,
                suffix='median'
            )

            log.info("Writing out MEDIAN image to: {}".format(
                    median_model.meta.filename))
            median_model.save(median_model.meta.filename)

            del median_model
        else:
            # since we're not saving intermediate results if the drizzled models
            # were written to disk, remove them
            if not self.in_memory:
                for fn in drizzled_models._models:
                    remove_file(fn)
                    log.info(f"Removing file {fn}")

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
