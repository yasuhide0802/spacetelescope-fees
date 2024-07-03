#! /usr/bin/env python
import gc
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import dq_initialization
from jwst.lib.basic_utils import use_datamodel, copy_datamodel


__all__ = ["DQInitStep"]


class DQInitStep(Step):
    """Initialize the Data Quality extension from the
    mask reference file.

    The dq_init step initializes the pixeldq attribute of the
    input datamodel using the MASK reference file.  For some
    FGS exp_types, initialize the dq attribute of the input model
    instead.  The dq attribute of the MASK model is bitwise OR'd
    with the pixeldq (or dq) attribute of the input model.
    """

    class_alias = "dq_init"

    spec = """
    """
    reference_file_types = ['mask']

    def process(self, input):
        """Perform the dq_init calibration step

        Parameters
        ----------
        input : JWST datamodel
            input jwst datamodel

        Returns
        -------
        output_model : JWST datamodel
            result JWST datamodel
        """

        # Try to open the input as a regular RampModel
        try:
            input_model = use_datamodel(input, model_class=datamodels.RampModel)
            # Check to see if it's Guider raw data
            if input_model.meta.exposure.type in dq_initialization.guider_list:
                # Reopen as a GuiderRawModel
                input_model.close()
                input_model = datamodels.GuiderRawModel(input)
                self.log.info("Input opened as GuiderRawModel")

        except (TypeError, ValueError):
            # If the initial open attempt fails,
            # try to open as a GuiderRawModel
            try:
                input_model = datamodels.GuiderRawModel(input)
                self.log.info("Input opened as GuiderRawModel")
            except (TypeError, ValueError):
                self.log.error("Unexpected or unknown input model type")
        except Exception:
            self.log.error("Can't open input")
            raise

        result, input_model = copy_datamodel(input_model, self.modify_input)
        del input

        # Retrieve the mask reference file name
        self.mask_filename = self.get_reference_file(result, 'mask')
        self.log.info('Using MASK reference file %s', self.mask_filename)

        # Check for a valid reference file
        if self.mask_filename == 'N/A':
            self.log.warning('No MASK reference file found')
            self.log.warning('DQ initialization step will be skipped')
            result.meta.cal_step.dq_init = 'SKIPPED'
            return result

        # Load the reference file
        mask_model = datamodels.MaskModel(self.mask_filename)

        # Apply the step
        result = dq_initialization.correct_model(result, mask_model)

        # Close the data models for the input and ref file
        del mask_model
        gc.collect()

        return result
