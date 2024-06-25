from stdatamodels.jwst import datamodels

from ..stpipe import Step
from stcal.dark_current import dark_sub
from jwst.lib.basic_utils import use_datamodel
import numpy as np


__all__ = ["DarkCurrentStep"]


class DarkCurrentStep(Step):
    """
    DarkCurrentStep: Performs dark current correction by subtracting
    dark current reference data from the input science data model.
    """

    class_alias = "dark_current"

    spec = """
        dark_output = output_file(default = None) # Dark model or averaged dark subtracted
        average_dark_current = float(default=None) # The average dark current for this detector in units of e-/sec.
    """

    reference_file_types = ['dark']

    def process(self, input):

        # Open the input data model
        with use_datamodel(input, model_class=datamodels.RampModel) as input_model:

            # Get the name of the dark reference file to use
            self.dark_name = self.get_reference_file(input_model, 'dark')
            self.log.info('Using DARK reference file %s', self.dark_name)

            # Check for a valid reference file
            if self.dark_name == 'N/A':
                self.log.warning('No DARK reference file found')
                self.log.warning('Dark current step will be skipped')
                result = input_model
                result.meta.cal_step.dark = 'SKIPPED'
                return result

            # Create name for the intermediate dark, if desired.
            dark_output = self.dark_output
            if dark_output is not None:
                dark_output = self.make_output_path(
                    basepath=dark_output,
                    suffix=False
                )

            # Open the dark ref file data model - based on Instrument
            instrument = input_model.meta.instrument.name
            if instrument == 'MIRI':
                dark_model = datamodels.DarkMIRIModel(self.dark_name)
            else:
                dark_model = datamodels.DarkModel(self.dark_name)

            # Store user-defined average_dark_current in model, if provided
            # A user-defined value will take precedence over any value present
            # in dark reference file
            self.set_average_dark_current(input_model, dark_model)

            # Do the dark correction
            result = dark_sub.do_correction(
                input_model, dark_model, dark_output
            )

            out_data, dark_data = result

            if dark_data is not None and dark_data.save:
                save_dark_data_as_dark_model(dark_data, dark_model, instrument)
            dark_model.close()

            out_ramp = dark_output_data_2_ramp_model(out_data, input_model)

        return out_ramp

    def set_average_dark_current(self, input_model, dark_model):
        """Take the three possible locations specifying
        the average dark current and assign them to the
        input model, in priority order:
        1) Any value provided to the step parameter, either from
        the user or a parameter reference file
        2) The 2-D array stored in dark_model.average_dark_current
        3) The scalar value stored in dark_model.meta.exposure.average_dark_current

        Inputs
        ------
        input_model : stdatamodels.jwst.datamodels.RampModel
            The input datamodel containing the 4-D ramp array
        dark_model : Union[stdatamodels.jwst.datamodels.DarkModel, stdatamodels.jwst.datamodels.DarkMIRIModel]
            The dark reference file datamodel
        """
        if self.average_dark_current is not None:
            input_model.average_dark_current[:, :] = self.average_dark_current
            self.log.info('Using Poisson noise from average dark current %s e-/sec', self.average_dark_current)
        else:
            # First prioritize a 2D average_dark_current, if defined in dark.
            # If not present, apply scalar to datamodel array, if scalar is present.
            if np.sum(dark_model.average_dark_current) == 0.0:
                input_model.average_dark_current[:, :] = dark_model.meta.exposure.average_dark_current
            elif np.shape(input_model.average_dark_current) != np.shape(dark_model.average_dark_current):
                self.log.warning("DarkModel average_dark_current does not match shape of data.\n"
                                 "Dark current from reference file cannot be applied.")
            else:
                input_model.average_dark_current = dark_model.average_dark_current


def save_dark_data_as_dark_model(dark_data, dark_model, instrument):
    """
    Save dark data from the dark current step as the appropriate dark model.

    Parameters
    ----------
    dark_data: DarkData
        Dark data used in the dark current step.

    dark_model: DarkMIRIModel or DarkModel
        The input dark model from reference.

    instrument: str
        The instrument name.
    """
    if instrument == "MIRI":
        out_dark_model = datamodels.DarkMIRIModel(
            data=dark_data.data,
            dq=dark_data.groupdq,
            err=dark_data.err)
    else:
        out_dark_model = datamodels.DarkModel(
            data=dark_data.data,
            dq=dark_data.groupdq,
            err=dark_data.err)
    out_dark_model.update(dark_model)

    out_dark_model.meta.exposure.nframes = dark_data.exp_nframes
    out_dark_model.meta.exposure.ngroups = dark_data.exp_ngroups
    out_dark_model.meta.exposure.groupgap = dark_data.exp_groupgap
    out_dark_model.save(dark_data.output_name)
    out_dark_model.close()


def dark_output_data_2_ramp_model(out_data, input_model):
    """
    Convert computed output data from the dark step to a RampModel.

    Parameters
    ----------
    out_data: ScienceData
        Computed science data from the dark current step.

    input_model: RampModel
        The input ramp model from which to subtract the dark current.

    Return
    ------
    out_model: RampModel
        The output ramp model from the dark current step.
    """

    if out_data.cal_step == "SKIPPED":
        # If processing was skipped in the lower-level routines,
        # just return the unmodified input model
        input_model.meta.cal_step.dark_sub = "SKIPPED"
        return input_model
    else:
        out_model = input_model.copy()
        out_model.meta.cal_step.dark_sub = out_data.cal_step
        out_model.data = out_data.data
        out_model.groupdq = out_data.groupdq
        out_model.pixeldq = out_data.pixeldq
        out_model.err = out_data.err
        return out_model
