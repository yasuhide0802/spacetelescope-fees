#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from . import flat_field


class FlatFieldStep(Step):
    """
    FlatFieldStep: Flat-field a science image using a flatfield reference image.
    """

    spec = """
        # Suffix for optional output file for interpolated flat fields.
        flat_suffix = string(default=None)
    """

    # NIRSpec image data is currently not supported.  There's a check for
    # this case in the `process` method, which will then log a warning and
    # set the step switch to SKIPPED.  The line for reference_file_types
    # is commented out to prevent CRDS from raising an exception before
    # we can handle this case ourselves.  This is a temporary workaround
    # until we get actual reference files for NIRSpec image data.
    # xxx reference_file_types = ["flat", "fflat", "sflat", "dflat"]

    def process(self, input):

        input_model = datamodels.open(input)

        # Figure out what kind of input data model is in use.
        if isinstance(input_model, datamodels.CubeModel):
            # multi-integration dataset
            self.log.debug('Input is a CubeModel')
        elif isinstance(input_model, datamodels.ImageModel):
            self.log.debug('Input is an ImageModel')
        elif isinstance(input_model, datamodels.MultiSlitModel):
            self.log.debug('Input is a MultiSlitModel')

        is_NIRSpec = (input_model.meta.instrument.name == "NIRSPEC")

        # Retrieve the reference file name or names
        if is_NIRSpec:
            exp_type = input_model.meta.exposure.type
            if exp_type not in ["NRS_FIXEDSLIT", "NRS_IFU", "NRS_MSASPEC"]:
                self.log.warning("Exposure type is %s; flat-fielding will be "
                                 "skipped because it is not currently "
                                 "supported for this mode.", exp_type)
                result = input_model.copy()
                result.meta.cal_step.flat_field = "SKIPPED"
                input_model.close()
                return result
            self.f_flat_filename = self.get_reference_file(input_model,
                                        'fflat')
            self.s_flat_filename = self.get_reference_file(input_model,
                                        'sflat')
            self.d_flat_filename = self.get_reference_file(input_model,
                                        'dflat')
            self.log.debug('Using FFLAT reference file: %s',
                           self.f_flat_filename)
            self.log.debug('Using SFLAT reference file: %s',
                           self.s_flat_filename)
            self.log.debug('Using DFLAT reference file: %s',
                           self.d_flat_filename)
        else:
            self.flat_filename = self.get_reference_file(input_model, 'flat')
            self.log.debug('Using FLAT reference file: %s', self.flat_filename)

        # Check for a valid reference file
        missing = False
        if is_NIRSpec:
            if self.f_flat_filename == 'N/A' or \
               self.s_flat_filename == 'N/A' or \
               self.d_flat_filename == 'N/A':
                self.log.warning('One or more flat-field reference files'
                                 ' were missing')
                missing = True
        else:
            if self.flat_filename == 'N/A':
                self.log.warning('No FLAT reference file found')
                missing = True
        if missing:
            self.log.warning('Flat-field step will be skipped')
            result = input_model.copy()
            result.meta.cal_step.flat_field = 'SKIPPED'
            input_model.close()
            return result

        # Find out what model to use for the flat field reference file(s).
        if is_NIRSpec:
            flat_model = None
            if input_model.meta.exposure.type == "NRS_MSASPEC":
                f_flat_model = \
                    datamodels.NirspecQuadFlatModel(self.f_flat_filename)
            else:
                f_flat_model = \
                    datamodels.NirspecFlatModel(self.f_flat_filename)
            s_flat_model = datamodels.NirspecFlatModel(self.s_flat_filename)
            d_flat_model = datamodels.NirspecFlatModel(self.d_flat_filename)
        else:
            self.log.debug('Opening flat as FlatModel')
            flat_model = datamodels.FlatModel(self.flat_filename)
            f_flat_model = None
            s_flat_model = None
            d_flat_model = None

        # Do the flat-field correction
        (output_model, interpolated_flats)  = \
                flat_field.do_correction(input_model, flat_model,
                                         f_flat_model, s_flat_model,
                                         d_flat_model, self.flat_suffix)

        # Close the inputs
        input_model.close()
        if is_NIRSpec:
            f_flat_model.close()
            s_flat_model.close()
            d_flat_model.close()
        else:
            flat_model.close()

        if interpolated_flats is not None:
            self.log.info("Writing interpolated flat fields.")
            self.save_model(interpolated_flats, self.flat_suffix)
            interpolated_flats.close()

        return output_model

if __name__ == '__main__':
    cmdline.step_script(flat_field_step)
