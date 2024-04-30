#! /usr/bin/env python
from functools import partial

from ..stpipe import Step
from .. import datamodels
from .pixel_replace import PixelReplacement

__all__ = ["PixelReplaceStep"]


class PixelReplaceStep(Step):
    """
    PixelReplaceStep: Module for replacing flagged bad pixels with an estimate
    of their flux, prior to spectral extraction.

    Attributes
    ----------
    algorithm : str
        Method used to estimate flux values for bad pixels. Currently only one option is
        implemented, using a profile fit to adjacent column values.

    n_adjacent_cols : int
        Number of adjacent columns (on either side of column containing a bad pixel) to use in
        creation of source profile, in cross-dispersion direction. The total number of columns
        used in the profile will be twice this number; on array edges, take adjacent columns until
        this number is reached.
    """

    class_alias = "pixel_replace"

    spec = """
        algorithm = option("fit_profile", "mingrad", "N/A", default="fit_profile")
        n_adjacent_cols = integer(default=3)    # Number of adjacent columns to use in creation of profile
        skip = boolean(default=True) # Step must be turned on by parameter reference or user
    """

    def process(self, input):
        """Execute the step.

        Parameters
        ----------
        input : JWST DataModel

        Returns
        -------
        JWST DataModel
            This will be `input` if the step was skipped; otherwise,
            it will be a model containing data arrays with estimated fluxes
            for any bad pixels, now flagged as TO-BE-DETERMINED (DQ bit 7?).
        """

        with datamodels.open(input) as input_model:
            # Make copy of input to prevent overwriting
            result = input_model.copy()
            # If more than one 2d spectrum exists in input, call replacement
            # for each spectrum - NRS_FIXEDSLIT, WFSS
            if isinstance(input_model, datamodels.MultiSlitModel):
                self.log.debug('Input is a MultiSlitModel.')
            # NRS_BRIGHTOBJ provides a SlitModel
            elif isinstance(input_model, datamodels.SlitModel):
                self.log.debug('Input is a SlitModel.')
            # MIRI_LRS-FIXEDSLIT comes in ImageModel - any others?
            elif isinstance(input_model, datamodels.ImageModel):
                self.log.debug('Input is an ImageModel.')
            elif isinstance(input_model, datamodels.IFUImageModel):
                self.log.debug('Input is an IFUImageModel.')
            # SOSS & LRS-SLITLESS have CubeModel, along with IFU modes WIP
            elif isinstance(input_model, datamodels.CubeModel):
                # It's a 3-D multi-integration model
                self.log.debug('Input is a CubeModel for a multiple integration file.')
            elif isinstance(input_model, datamodels.ModelContainer):
                self.log.debug('Input is an ModelContainer.')
            else:
                self.log.error(f'Input is of type {str(type(input_model))} for which')
                self.log.error('pixel_replace does not have an algorithm.\n')
                self.log.error('Pixel replacement will be skipped.')
                input_model.meta.cal_step.pixel_replace = 'SKIPPED'
                return input_model

            pars = {
                'algorithm': self.algorithm,
                'n_adjacent_cols': self.n_adjacent_cols,
            }

            #___________________
            # calewbb_spec3 case
            #___________________
            if isinstance(input_model, datamodels.ModelContainer):
                output_model = input_model.copy()                
                # Setup output path naming if associations are involved.
                asn_id = None
                try:
                    asn_id = input_model.meta.asn_table.asn_id
                except (AttributeError, KeyError):
                    pass
                if asn_id is None:
                    asn_id = self.search_attr('asn_id')
                if asn_id is not None:
                    _make_output_path = self.search_attr(
                        '_make_output_path', parent_first=True
                    )
                    self._make_output_path = partial(
                        _make_output_path,
                        asn_id=asn_id
                    )
                # Check models to confirm they are the correct type
                for i, model in enumerate(output_model):
                    run_pixel_replace = True
                    # for each spectrum - NRS_FIXEDSLIT, WFSS
                    if isinstance(model, datamodels.MultiSlitModel):
                        self.log.debug('Input is a MultiSlitModel.')
                        # NRS_BRIGHTOBJ provides a SlitModel
                    elif isinstance(model, datamodels.SlitModel):
                        self.log.debug('Input is a SlitModel.')
                        # MIRI_LRS-FIXEDSLIT comes in ImageModel - any others?
                    elif isinstance(model, datamodels.ImageModel):
                        self.log.debug('Input is an ImageModel.')
                    elif isinstance(model, datamodels.IFUImageModel):
                        self.log.debug('Input is an IFUImageModel.')
                        # SOSS & LRS-SLITLESS have CubeModel, along with IFU modes WIP
                    elif isinstance(model, datamodels.CubeModel):
                        # It's a 3-D multi-integration model
                        self.log.debug('Input is a CubeModel for a multiple integration file.')
                    else:
                        self.log.error(f'Input is of type {str(type(input_model))} for which')
                        self.log.error('pixel_replace does not have an algorithm.\n')
                        self.log.error('Pixel replacement will be skipped.')
                        model.meta.cal_step.pixel_replace = 'SKIPPED'
                        run_pixel_replace = False

                    # all checks on model have passed. Now run pixel replacement
                    if run_pixel_replace:
                        replacement = PixelReplacement(model, **pars)
                        replacement.replace()
                        model = replacement.output
                        self.record_step_status(replacement.output, 'pixel_replace', success=True)

                return output_model
            #________________________________________
            # calewbb_spec2 case - single input model
            #________________________________________
            else:
                replacement = PixelReplacement(result, **pars)
                replacement.replace()
                self.record_step_status(replacement.output, 'pixel_replace', success=True)
                result = replacement.output
                return result
