"""Master Background Pipeline for applying Master Background to NIRSpec Slit-like data"""

from ..barshadow import barshadow_step
from .. import datamodels
from ..flatfield import flat_field_step
from ..master_background import nirspec_utils
from ..pathloss import pathloss_step
from ..photom import photom_step
from ..stpipe import Pipeline

# Step parameters to generally ignore when copying from the parent steps.
GLOBAL_PARS_TO_IGNORE = ['output_ext', 'output_file', 'output_use_model', 'output_use_index',
                         'inverse', 'pre_hooks', 'post_hooks', 'save_results', 'suffix']


class MasterBackgroundNRSSlitsPipe(Pipeline):
    """Apply master background processing to NIRSpec Slit-like data

    For MOS, and ignoring FS, the calibration process needs to occur
    twice: Once to calibrate background slits and create a master background.
    Then a second time to calibrate science using the master background.

    Notes
    -----
    The algorithm is as follows:

    - Calibrate all slits
      - For each step:
        - Force the source type to be extended source for all slits.
        - Return the correction array used.
    - Create the 1D master background
    - For each slit
      - Expand out the 1D master background to match the 2D wavelength grid of the slit
      - Reverse-calibrate the 2D background, using the correction arrays calculated above.
      - Subtract the background from the input slit data
    """

    spec = """
        force_subtract = boolean(default=False)  # Force subtracting master background
        save_background = boolean(default=False) # Save computed master background
        user_background = string(default=None)   # Path to user-supplied master background
        inverse = boolean(default=False)    # Invert the operation
        output_use_model = boolean(default=True)
    """

    # Define aliases to steps
    step_defs = {
        'flat_field': flat_field_step.FlatFieldStep,
        'pathloss': pathloss_step.PathLossStep,
        'barshadow': barshadow_step.BarShadowStep,
        'photom': photom_step.PhotomStep,
    }

    def process(self, data):
        """Compute and subtract a master background spectrum

        Parameters
        ----------
        data : `~jwst.datamodels.MultiSlitModel`
            The data to operate on.

        user_background : None, string, or `~jwst.datamodels.MultiSpecModel`
            Optional user-supplied master background 1D spectrum, path to file
            or opened datamodel

        save_background : bool, optional
            Save computed master background.

        force_subtract : bool, optional
            Optional user-supplied flag that overrides step logic to force subtraction of the
            master background.
            Default is False, in which case the step logic determines if the calspec2 background step
            has already been applied and, if so, the master background step is skipped.
            If set to True, the step logic is bypassed and the master background is subtracted.

        Attributes
        ----------
        correction_pars : dict
            The master background information from a previous invocation of the step.
            Keys are:

            - "masterbkg_1d": `~jwst.datamodels.CombinedSpecModel`
                The 1D version of the master background.

            - "masterbkg_2d": `~jwst.datamodels.MultiSlitModel`
                The 2D slit-based version of the master background.

        Returns
        -------
        result : `~jwst.datamodels.MultiSlitModel`
        """

        with datamodels.open(data) as data_model:
            # If some type of background processing had already been done. Abort.
            # UNLESS forcing is enacted.
            if not self.force_subtract and \
               'COMPLETE' in [data_model.meta.cal_step.back_sub, data_model.meta.cal_step.master_background]:
                self.log.info('Background subtraction has already occurred. Skipping.')
                data.meta.cal_step.master_background = 'SKIP'
                return data

            if self.use_correction_pars:
                self.log.info('Using pre-calculated correction parameters.')
                master_background = self.correction_pars['masterbkg_1d']
                mb_multislit = self.correction_pars['masterbkg_2d']
            else:
                self.log.info('Calculating master background')
                master_background, mb_multislit = self._calc_master_background(data_model)

            # Now apply the de-calibrated background to the original science
            result = nirspec_utils.apply_master_background(data_model, mb_multislit, inverse=self.inverse)

            # Mark as completed and setup return data
            result.meta.cal_step.master_background = 'COMPLETE'
            self.correction_pars = {
                'masterbkg_1d': master_background,
                'masterbkg_2d': mb_multislit
            }

        return result

    def set_step_pars(self):
        """Set substep parameters from the parents substeps"""
        if not self.parent:
            return

        steps = ['barshadow', 'flat_field', 'pathloss', 'photom']
        pars_to_ignore = {
            'barshadow': ['source_type'],
            'flat_field': ['save_interpolated_flat'],
            'pathloss': ['source_type'],
            'photom': ['source_type']
        }

        for step in steps:
            pars = getattr(self.parent, step).get_pars()
            for par in pars_to_ignore[step] + GLOBAL_PARS_TO_IGNORE:
                del pars[par]
            getattr(self, step).update(pars)

    def _calc_master_background(self, data):
        """Calculate master background from background slits

        Parameters
        ----------
        data : `~jwst.datamodels.MultiSlitModel`
            The data to operate on.

        Returns
        -------
        masterbkg_1d, masterbkg_2d : `~jwst.datamodels.CombinedSpecModel`, `~jwst.datamodels.MultiSlitModel`
            The master background in 1d and 2d, multislit formats.
        """

        # Set relevant parameters from the parent version of the steps
        self.set_step_pars()
        saved_pars = self.get_pars()

        # First pass: just do the calibration to determine the correction
        # arrays. However, force all slits to be processed as extended sources.
        self.pathloss.source_type = 'EXTENDED'
        self.barshadow.source_type = 'EXTENDED'
        self.photom.source_type = 'EXTENDED'

        pre_calibrated = self.flat_field(data)
        pre_calibrated = self.pathloss(pre_calibrated)
        pre_calibrated = self.barshadow(pre_calibrated)
        pre_calibrated = self.photom(pre_calibrated)

        # Create the 1D, fully calibrated master background.
        master_background = nirspec_utils.create_background_from_multislit(pre_calibrated)
        if master_background is None:
            self.log.info('No master background could be calculated. Skipping.')
            data.meta.cal_step.master_background = 'SKIP'
            return data

        # Now decalibrate the master background for each individual science slit.
        # First step is to map the master background into a MultiSlitModel
        # where the science slits are replaced by the master background.
        # Here the broadcasting from 1D to 2D need also occur.
        mb_multislit = nirspec_utils.map_to_science_slits(pre_calibrated, master_background)

        # Now that the master background is pretending to be science,
        # walk backwards through the steps to uncalibrate, using the
        # calibration factors carried from `pre_calibrated`.
        self.photom.use_correction_pars = True
        self.photom.inverse = True
        self.barshadow.use_correction_pars = True
        self.barshadow.inverse = True
        self.pathloss.use_correction_pars = True
        self.pathloss.inverse = True
        self.flat_field.use_correction_pars = True
        self.flat_field.inverse = True

        mb_multislit = self.photom(mb_multislit)
        mb_multislit = self.barshadow(mb_multislit)
        mb_multislit = self.pathloss(mb_multislit)
        mb_multislit = self.flat_field(mb_multislit)

        self.update(saved_pars)
        return master_background, mb_multislit
