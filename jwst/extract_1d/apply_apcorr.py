from typing import Union, Tuple

import numpy as np

from astropy.io import fits
from scipy.interpolate import interp2d

from ..assign_wcs.util import compute_scale
from ..datamodels import DataModel, ReferenceFileModel


class ApCorr:
    """Perform the appropriate aperture correction on input extraction data.

    """
    def __init__(self, input_model: DataModel, appcorr_model: ReferenceFileModel, slitname: str = None,
                 location: Tuple[float, float] = None):
        # TODO: Need to keep the input_model to access wcs and whatnot for calculating pixel scale if needed
        self.reference = None
        self.correction = None

        self.model = input_model
        self._reference_table = appcorr_model.appcorr_table
        self.slitname = slitname
        self.location = location
        self.match_keys = self._get_match_keys()

        self.match_pars = {
            key: getattr(input_model.meta, key) for key in self.match_keys if hasattr(input_model.meta, key)
        }

        if self.slitname:
            self.match_pars['slit'] = self.slitname

        self._reduce_reftable()

        if self.reference['size'].unit.startswith('arcsec'):
            if self.location is not None:
                self.reference['size'] /= compute_scale(self.model.wcs, location)
            else:
                raise ValueError(
                    f'If the size column for the input APCORR reference file are in units with arcseconds, a location '
                    f'(RA, DEC) must be provided in order to convert to pixels.'
                )

        self.apcorr_func, self.apcorr_err_func = self._approx_apcorr_fn()

    def _get_match_keys(self):
        match_pars = {
            'MIRI': {
                'LRS': ['subarray'],
                'MRS': ['radius', 'axis_ratio', 'axis_pa']
            },
            'NIRSPEC': {
                'IFU/MOS': ['filter', 'grating', 'pixphase'],
                'FS': ['filter', 'grating', 'slit', 'pixphase'],
            },
            'NIRCAM': {
                'WFSS': ['filter', 'pupil']
            },
            'NIRISS': {
                'WFSS': ['filter', 'pupil']
            }
        }

        instrument = self.model.meta.instrument.name.upper()
        mode = self.model.meta.mode.name.upper()

        relevant_pars = match_pars[instrument]

        if mode in ('IFU', 'MOS'):
            return relevant_pars['IFU/MOS']

        return relevant_pars[mode]

    def _reduce_reftable(self):
        table = self._reference_table.copy()

        for key, value in self.match_pars:
            table = table[table['key'] == value]

        self.reference = table

    def _approx_apcorr_fn(self):
        return (
            interp2d(self.reference['wavelength'], self.reference['size'], self.reference['apcorr']),
            interp2d(self.reference['wavelength'], self.reference['size'], self.reference['apcorr_err'])
        )

    def apply_apcorr(self, spec_table):
        # TODO: multiply derived correction  with relevant cols such as flux, err, surf_bright, sb_error, etc
        ...

# def inperpolate_apcorr(wavelength, apcorr_wl, size, size_units, **compute_scale_kwargs):
#
#     if size_units.startswith('arcsec'):
#         # This assumes that the image scale is the same in both the dispersion
#         # and cross-dispersion directions, which is probably not valid.
#         # This is the pixel scale in arcseconds per pixel.
#         # pixel_scale = compute_scale(**compute_scale_kwargs)
#         # size /= pixel_scale                     # convert to pixels
#
#     wavelength[np.isnan(wavelength)] = -1.
#     extr_size = np.interp(wavelength, apcorr_wl, size,
#                           left=np.NaN, right=np.NaN)
#     no_cal = np.isnan(extr_size)
