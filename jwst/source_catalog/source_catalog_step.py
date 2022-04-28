"""
Module for the source catalog step.
"""

import os

from .detection import JWSTSourceFinder
from .reference_data import ReferenceData
from .source_catalog import JWSTSourceCatalog
from .. import datamodels
from ..stpipe import Step

__all__ = ["SourceCatalogStep"]


class SourceCatalogStep(Step):
    """
    Create a final catalog of source photometry and morphologies.

    Parameters
    -----------
    input : str or `ImageModel`
        A FITS filename or an `ImageModel` of a drizzled image.
    """

    spec = """
        bkg_boxsize = integer(default=100)    # background mesh box size in pixels
        kernel_fwhm = float(default=2.0)      # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=3.0)    # SNR threshold above the bkg
        npixels = integer(default=5)          # min number of pixels in source
        deblend = boolean(default=False)      # deblend sources?
        aperture_ee1 = integer(default=30)    # aperture encircled energy 1
        aperture_ee2 = integer(default=50)    # aperture encircled energy 2
        aperture_ee3 = integer(default=70)    # aperture encircled energy 3
        ci1_star_threshold = float(default=2.0)  # CI 1 star threshold
        ci2_star_threshold = float(default=1.8)  # CI 2 star threshold
        suffix = string(default='cat')        # Default suffix for output files
    """

    reference_file_types = ['apcorr', 'abvegaoffset']

    def _get_reffile_paths(self, model):
        filepaths = []
        for reffile_type in self.reference_file_types:
            filepath = self.get_reference_file(model, reffile_type)
            self.log.info(f'Using {reffile_type.upper()} reference file: '
                          f'{filepath}')
            filepaths.append(filepath)
        return filepaths

    def process(self, input_model):
        with datamodels.open(input_model) as model:
            reffile_paths = self._get_reffile_paths(model)
            aperture_ee = (self.aperture_ee1, self.aperture_ee2,
                           self.aperture_ee3)
            refdata = ReferenceData(input_model, reffile_paths, aperture_ee)
            aperture_params = refdata.aperture_params
            abvega_offset = refdata.abvega_offset

            finder = JWSTSourceFinder(self.bkg_boxsize, self.kernel_fwhm,
                                      self.snr_threshold, self.npixels,
                                      deblend=self.deblend)
            segment_img = finder(model)
            if segment_img is None:
                return

            ci_star_thresholds = (self.ci1_star_threshold,
                                  self.ci2_star_threshold)
            catobj = JWSTSourceCatalog(model, segment_img,
                                       kernel_fwhm=self.kernel_fwhm,
                                       aperture_params=aperture_params,
                                       abvega_offset=abvega_offset,
                                       ci_star_thresholds=ci_star_thresholds)
            catalog = catobj.catalog

            if self.save_results:
                cat_filepath = self.make_output_path(ext='.ecsv')
                catalog.write(cat_filepath, format='ascii.ecsv',
                              overwrite=True)
                model.meta.source_catalog = os.path.basename(cat_filepath)
                self.log.info(f'Wrote source catalog: {cat_filepath}')

                segm_model = datamodels.SegmentationMapModel(segment_img.data)
                segm_model.update(model, only="PRIMARY")
                segm_model.meta.wcs = model.meta.wcs
                segm_model.meta.wcsinfo = model.meta.wcsinfo
                self.save_model(segm_model, suffix='segm')
                model.meta.segmentation_map = segm_model.meta.filename
                self.log.info('Wrote segmentation map: '
                              f'{segm_model.meta.filename}')

        return catalog
