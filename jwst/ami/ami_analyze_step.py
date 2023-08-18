from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import ami_analyze

__all__ = ["AmiAnalyzeStep"]


class AmiAnalyzeStep(Step):
    """Performs analysis of an AMI mode exposure by applying the LG algorithm."""

    class_alias = "ami_analyze"

    spec = """
        oversample = integer(default=3, min=1)  # Oversampling factor
        rotation = float(default=0.0)           # Rotation initial guess [deg]
        psf_offset = string(default='0.0 0.0') # PSF offset values to use to create the model array
        rotation_search = string(default='-3 3 1') # Rotation search parameters: start, stop, step
        affine2d = None # User-defined Affine2d object
        src = string(default='A0V') # Source spectral type for model
        bandpass = None # Synphot spectrum or numpy array to override filter/source
        usebp = bool(default=True) # If True, exclude pixels marked DO_NOT_USE from fringe fitting
        firstfew = None # analyze only first few integrations
        chooseholes = None # fit only certain fringes e.g. ['B4','B5','B6','C2']
    """

    #reference_file_types = ['throughput']

    def process(self, input):
        """
        Performs analysis of an AMI mode exposure by applying the LG algorithm.

        Parameters
        ----------
        input: string
            input file name

        Returns
        -------
        oifitsmodel: AmiOIModel object
            AMI tables of median observables from LG algorithm fringe fitting in OIFITS format
        oifitsmodel_multi: AmiOIModel object
            AMI tables of observables for each integration from LG algorithm fringe fitting in OIFITS format
        amilgmodel: AmiLGFitModel object
            AMI cropped data, model, and residual data from LG algorithm fringe fitting
        """
        # Retrieve the parameter values
        oversample = self.oversample
        rotate = self.rotation
        src = self.src
        bandpass = self.bandpass
        usebp = self.usebp
        firstfew = self.firstfew
        chooseholes = self.chooseholes
        affine2d = self,affine2d

        # pull out parameters that are strings and change to floats
        psf_offset = [float(a) for a in self.psf_offset.split()]
        rotsearch_parameters = [float(a) for a in self.rotation_search.split()]

        self.log.info(f"Oversampling factor = {oversample}")
        self.log.info(f"Initial rotation guess = {rotate} deg")
        self.log.info(f"Initial values to use for psf offset = {psf_offset}")

        # Open the input data model. Can be 2D or 3D image, so use general DataModel
        try:
            input_model = datamodels.DataModel(input)
        except ValueError as err:
            raise RuntimeError(f"{err}. Input unable to be read into a DataModel.")

        # Apply the LG+ methods to the data
        result = ami_analyze.apply_LG_plus(input_model,
                                           oversample, rotate,
                                           psf_offset,
                                           rotsearch_parameters,
                                           src, bandpass, usebp, firstfew, chooseholes, affine2d
                                           )

        # Close the reference file and update the step status
        result.meta.cal_step.ami_analyze = 'COMPLETE'

        return result
