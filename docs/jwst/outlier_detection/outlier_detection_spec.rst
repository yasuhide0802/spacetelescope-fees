.. _outlier-detection-spec:

Outlier Detection for Slit-like Spectroscopic Data
==================================================

:Classes: `jwst.outlier_detection.OutlierDetectionSpecStep`
:Aliases: outlier_detection_spec

This module serves as the interface for applying ``outlier_detection`` to slit-like
spectroscopic observations.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST
spectroscopic observations.

Specifically, this routine performs the following operations (modified from the
:ref:`Default Outlier Detection Algorithm <outlier-detection-imaging>` ):

#. Extract parameter settings from input model and merge them with any user-provided values

   - the same set of parameters available to:
     ref:`Default Outlier Detection Algorithm <outlier-detection-imaging>`
     also applies to this code
#. Convert input data, as needed, to make sure it is in a format that can be processed

   - A :py:class:`~jwst.datamodels.ModelContainer` serves as the basic format 
     for all processing performed by
     this step, as each entry will be treated as an element of a stack of images
     to be processed to identify bad pixels, cosmic-rays and other artifacts
   - If the input data is a :py:class:`~jwst.datamodels.CubeModel`, convert it into a 
     :py:class:`~jwst.datamodels.ModelContainer`.
     This allows each plane of the cube to be treated as a separate 2D image
     for resampling (if done) and for combining into a median image.
#. Resample all input images into a :py:class:`~jwst.datamodels.ModelContainer` using
   :py:class:`~jwst.resample.resample_spec.ResampleSpecData`

   - Resampled images are written out to disk if the ``save_intermediate_results``
     parameter is set to `True`
   - **If resampling is turned off**, the original unrectified inputs are used to create
     the median image for cosmic-ray detection
#. Create a median image from (possibly) resampled :py:class:`~jwst.datamodels.ModelContainer`

   - The median image is written out to disk if the ``save_intermediate_results``
     parameter is set to `True`
#. Blot median image to match each original input image

   - **If resampling is turned off**, the median image is used for comparison
     with the original input models for detecting outliers
#. Perform statistical comparison between blotted image and original image to identify outliers
#. Update input data model DQ arrays with mask of detected outliers


Step Arguments
==============

``--weight_type``
  The type of data weighting to use during resampling.

``--pixfrac``
  The pixel fraction used during resampling;
  valid values go from 0.0 to 1.0.

``--kernel``
  The form of the kernel function used to distribute flux onto a
  resampled image.

``--fillval``
  The value to assign to resampled image pixels that have zero weight or
  do not receive any flux from any input pixels during drizzling.
  Any floating-point value, given as a string, is valid.
  A value of 'INDEF' will use the last zero weight flux.

``--maskpt``
  The percent of maximum weight to use as lower-limit for valid data;
  valid values go from 0.0 to 1.0.

``--snr``
  The signal-to-noise values to use for bad pixel identification.
  Since cosmic rays often extend across several pixels the user
  must specify two cut-off values for determining whether a pixel should
  be masked: the first for detecting the primary cosmic ray, and the
  second (typically lower threshold) for masking lower-level bad pixels
  adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example "5.0 4.0").

``--scale``
  The scaling factor applied to derivative used to identify bad pixels.
  Since cosmic rays often extend across several pixels the user
  must specify two cut-off values for determining whether a pixel should
  be masked: the first for detecting the primary cosmic ray, and the
  second (typically lower threshold) for masking lower-level bad pixels
  adjacent to those found in the first pass.  Valid values are a pair of
  floating-point values in a single string (for example "1.2 0.7").

``--backg``
  User-specified background value to apply to the median image.

``--save_intermediate_results``
  Specifies whether or not to save any intermediate products created
  during step processing.

``--resample_data``
  Specifies whether or not to resample the input images when
  performing outlier detection.

``--good_bits``
  The DQ bit values from the input image DQ arrays
  that should be considered 'good' when building the weight mask. See
  DQ flag :ref:`dq_parameter_specification` for details.

``--in_memory``
  Specifies whether or not to load and create all images that are used during
  processing into memory. If ``False``, input files are loaded from disk when
  needed and all intermediate files are stored on disk, rather than in memory.


Reference Files
===============

The ``outlier_detection_spec`` step uses the PARS-OUTLIERDETECTIONSPECSTEP parameter reference file.

.. include:: ../references_general/pars-outlierdetectionspecstep_reffile.inc

.. automodapi:: jwst.outlier_detection.outlier_detection_spec_step
