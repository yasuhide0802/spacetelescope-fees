.. _outlier-detection-imaging:

Outlier Detection Algorithm for Imaging Data
============================================

:Classes: `jwst.outlier_detection.OutlierDetectionImagingStep`
:Aliases: outlier_detection_imaging

This module serves as the interface for applying outlier detection to direct
image observations, like those taken with MIRI, NIRCam and NIRISS.  The code implements the
basic outlier detection algorithm used with HST data, as adapted to JWST.

Specifically, this routine performs the following operations:

#. Extract parameter settings from input model and merge them with any user-provided values.

#. Convert input data, as needed, to make sure it is in a format that can be processed.

   * A :py:class:`~jwst.datamodels.ModelContainer` serves as the basic format for
     all processing performed by
     this step, as each entry will be treated as an element of a stack of images
     to be processed to identify bad-pixels/cosmic-rays and other artifacts.
   * If the input data is a :py:class:`~jwst.datamodels.CubeModel`, convert it into a ModelContainer.
     This allows each plane of the cube to be treated as a separate 2D image
     for resampling (if done) and for combining into a median image.

#. By default, resample all input images.

   * The resampling step starts by computing an output WCS that is large enoug
     to encompass all the input images.
   * All images from the *same exposure* will get resampled onto this output
     WCS to create a mosaic of all the chips for that exposure.  This product
     is referred to as a "grouped mosaic" since it groups all the chips from
     the same exposure into a single image.
   * Each dither position will result in a separate grouped mosaic, so only
     a single exposure ever contributes to each pixel in these mosaics.
   * An explanation of how all NIRCam multiple detector group mosaics are
     defined from `a single exposure or from a dithered set of exposures
     <https://jwst-docs.stsci.edu/near-infrared-camera/nircam-operations/nircam-dithers-and-mosaics>`_
     can be found here.
   * The ``fillval`` parameter specifies what value to use in the ouptut
     resampled image for any pixel which has no valid contribution from any
     input exposure.  The default value of ``INDEF`` indicates that the value
     from the last exposure will be used, while a value of 0 would result in
     holes.
   * The resampling can be controlled with the ``pixfrac``, ``kernel`` and
     ``weight_type`` parameters.
   * The ``pixfrac`` indicates the fraction by
     which input pixels are "shrunk" before being drizzled onto the
     output image grid, given as a real number between 0 and 1. This specifies
     the size of the footprint, or "dropsize", of a pixel in units of the input
     pixel size.
   * The ``kernel`` specifies the form of the kernel function used to distribute flux onto
     the separate output images.
   * The ``weight_type`` indicates the type of weighting image to apply with the bad pixel mask.
     Available options are ``ivm`` (default) for computing and using an inverse-variance map
     and ``exptime`` for weighting by the exposure time.
   * The ``good_bits`` parameter specifies what DQ values from the input exposure
     should be used when resampling to create the output mosaic.  Any pixel with a
     DQ value not included in this value (or list of values) will be ignored when
     resampling.
   * Resampled images will be written out to disk with the suffix ``_<asn_id>_outlier_i2d.fits``
     if the input model container has an <asn_id>, otherwise the suffix will be ``_outlier_i2d.fits``
     by default.
   * **If resampling is turned off** through the use of the ``resample_data`` parameter,
     a copy of the unrectified input images (as a ModelContainer)
     will be used for subsequent processing.

#. Create a median image from all grouped observation mosaics.

   * The median image is created by combining all grouped mosaic images or
     non-resampled input data (as planes in a ModelContainer) pixel-by-pixel.
   * The ``maskpt`` parameter sets the percentage of the weight image values to
     use, and any pixel with a weight below this value gets flagged as "bad" and
     ignored when resampled.
   * The median image is written out to disk as `_<asn_id>_median.fits` by default.

#. By default, the median image is blotted back (inverse of resampling) to
   match each original input image.

   * **If resampling is turned off**, the median image is compared directly to
     each input image.

#. Perform statistical comparison between blotted image and original image to identify outliers.

   * If resampling is disabled (``resample_data == False``) a large number of parameters
     will be ignored and instead the outlier mask will be computed using the following
     formula:

       .. math:: | image\_input - image\_median | > SNR * input\_err

   * When resampling is enabled, the comparison uses the original input images, the blotted
     median image, and the derivative of the blotted image to
     create a cosmic ray mask for each input image.
   * The derivative of the blotted image gets created using the blotted
     median image to compute the absolute value of the difference between each pixel and
     its four surrounding neighbors with the largest value being the recorded derivative.
   * These derivative images are used to flag cosmic rays
     and other blemishes, such as satellite trails. Where the difference is larger
     than can be explained by noise statistics, the flattening effect of taking the
     median, or an error in the shift (the latter two effects are estimated using
     the image derivative), the suspect pixel is masked.
   * The ``backg`` parameter specifies a user-provided value to be used as the
     background estimate.  This gets added to the background-subtracted
     blotted image to attempt to match the original background levels of the
     original input mosaic so that cosmic-rays (bad pixels) from the input
     mosaic can be identified more easily as outliers compared to the blotted
     mosaic.
   * Cosmic rays are flagged using the following rule:

     .. math:: | image\_input - image\_blotted | > scale*image\_deriv + SNR*noise

   * The ``scale`` is defined as the multiplicative factor applied to the
     derivative which is used to determine if the difference between the data
     image and the blotted image is large enough to require masking.
   * The ``noise`` is calculated using a combination of the detector read
     noise and the poisson noise of the blotted median image plus the sky background.
   * The user must specify two cut-off signal-to-noise values using the
     ``snr`` parameter for determining whether a pixel should be masked:
     the first for detecting the primary cosmic ray, and the second for masking
     lower-level bad pixels adjacent to those found in the first pass. Since
     cosmic rays often extend across several pixels, the adjacent pixels make
     use of a slightly lower SNR threshold.

#. Update input data model DQ arrays with mask of detected outliers.

Memory Model for Outlier Detection Algorithm
---------------------------------------------
The outlier detection algorithm can end up using massive amounts of memory
depending on the number of inputs, the size of each input, and the size of the
final output product.  Specifically,

#. The input :py:class:`~jwst.datamodels.ModelContainer` or
   :py:class:`~jwst.datamodels.CubeModel`
   for IFU data, by default, all input exposures would have been kept open in memory to make
   processing more efficient.

#. The initial resample step creates an output product for EACH input that is the
   same size as the final
   output product, which for imaging modes can span all chips in the detector while
   also accounting for all dithers.  For some Level 3 products, each resampled image can
   be on the order of 2Gb or more.

#. The median combination step then needs to have all pixels at the same position on
   the sky in memory in order to perform the median computation.  The simplest implementation
   for this step requires keeping all resampled outputs fully in memory at the same time.

Many Level 3 products only include a modest number of input exposures that can be
processed using less than 32Gb of memory at a time.  However, there are a number of
ways this memory limit can be exceeded.  This has been addressed by implementing an
overall memory model for the outlier detection that includes options to minimize the
memory usage at the expense of file I/O.  The control over this memory model happens
with the use of the ``in_memory`` parameter.  The full impact of this parameter
during processing includes:

#. The ``in_memory`` parameter gets passed to the :py:class:`~jwst.resample.ResampleStep`
   to set whether or not to keep the resampled images in memory or not.  By default,
   the outlier detection processing sets this parameter to `False` so that each resampled
   image gets written out to disk.

#. Computing the median image works section-by-section by only keeping 1Mb of each input
   in memory at a time.  As a result, only the final output product array for the final
   median image along with a stack of 1Mb image sections are kept in memory.

#. The final resampling step also avoids keeping all inputs in memory by only reading
   each input into memory 1 at a time as it gets resampled onto the final output product.

These changes result in a minimum amount of memory usage during processing at the obvious
expense of reading and writing the products from disk.


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

``--allowed_memory``
  Specifies the fractional amount of
  free memory to allow when creating the resampled image. If ``None``, the
  environment variable ``DMODEL_ALLOWED_MEMORY`` is used. If not defined, no
  check is made. If the resampled image would be larger than specified, an
  ``OutputTooLargeError`` exception will be generated.

  For example, if set to ``0.5``, only resampled images that use less than half
  the available memory can be created.

``--in_memory``
  Specifies whether or not to load and create all images that are used during
  processing into memory. If ``False``, input files are loaded from disk when
  needed and all intermediate files are stored on disk, rather than in memory.


Reference Files
===============

The ``outlier_detection_imaging`` step uses the PARS-OUTLIERDETECTIONIMAGINGSTEP parameter reference file.

.. include:: ../references_general/pars-outlierdetectionimagingstep_reffile.inc

.. automodapi:: jwst.outlier_detection.outlier_detection_imaging_step
