.. _outlier-detection-coron:

Outlier Detection for Coronagraphic Data
========================================

:Classes: `jwst.outlier_detection.OutlierDetectionCoronStep`
:Aliases: outlier_detection_coron

This module serves as the interface for applying outlier detection to coronagraphic
image observations.

Specifically, this routine performs the following operations:

#. Extract parameter settings from input model and merge them with any user-provided values.

#. Convert input data, as needed, to make sure it is in a format that can be processed.
   A :py:class:`~jwst.datamodels.CubeModel` serves as the basic format for all processing
   performed by this step.

#. Create a median image preserving the spatial dimensions of the cube.

   * The ``maskpt`` parameter sets the percentage of the weight image values to
     use, and any pixel with a weight below this value gets flagged as "bad" and
     ignored when resampled.

#. Perform statistical comparison between median image and original image to identify outliers.

   The core detection algorithm uses the following to generate an outlier mask

   .. math:: | image\_input - image\_median | > SNR*input\_err

#. Update input data model DQ arrays with mask of detected outliers.


Step Arguments
==============

``--maskpt``
  The percent of maximum weight to use as lower-limit for valid data;
  valid values go from 0.0 to 1.0.

``--snr``
  The signal-to-noise values to use for bad pixel identification.
  Valid values are any positive float.

``--save_intermediate_results``
  Specifies whether or not to save any intermediate products created
  during step processing.

``--good_bits``
  The DQ bit values from the input image DQ arrays
  that should be considered 'good' when building the weight mask. See
  DQ flag :ref:`dq_parameter_specification` for details.


Reference Files
===============

The ``outlier_detection_coron`` step uses the PARS-OUTLIERDETECTIONCORONSTEP parameter reference file.

.. include:: ../references_general/pars-outlierdetectioncoronstep_reffile.inc

.. automodapi:: jwst.outlier_detection.outlier_detection_coron_step
