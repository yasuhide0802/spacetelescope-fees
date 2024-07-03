"""General utility objects"""

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags
import numpy as np


def use_datamodel(input, model_class=None):
    """Determine if input is a datamodel, if so return it, else open it

    Parameters
    ----------
    input : string or datamodel
        Either the name of the file to open or a datamodel

    model_class : jwst datamodel class
        Specific class desired for the output datamodel, e.g. RampModel

    Returns
    -------
    model : datamodel
        The datamodel object
    """
    if isinstance(input, datamodels.JwstDataModel) or isinstance(input, datamodels.RampModel):
        model = input
    else:
        model = datamodels.open(input)
    if model_class is not None:
        model = model_class(model)
    return model

def copy_datamodel(input, modify_input):
    """Return a copy of the datamodel and set the input to None to recover memory or simply
    return the input if the data is to be modified in-place.

    Parameters
    ----------
    input : jwst datamodel
        The datamodel to copy

    modify_input = boolean
        If True modify the input datamodel in-place, don't make a copy

    Returns
    -------
    input : datamodel
        The original datamodel object

    model_copy : datamodel
        The datamodel copy
    """
    if not modify_input:
        model_copy = input.copy()
        return model_copy, None
    else:
        return input, None


def set_nans_to_donotuse(data, dq):
    """Set all NaN values in the data that have an even value to
    DO_NOT_USE.

    Parameters
    ----------
    data : numpy array
        The science data array to find NaN values and
        check of these have a DQ flag=DO_NOT_USE, or
        set it if not.

    dq : numpy array
        The DQ array to be checked.

    Returns
    -------
    dq : numpy array
        The updated DQ array.
    """
    dq[np.isnan(data)] |= dqflags.pixel['DO_NOT_USE']
    return dq


class LoggingContext:
    """Logging context manager

    Keep logging configuration within a context

    Based on the Python 3 Logging Cookbook example

    Parameters
    ==========
    logger: logging.Logger
        The logger to modify.

    level: int
        The log level to set.

    handler: logging.Handler
        The handler to use.

    close: bool
        Close the handler when done.
    """

    def __init__(self, logger, level=None, handler=None, close=True):
        self.logger = logger
        self.level = level
        self.handler = handler
        self.close = close

        self.old_level = None

    def __enter__(self):
        if self.level is not None:
            self.old_level = self.logger.level
            self.logger.setLevel(self.level)
        if self.handler:
            self.logger.addHandler(self.handler)

    def __exit__(self, et, ev, tb):
        if self.level is not None:
            self.logger.setLevel(self.old_level)
        if self.handler:
            self.logger.removeHandler(self.handler)
        if self.handler and self.close:
            self.handler.close()
        # implicit return of None => don't swallow exceptions
