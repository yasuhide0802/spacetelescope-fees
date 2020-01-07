from .image import ImageModel
from .model_base import DataModel


__all__ = ['IFUImageModel']


class IFUImageModel(DataModel):
    """
    A data model for 2D IFU images.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array

    zeroframe : numpy float32 array
         Zeroframe array

    area : numpy float32 array
         Pixel area map array

    var_poisson : numpy float32 array
         variance due to poisson noise

    var_rnoise : numpy float32 array
         variance due to read noise

    wavelength : numpy float32 array
         wavelength

    pathloss : numpy float32 array
         pathloss correction
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/ifuimage.schema"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, ImageModel):
            super(IFUImageModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            if init.hasattr('area'):
                self.area = init.area
            if init.hasattr('var_poisson'):
                self.var_poisson = init.var_poisson
            if init.hasattr('var_rnoise'):
                self.var_rnoise = init.var_rnoise
            return

        super(IFUImageModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
