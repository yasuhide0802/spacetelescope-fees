from jwst import datamodels
from jwst.white_light.white_light import white_light

from astropy.time import Time, TimeDelta
import numpy as np
import pytest

@pytest.fixture(scope='module')
def make_datamodel():
    """Make data for white light tests"""

    model = datamodels.MultiSpecModel()
    model.meta.exposure.group_time = 0.15904
    model.meta.exposure.ngroups = 60
    model.meta.exposure.start_time = 58627.0
    model.meta.exposure.integration_start = 1
    model.meta.exposure.integration_end = 2

    # Make data arrays
    flux = np.random.rand(200)*1e-9
    wavelength = np.arange(11,13, step=0.01)
    error = np.random.rand(200)
    surf_bright = np.zeros(200)
    sb_error = np.zeros(200)
    dq = np.ones(200)
    background = np.zeros(200)
    berror = np.zeros(200)
    npixels = np.zeros(200)

    data = [(i,j,k,l,m,n,o,p,q) for i,j,k,l,m,n,o,p,q in zip(wavelength, flux, error,
                                                            surf_bright, sb_error, dq,
                                                            background, berror, npixels)]

    spec_table = np.array(data, dtype=[('WAVELENGTH', 'f8'),('FLUX', 'f8'),
                                    ('ERROR', 'f8'), ('SURF_BRIGHT', 'f8'),
                                    ('SB_ERROR', 'f8'),('DQ', 'u4'), ('BACKGROUND', 'f8'),
                                    ('BERROR', 'f8'), ('NPIXELS', 'f8')])



    spectrum = datamodels.SpecModel()
    spectrum.spec_table = spec_table

    model.spec.append(spectrum)

    integrations = [(  1, 58627.53891071, 58627.53896565, 58627.5390206 , 0., 0., 0.),
                    (  2, 58627.5390206 , 58627.53907555, 58627.5391305 , 0., 0., 0.),
                    (  3, 58627.5391305 , 58627.53918544, 58627.53924039, 0., 0., 0.)]

    integration_table = np.array(integrations, dtype=[('integration_number', 'i4'),
                                                      ('int_start_MJD_UTC', 'f8'),
                                                      ('int_mid_MJD_UTC', 'f8'),
                                                      ('int_end_MJD_UTC', 'f8'),
                                                      ('int_start_BJD_TDB', 'f8'),
                                                      ('int_mid_BJD_TDB', 'f8'),
                                                      ('int_end_BJD_TDB', 'f8')])
    model.int_times = (integration_table)

    return model


def test_white_light_with_int_tables(make_datamodel):
    data = make_datamodel
    result = white_light(data,None,None)

    int_start = data.meta.exposure.integration_start

    # We know there is only one table, so set we are hardcoding.
    ntables = 1
    int_num = data.int_times['integration_number']
    mid_utc = data.int_times['int_mid_MJD_UTC']

    offset = int_start - int_num[0]
    time_arr = np.zeros(ntables, dtype=np.float64)
    time_arr[0 : 1] = mid_utc[offset : offset + ntables]
    int_times = Time(time_arr, format='mjd', scale='utc')

    # Sum the fluxes
    fluxsums = data.spec[0].spec_table['FLUX'].sum()

    assert result['MJD'] == int_times.mjd
    assert result['whitelight_flux'] ==  fluxsums

def test_white_light_with_expstart(make_datamodel):
    data = make_datamodel

    # Make the integration_end larger than the number of
    # integration. This forces the algorithm to use EXPSTART
    # and TGROUP
    data.meta.exposure.integration_end = 4

    result = white_light(data,None,None)

    dt_arr = np.zeros(1, dtype=np.float64)
    dt = (data.meta.exposure.group_time *
         (data.meta.exposure.ngroups + 1))

    # We know there is only one table, so set we are hardcoding.
    ntables_current = 1
    dt_arr[0 : 1] = np.arange(1, 1 + ntables_current) * dt - (dt / 2.)
    int_dt = TimeDelta(dt_arr, format='sec')

    int_times = (Time(data.meta.exposure.start_time, format='mjd')
                 + int_dt)
    # Sum the fluxes
    fluxsums = data.spec[0].spec_table['FLUX'].sum()

    assert result['MJD'][0] == int_times.mjd[0]
    assert result['whitelight_flux'] == fluxsums
