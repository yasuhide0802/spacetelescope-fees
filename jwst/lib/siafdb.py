"""SIAF Database Access

Provide a common interface to different versions of the SIAF.

Under operations, the SIAF is found in a sqlite database.
Otherwise, use the standard interface defined by the `pysiaf` package
"""
from collections import namedtuple
from datetime import date
import logging
import os
from pathlib import Path

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Attempt import of sqlite3 and flag whether there are fails.
try:
    import sqlite3
except ModuleNotFoundError:
    sqlite3 = False

# Map instrument three character mnemonic to full name
INSTRUMENT_MAP = {
    'fgs': 'fgs',
    'mir': 'miri',
    'nis': 'niriss',
    'nrc': 'nircam',
    'nrs': 'nirspec'
}

# SIAF container
# The names should correspond to the names in the ``wcsinfo`` schema.
# It is populated by the SIAF values in the PRD database based
# on APERNAME and UseAfterDate and used to populate the keywords
# in Level1bModel data models.
SIAF = namedtuple("SIAF", ["v2_ref", "v3_ref", "v3yangle", "vparity",
                           "crpix1", "crpix2", "cdelt1", "cdelt2",
                           "vertices_idl"])
# Set default values for the SIAF.
# Values which are needed by the pipeline are set to None which
# triggers a ValueError if missing in the SIAF database.
# Quantities not used by the pipeline get a default value -
# FITS keywords and aperture vertices.
SIAF.__new__.__defaults__ = (None, None, None, None, 0, 0, 3600, 3600,
                             (0, 1, 1, 0, 0, 0, 1, 1))

SIAF_REQUIRED = ['V2Ref', 'V3Ref', 'V3IdlYAngle', 'VIdlParity']
SIAF_OPTIONAL = ['XSciRef', 'YSciRef', 'XSciScale', 'YSciScale']
SIAF_VERTICIES = ['XIdlVert1', 'XIdlVert2', 'XIdlVert3', 'XIdlVert4',
                  'YIdlVert1', 'YIdlVert2', 'YIdlVert3', 'YIdlVert4']
SIAF_MAP = {'V2Ref': 'v2_ref', 'V3Ref': 'v3_ref', 'V3IdlYAngle': 'v3yangle', 'VIdlParity': 'vparity',
            'XSciRef': 'crpix1', 'YSciRef': 'crpix2', 'XSciScale': 'cdelt1', 'YSciScale': 'cdelt2'}


class SiafDb:
    """SIAF Database Access

    Provide a common interface to different versions of the SIAF.

    Under operations, the SIAF is found in a sqlite database.
    Otherwise, use the standard interface defined by the `pysiaf` package

    Parameters
    ----------
    source : None, str, or a path-like object
        The SIAF database source. See notes for more details.

    Notes
    -----
    The interpretation of `source` is as follows:

    If None, the environmental 'XML_DATA' is queried for a value. If the folder
    'SIAFXML' exists, `pysiaf` is used with that folder. If not, the 'prd.db'
    file is tried. If 'XML_DATA' is not defined, then the `pysiaf` package is
    used.

    Otherwise `source` will be treated as a path to some SIAF database.
    `SiafDbSqlite` will first be attempted, then `SiafDbPysiaf`.

    """
    def __init__(self, source=None):
        self._db = None
        self._source = None

        # If no source, retrieve the environmental XML_DATA
        if source is None:
            source = os.environ.get('XML_DATA', None)
            if source is not None:
                source = Path(source)
                if (source / 'SIAFXML').exists():
                    source = source / 'SIAFXML'
                else:
                    source = source / 'prd.db'
        self._source = source

        # Attempt to access source as an sqlite database.
        try:
            db = SiafDbSqlite(source)
        except ValueError:
            # Source is incompatible.
            logger.debug('Could not open as a sqlite object: %s', source)
        else:
            self._db = db
            return

        # Attempt to access source as a pysiaf source
        self._db = SiafDbPySiaf(source)

    def close(self):
        self._db.close()

    def get_wcs(self, aperture, useafter=None):
        """
        Query the SIAF database file and get WCS values.

        Given an ``APERTURE_NAME`` and a ``USEAFTER`` date query the SIAF database
        and extract the following keywords:
        ``V2Ref``, ``V3Ref``, ``V3IdlYAngle``, ``VIdlParity``,
        ``XSciRef``, ``YSciRef``, ``XSciScale``, ``YSciScale``,
        ``XIdlVert1``, ``XIdlVert2``, ``XIdlVert3``, ``XIdlVert4``,
        ``YIdlVert1``, ``YIdlVert2``, ``YIdlVert3``, ``YIdlVert4``

        Parameters
        ----------
        aperture_name : str
            The name of the aperture to retrieve.
        useafter : str
            The date of observation (``model.meta.date``).
            If None, set to the current day.

        Returns
        -------
        siaf : namedtuple
            The SIAF namedtuple with values from the PRD database.
        """
        if not useafter:
            useafter = date.today().strftime('%Y-%m-%d')
        return self._db.get_wcs(aperture, useafter)

    def __deepcopy__(self, memo):
        """Do not copy, just return."""
        return self

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self._db.close()


class SiafDbPySiaf:
    """Use pysiaf as the source of siaf information

    Parameters
    ----------
    source : None, str, or a file-like object
        The SIAF database source. See notes for more details.

    Notes
    -----
    The interpretation of `source` is as follows:

    If None, then the `pysiaf` package is used.
    If a string, the string is treated as a path.
    If that path is to a folder, the `pysiaf` package is used with the folder
        as the XML source folder. See the `pysiaf` package for more information.
    Otherwise, fail.
    """
    def __init__(self, source=None):
        try:
            import pysiaf
        except ImportError:
            raise ValueError('Package "pysiaf" is not installed. Cannot use the pysiaf api')
        self.pysiaf = pysiaf

        if source is not None:
            source = Path(source)
            if not source.is_dir():
                raise ValueError('Source %s: Needs to be a folder for use with pysiaf')
        self._source = source

    def close(self):
        pass

    def get_wcs(self, aperture, useafter):
        """
        Query the SIAF database file and get WCS values.

        Given an ``APERTURE_NAME`` and a ``USEAFTER`` date query the SIAF database
        and extract the following keywords:
        ``V2Ref``, ``V3Ref``, ``V3IdlYAngle``, ``VIdlParity``,
        ``XSciRef``, ``YSciRef``, ``XSciScale``, ``YSciScale``,
        ``XIdlVert1``, ``XIdlVert2``, ``XIdlVert3``, ``XIdlVert4``,
        ``YIdlVert1``, ``YIdlVert2``, ``YIdlVert3``, ``YIdlVert4``

        Parameters
        ----------
        aperture : str
            The name of the aperture to retrieve.
        useafter : str
            The date of observation (``model.meta.date``)

        Returns
        -------
        siaf : namedtuple
            The SIAF namedtuple with values from the PRD database.
        """
        instrument = INSTRUMENT_MAP[aperture[:3].lower()]
        siaf = self.pysiaf.Siaf(instrument, basepath=self._source)
        aperture = siaf[aperture.upper()]

        # Build the SIAF entry. Missing required values is an error.
        # Otherwise, use defaults.
        default_siaf = SIAF()
        values = {SIAF_MAP[key]: getattr(aperture, key) for key in SIAF_REQUIRED}
        if not all(values):
            raise RuntimeError(f'Required SIAF entries for {instrument}:{aperture} are not all defined: {values}')
        for key in SIAF_OPTIONAL:
            value = getattr(aperture, key)
            value = value if value else getattr(default_siaf, SIAF_MAP[key])
            values[SIAF_MAP[key]] = value
        vertices = list()
        for key in SIAF_VERTICIES:
            value = getattr(aperture, key)
            value = value if value else getattr(default_siaf, SIAF_MAP[key])
            vertices.append(value)
        vertices = tuple(vertices)

        # Fill out the Siaf
        siaf = SIAF(**values, vertices_idl=vertices)

        return siaf


class SiafDbSqlite:
    """Use a sqlite db as the source of siaf information

    Parameters
    ----------
    source : str, or a file-like object
        The SIAF database source.
    """
    def __init__(self, source):
        self._cursor = None
        self._db = None
        self._source = None

        if not sqlite3:
            raise ValueError('Package "sqlite3" is not install. Cannot use the sqlite api.')

        if source is None:
            raise ValueError('Source: %s: Cannot be None. A sqlite path needs to be specified.', source)

        source = Path(source)
        if not source.exists():
            raise ValueError('Source: %s does not exist.', source)
        self._source = f'file:{str(source)}?mode=ro'

        self._db = sqlite3.connect(self._source, uri=True)
        self._cursor = self._db.cursor()
        logger.info("Using SIAF database from %s", source)

    def close(self):
        self._db.close()

    def get_wcs(self, aperture, useafter):
        """
        Query the SIAF database file and get WCS values.

        Given an ``APERTURE_NAME`` and a ``USEAFTER`` date query the SIAF database
        and extract the following keywords:
        ``V2Ref``, ``V3Ref``, ``V3IdlYAngle``, ``VIdlParity``,
        ``XSciRef``, ``YSciRef``, ``XSciScale``, ``YSciScale``,
        ``XIdlVert1``, ``XIdlVert2``, ``XIdlVert3``, ``XIdlVert4``,
        ``YIdlVert1``, ``YIdlVert2``, ``YIdlVert3``, ``YIdlVert4``

        Parameters
        ----------
        aperture : str
            The name of the aperture to retrieve.
        useafter : str
            The date of observation (``model.meta.date``)

        Returns
        -------
        siaf : namedtuple
            The SIAF namedtuple with values from the PRD database.
        """
        logger.info("Querying SIAF for aperture "
                    "%s with USEAFTER %s", aperture, useafter)
        RESULT = {}
        try:
            self._cursor.execute("SELECT Apername, V2Ref, V3Ref, V3IdlYAngle, VIdlParity, "
                                 "XSciRef, YSciRef, XSciScale, YSciScale, "
                                 "XIdlVert1, XIdlVert2, XIdlVert3, XIdlVert4, "
                                 "YIdlVert1, YIdlVert2, YIdlVert3, YIdlVert4 "
                                 "FROM Aperture WHERE Apername = ? "
                                 "and UseAfterDate <= ? ORDER BY UseAfterDate LIMIT 1",
                                 (aperture, useafter))
            for row in self._cursor:
                RESULT[row[0]] = tuple(row[1:17])
            self._db.commit()
        except (sqlite3.Error, sqlite3.OperationalError) as err:
            print("Error: " + err.args[0])
            raise

        logger.info("loaded %s table rows", len(RESULT))
        default_siaf = SIAF()
        if RESULT:
            # This populates the SIAF tuple with the values from the database.
            # The last 8 values returned from the database are the vertices.
            # They are wrapped in a list and assigned to SIAF.vertices_idl.
            values = list(RESULT.values())[0]
            vert = values[-8:]
            values = list(values[: - 8])
            values.append(vert)
            # If any of "crpix1", "crpix2", "cdelt1", "cdelt2", "vertices_idl" is None
            # reset ot to the default value.
            for i in range(4, 8):
                if values[i] is None:
                    values[i] = default_siaf[i]
            siaf = SIAF(*values)
            return siaf
        else:
            raise RuntimeError(f'No SIAF entries found for {aperture}')

    def __del__(self):
        """Close out any connections"""
        if self._db:
            self._db.close()
