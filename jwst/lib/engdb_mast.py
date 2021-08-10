"""
Access the JWST Engineering Mnemonic Database through MAST
"""
import logging
from os import getenv
import requests

from astropy.table import Table
from astropy.time import Time

# Default MAST info.
MAST_BASE_URL = 'https://mast.stsci.edu'
API_URI = 'api/v0.1/Download/file'
SERVICE_URI = 'mast:jwstedb/'

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class EngdbMast():
    """
    Access the JWST Engineering Database through MAST

    Parameters
    ----------
    base_url : str
        The base url for the engineering RESTful service. If not defined,
        the environmental variable ENG_BASE_URL is queried. Otherwise
        the default MAST website is used.

    token : str
        The MAST access token. If not defined, the environmental variable
        MAST_API_TOKEN is queried. A token is required.
        For more information, see 'https://auth.mast.stsci.edu/'

    db_kwargs : dict
        Other keyword arguments that may be needed by other versions of the service
        but are not used by this service.

    Attributes
    ----------
    req: `requests.Request`
        The pre-built Request, with authorization and url defined.

    response: `requests.response`
        The results of the last query.

    session: `requests.Session`
        The session used for the connection.

    starttime: `astropy.time.Time`
        The start time of the last query.

    endtime: `astropy.time.Time`

    Raises
    ------
    RuntimeError
        Any and all failures with connecting with the MAST server.
    """
    def __init__(self, base_url=None, token=None, **db_kwargs):

        # Determine the database to use
        if base_url is None:
            base_url = getenv('ENG_BASE_URL', MAST_BASE_URL)
        if base_url[-1] != '/':
            base_url += '/'

        # Get the token
        if token is None:
            token = getenv('MAST_API_TOKEN', None)
        if token is None:
            raise RuntimeError('No MAST token provided but is required. See https://auth.mast.stsci.edu/ for more information.')

        # Check for basic aliveness.
        try:
            resp = requests.get(base_url + 'api/')
        except requests.exceptions.ConnectionError as exception:
            raise RuntimeError(f'MAST url: {base_url} is unreachable.') from exception
        if resp.status_code != 200:
            raise RuntimeError(f'MAST url: {base_url} is not available. Returned HTTPS status {resp.status_code}')

        # Basics are covered. Finalize initialization.
        self.req = requests.Request(method='GET',
                                    url=base_url + API_URI,
                                    headers={'Authorization': f'token {token}'})
        self.session = requests.Session()

    def get_records(
            self,
            mnemonic,
            starttime,
            endtime,
            time_format=None,
            **other_kwargs
    ):
        """
        Retrieve all results for a mnemonic in the requested time range.

        Parameters
        ----------
        mnemonic: str
            The engineering mnemonic to retrieve

        starttime: str or astropy.time.Time
            The, inclusive, start time to retireve from.

        endttime: str or astropy.time.Time
            The, inclusive, end time to retireve from.

        result_format: str
            The format to request from the service.
            If None, the `default_format` is used.

        time_format: str
            The format of the input time used if the input times
            are strings. If None, a guess is made.

        other_kwargs : dict
            Keyword arguments not relevant to this implementation.

        Returns
        -------
        records: `astropy.Table`
            Returns the resulting table.

        Notes
        -----
        The engineering service always returns the bracketing entries
        before and after the requested time range.
        """
        if not isinstance(starttime, Time):
            starttime = Time(starttime, format=time_format)
        if not isinstance(endtime, Time):
            endtime = Time(endtime, format=time_format)

        # Make the request
        mnemonic = mnemonic.strip()
        mnemonic = mnemonic.upper()
        starttime_fmt = starttime.strftime('%Y%m%dT%H%M%S')
        endtime_fmt = endtime.strftime('%Y%m%dT%H%M%S')
        uri = f'{mnemonic}-{starttime_fmt}-{endtime_fmt}.csv'
        self.req.params = {'uri': SERVICE_URI + uri}
        prepped = self.session.prepare_request(self.req)
        settings = self.session.merge_environment_settings(prepped.url, {}, None, None, None)
        response = self.session.send(prepped, **settings)

        # Convert to table.
        r_list = response.text.split('\r\n')
        table = Table.read(r_list, format='ascii.csv')

        return table
