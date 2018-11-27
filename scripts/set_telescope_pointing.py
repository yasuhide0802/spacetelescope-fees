#!/usr/bin/env python

# Copyright (C) 2010-2011 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

'''
This script adds absolute pointing information to the FITS files provided
to it on the command line (one or more).

Currently it only uses a constant value for the engineering keywords
since the Engineering Database does not yet contain them.

It assumes the following keywords are present in the file header:

V2_REF (arcseconds)
V3_REF (arcseconds)
VPARITY (+1 or -1)
V3I_YANG (decimal degrees)

The keywords added are:

RA_V1
DEC_V1
PA_V3
CRVAL1
CRVAL2
PC1_1
PC1_2
PC2_1
PC2_2

It does not currently place the new keywords in any particular location
in the header other than what is required by the standard.
'''
import argparse
import logging
import sys

from jwst.lib.set_telescope_pointing import add_wcs

logger = logging.getLogger('jwst')
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Update basic WCS information in JWST exposures from the engineering database.'
    )
    parser.add_argument(
        'exposures', type=str, nargs='+',
        help='List of JWST exposures to update.'
    )
    parser.add_argument(
        '--strict-time', action='store_true',
        help='Force pointings to be within the observation time.'
    )
    args = parser.parse_args()

    for filename in args.exposures:
        logger.info('Setting pointing for {}'.format(filename))
        add_wcs(filename, strict_time=args.strict_time)
