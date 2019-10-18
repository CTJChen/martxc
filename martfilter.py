#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package (internal).
This script filters an eventlist and fill in necessary keywords. 
"""


import sys
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
from martxclib.martxcfun import *

parser = Parser(
    description=__doc__,
    epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


# Required arguments
parser.add_argument(
    '-evt', type=str, required=True,
    help='Name of the event file.')

parser.add_argument(
    '-out', type=str, required=False, default='expmap.fits',
    help='Name of the output file.')

parser.add_argument(
    '-gti', nargs='+', type=float, 
    help="""TSTART and TSTOP for selecting the events. \n
    Example: -gti TSTART TSTOP""")

parser.add_argument(
    '-nt', type=int, default=1, 
    help="""Telescope number.""")

parser.add_argument('-verbose', type=bool, required=False, default=True,
    help='Add a time constraint to the attitude file. Only the periods within the specified time would be considered.')

args = parser.parse_args()

verbose = args.verbose
vprint = verboseprint(verbose)



evt = args.evt
out = args.out
if args.gti is not None:
    tstart, tstop = args.gti
ntele = str(args.nt)
arthdu = fits.open(evt)
evt = Table(arthdu[1].data)

# filter - GTI should be the following (based on RA/DEC distributions)

if args.gti is not None:
    ix = np.where((evt['TIME'] <= tstop) & (evt['TIME']>= tstart))[0]
    evt = evt[ix]
else:
    tstop = evt['TIME'][-1]
    tstart = evt['TIME'][0]
igoodevt = (evt['RAW_X'] >= 0) & (evt['RAW_X'] <= 47) & \
(evt['RAW_Y'] >= 0) & (evt['RAW_Y'] <= 47) & \
(evt['PHA_BOT'] >= 0.) & (evt['PHA_BOT'] <= 1023)
gtifrac = sum(igoodevt) / float(len(evt))
print('fraction of good eventx:', sum(igoodevt)/len(evt))
evt = evt[igoodevt]
for i in [0,1]:
    arthdu[i].header['DETNAME'] = 'M' + ntele
    arthdu[i].header['TELESCOP'] = 'SRG'
    arthdu[i].header['INSTRUME'] = 'ART-XC'
    arthdu[i].header['RADECSYS'] = 'FK5'
    arthdu[i].header['TIMESYS'] = 'TT' 
    arthdu[i].header['MJDREFF'] = 0.0
    arthdu[i].header['TIMEREF'] = 'LOCAL'
    arthdu[i].header['EQUINOX'] = 2000
    arthdu[i].header['EXPOSURE'] = tstop - tstart
    arthdu[i].header['LIVETIME'] = tstop - tstart
arthdu[1].data = fits.table_to_hdu(evt).data
arthdu.writeto(out, overwrite=True)
