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
    help="""TSTART and TSTOP for selecting the events. supercedes gtifits \n
    Example: -gti TSTART TSTOP""")

parser.add_argument(
    '-energy', nargs='+', type=float, 
    help="""enmin and enmax for selecting the events. \n
    Example: -energy enmin enmax""")

parser.add_argument(
    '-gtifits', type=int, 
    help="""specify the HDU number of the GTI table if it's part of the event list.""")

parser.add_argument(
    '-nt', type=int, default=1, 
    help="""Telescope number.""")

parser.add_argument('-verbose', type=bool, required=False, default=True,
    help='Add a time constraint to the attitude file. Only the periods within the specified time would be considered.')

parser.add_argument('-overwrite', type=bool, required=False, default=True,
    help='Overwrite if set as True.')

args = parser.parse_args()

verbose = args.verbose
vprint = verboseprint(verbose)



evt = args.evt
out = args.out
if args.energy is not None:
    enmin, enmax = args.energy
overwrite = args.overwrite
    
ntele = str(args.nt)
arthdu = fits.open(evt)
evt = Table(arthdu[1].data)

if args.gtifits is not None:
    tstart = arthdu[args.gtifits].data['START']
    tstop = arthdu[args.gtifits].data['STOP']
if args.gti is not None:
    tstart, tstop = args.gti

if args.gti is None and args.gtifits is None:
    tstop = evt['TIME'][-1]
    tstart = evt['TIME'][0]


# filtering events outside of GTI

if len(tstart) > 1:
    i = 0
    ix = (evt['TIME'] >= tstart[i]) & (evt['TIME'] <= tstop[i])
    exposure = tstop[i] - tstart[i]
    evtout = evt[ix]
    for i in np.arange(len(tstart) - 1) + 1:
        ix = (evt['TIME'] >= tstart[i]) & (evt['TIME'] <= tstop[i])
        evtout = np.hstack((evtout, evt[ix]))
        exposure += tstop[i] - tstart[i]
    evt = evtout
else:
    ix = (evt['TIME'] <= tstop) & (evt['TIME']>= tstart)
    evt = evt[ix]
    exposure = tstop - tstart
if args.energy is not None:
    ix = (evt['ENERGY'] <= enmax) & (evt['ENERGY']>= enmin)
    evt = evt[ix]

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
    arthdu[i].header['MJDREFF'] = 51543.875
    arthdu[i].header['TIMEREF'] = 'LOCAL'
    arthdu[i].header['EQUINOX'] = 2000
    arthdu[i].header['EXPOSURE'] = exposure
    arthdu[i].header['LIVETIME'] = exposure
arthdu[1].data = fits.table_to_hdu(evt).data
arthdu.writeto(out, overwrite=overwrite)
