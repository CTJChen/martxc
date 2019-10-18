#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package (internal).
This script takes the level 2 event file
and assign a PI for each event based on the given  an eventlist and fill in necessary keywords. 
"""


import sys
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table as tab
from astropy.table import Column
from martxclib.martxcfun import *

def get_pi(energy):
    if type(energy) is np.ndarray:
        pivalue = np.zeros_like(energy,dtype=int)
        for idn, en in enumerate(energy):
            pivalue[idn] = np.abs(enarr - en).argmin()
            if enarr[pivalue[idn]] > en:
                pivalue[idn] -= 1
        return pivalue
    else:
        pivalue = np.abs(enarr - energy).argmin()
        if enarr[pivalue] > energy:pivalue -= 1
        return pivalue
        

def add_column(hdu, coldata, colname, formatstr, unit):
    """
    Add a column to the given hdu
    coldata should be a numpy array
    Return:
        A fits hdu object with the new column.
    """
    table = hdu.data
    col = fits.Column(name=colname, format=formatstr, unit=unit, array=coldata)
    # NOTE: append the new time column to the *last*!
    # Otherwise the TLMIN??/TLMAX?? keyword pairs, which record the
    # minimum/maximum values of corresponding columns, will become
    # *out of order*. Therefore the output FITS file causes weird problems
    # with DS9 and DM tools.
    newtable = fits.BinTableHDU.from_columns(
            table.columns + fits.ColDefs([col]))
    hdu.data = newtable.data
    # update header
    hdu.header.update(newtable.header)
    return hdu 



parser = Parser(
    description=__doc__,
    epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


# Required arguments
parser.add_argument(
    '-evt', type=str, required=True,
    help='Name of the event file.')

parser.add_argument(
    '-arf', type=str, required=True,
    help='Name of the ARF file.')

parser.add_argument(
    '-out', type=str, required=False, default='pi.fits',
    help='Name of the output file.')

parser.add_argument('-verbose', type=bool, required=False, default=True,
    help='set for verbosity')

parser.add_argument('-overwrite', type=bool, required=False, default=True,
    help='Overwrite if set as True.')

args = parser.parse_args()

verbose = args.verbose
vprint = verboseprint(verbose)


evt = args.evt
arf = args.arf
out = args.out
overwrite = args.overwrite


arfhdu = fits.open(arf)
arf1 = tab(arfhdu[1].data)

enarr = arf1['ENERG_LO'].quantity.value
enarr = np.hstack([enarr, np.array(arf1['ENERG_HI'][-1])])

arthdu = fits.open(evt)
piarr = get_pi(arthdu[1].data['ENERGY'])

evthead = arthdu[1].header

arthdu[1] = add_column(arthdu[1], piarr, 'PI', 'I', '')
arthdu.writeto(out, overwrite=overwrite) 
print('wrote PI arrays into ' + out)

