#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package. 
This is a script for internal usage only.
Current version : 0.0.1
What this script does: 
Takes an RA/DEC and tells you which tile it belongs to.
Also check if it's in NEP
"""
from astropy.io import ascii
import pandas as pd
import numpy as np
from astropy.table import Table as tab
import argparse
from martxclib.martxcfun import Parser

parser = Parser(
    description=__doc__,
    epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '-ra', type=float, required=True,
    help='RA in degrees')

parser.add_argument(
    '-dec', type=float, required=True,
    help='DEC in degrees')

args = parser.parse_args()
src_ra = args.ra
src_dec = args.dec

data = ascii.read('/home/ctchen/lib/martxc_local/artxc_survey_grid.txt')
tiles = data['tile'].quantity.value.astype(str)
for idx, tt in enumerate(tiles):
    tiles[idx] = tt.zfill(6)
data['tile'] = tiles
data = tab(data).to_pandas()

raarr = data['ramin'].values
raarr = np.hstack([raarr, data['ramax'].values[-1]])

decarr = data['decmin'].values
decarr = np.hstack([decarr, data['decmax'].values[-1]])


def isintile(row, ra, dec):
    if (row.ramin <= ra < row.ramax) & (row.decmin <= dec < row.decmax):
        return True
    else:
        return False    

if sum(data.apply(isintile, axis=1,args=(src_ra, src_dec))) == 1:
    tile = data[data.apply(isintile, axis=1,args=(src_ra, src_dec))].values[0][0]
    nep = data[data.apply(isintile, axis=1,args=(src_ra, src_dec))].values[0][-1]
else:
    tile = 'N/A'
    nep = 'N/A'
if nep == 1:
    boolstr = 'IS in NEP'
elif nep == 0:
    boolstr = 'IS NOT in NEP'
else:
    boolstr = 'not clear if it is in NEP'
print(src_ra, src_dec, 'is in tile ', tile, 'and ' + boolstr) 