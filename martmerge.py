#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package. 
This is a script for internal useage only.
Update event list by combining the existing event with an additional event
Only events with the same DETN and URDN can be merged. 
The output event list will be sorted 
GTI table and STDGTI table will also be merged if present. 
The merged event will be written into history
History entries, by default, will be checked before merging to avoid duplication.
Note -- assuming all events were sorted by increasing time
"""

import sys, os
import numpy as np
import astropy.io.fits as fits
from martxclib.martxcfun import *
from matplotlib.colors import LogNorm
from astropy.table import Table as tab
from astropy.table import vstack
#caldb = check_caldb()


# ART-XC has 48 by 48 pixels
# with a pixel angular size of 0.724 arcmin

parser = Parser(description=__doc__,
    epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required arguments
parser.add_argument('-evtm', type=str, required=False, 
    help='Name of the master event file.')

parser.add_argument('-evta', type=str, required=False, 
    help='Name of the event file to be appended.')

parser.add_argument('-evtlist', type=str, required=False, 
    help='Name of the txt file including all the files to be merged. supercede evtm and evta.')

parser.add_argument('-out', type=str, required=False,
    help='Name of the output file, if empty, -evtm will be overwritten')

parser.add_argument('-verbose', type=bool, required=False, default=False,
    help='Seeting verbosity.')

parser.add_argument('-overwrite', type=bool, required=False, default=True,
    help='overwrite if set as True. Default is True!')

parser.add_argument('-force', type=bool, required=False, default=False,
    help='force if set as True. Warning! Make sure you no what you need before setting force == True.')

args = parser.parse_args()


evtm = args.evtm
evta = args.evta
evtlist = args.evtlist
force = args.force
out = args.out
if out is None:
    out = evtm
    overwrite = True
verbose = args.verbose
vprint = verboseprint(verbose)
overwrite = args.overwrite

if evtlist is not None:
    i = 0
    evtfiles = open(evtlist,'r').readlines()
    evtm = evtfiles[i].strip()
    fname_m = os.path.abspath(evtm).split('/')[-2] + '/' + os.path.abspath(evtm).split('/')[-1]
    hdum = fits.open(evtm)
    keylist = []
    for key in hdum[0].header.keys():keylist.append(key)
    if not 'HISTORY' in keylist:
        vprint('No history in the primary HDU')
        hdum[0].header['HISTORY'] = 'Initializing with ' + fname_m
    for i in np.arange(len(evtfiles) - 1) + 1:
        evta = evtfiles[i].strip()
        fname_a = os.path.abspath(evta).split('/')[-2] + '/' + os.path.abspath(evta).split('/')[-1]

        for card in hdum[0].header['HISTORY']:
            if fname_a in card:
                print('no merging needed, ' + fname_a + ' has been merged before.')
        hdua = fits.open(evta)
        tabm = tab(hdum[1].data)
        taba = tab(hdua[1].data)
        if tabm['TIME'][-1] <= taba['TIME'][0]:
            # if all events the list to be appended are later than the last event of primary
            newtab = vstack([tabm, taba])
        elif taba['TIME'][-1] <= tabm['TIME'][0]:
            # if the first event of the primary list is later than all events in the list to be appended
            newtab = vstack([taba, tabm])
        elif force:
            # something is wrong, sort the output
            vprint('The time interval of evta is already part of evtm')
            vprint('force == True, continuing now.')
            newtab = vstack([tabm, taba])
            newtab.sort('TIME')
        else:
            # something is wrong, sort the output 
            sys.stderr.write('error: The time interval of ' + fname_a + ' is already part of evtm, stopping now')
        hdum[1].data = fits.table_to_hdu(newtab).data
        tabm = tab(hdum[3].data)
        taba = tab(hdua[3].data)
        newtab = vstack([tabm, taba])
        newtab.sort('START')
        hdum[3].data = fits.table_to_hdu(newtab).data    

        # update the STDGTI table
        tabm = tab(hdum[4].data)
        taba = tab(hdua[4].data)
        newtab = vstack([tabm, taba])
        newtab.sort('START')
        hdum[4].data = fits.table_to_hdu(newtab).data
        print('Appended ' + fname_a + ' to  ', fname_m)
        hdum[0].header['HISTORY'] = 'Appended ' + fname_a
    hdum.writeto(out, overwrite=overwrite)
else:
    fname_m = os.path.abspath(evtm).split('/')[-2] + '/' + os.path.abspath(evtm).split('/')[-1]
    fname_a = os.path.abspath(evta).split('/')[-2] + '/' + os.path.abspath(evta).split('/')[-1]
    hdum = fits.open(evtm)
    hdua = fits.open(evta)
    keylist = []
    for key in hdum[0].header.keys():keylist.append(key)
    if not 'HISTORY' in keylist:
        vprint('No history in the primary HDU')
        hdum[0].header['HISTORY'] = 'Initil event list ' + fname_m
    else:
        for card in hdum[0].header['HISTORY']:
            if fname_a in card:
                print('no merging needed, evta has been merged before.')
                sys.exit(2)
    tabm = tab(hdum[1].data)
    taba = tab(hdua[1].data)
    if tabm['TIME'][-1] <= taba['TIME'][0]:
        # if all events the list to be appended are later than the last event of primary
        newtab = vstack([tabm, taba])
    elif taba['TIME'][-1] <= tabm['TIME'][0]:
        # if the first event of the primary list is later than all events in the list to be appended
        newtab = vstack([taba, tabm])
    elif force:
        # something is wrong, sort the output
        vprint('The time interval of evta is already part of evtm')
        vprint('force == True, continuing now.')
        newtab = vstack([tabm, taba])
        newtab.sort('TIME')
    else:
        # something is wrong, sort the output 
        sys.stderr.write('error: The time interval of evta is already part of evtm, stopping now')
        sys.exit(2)
    # update the event table
    hdum[1].data = fits.table_to_hdu(newtab).data
    # update the GTI table
    tabm = tab(hdum[3].data)
    taba = tab(hdua[3].data)
    newtab = vstack([tabm, taba])
    newtab.sort('START')
    hdum[3].data = fits.table_to_hdu(newtab).data    
    # update the STDGTI table
    tabm = tab(hdum[4].data)
    taba = tab(hdua[4].data)
    newtab = vstack([tabm, taba])
    newtab.sort('START')
    hdum[4].data = fits.table_to_hdu(newtab).data
    print('Appended ' + fname_a + ' to  ', fname_m)
    hdum[0].header['HISTORY'] = 'Appended ' + fname_a
    hdum.writeto(out, overwrite=overwrite)



