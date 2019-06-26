#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package.
This script computes an image based on a given event file. \n
Currently only works with a simulated event list.
"""

import sys
import numpy as np
import astropy.io.fits as fits
import argparse
import matplotlib.pyplot as plt
from astropy import wcs
from martxclib.martxcfun import *


class HelpfulParser(argparse.ArgumentParser):
	'''
	Handling errors
	'''
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		sys.stderr.write('run \"python martev2img.py -h\" for helps\n')
		sys.exit(2)


parser = HelpfulParser(
	description=__doc__,
	epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)


# Required arguments
parser.add_argument(
	'-evtname', type=str, required=True,
	help='Name of the event file.')

parser.add_argument(
	'-outname', type=str, required=True, default='expmap.fits',
	help='Name of the output file.')

parser.add_argument(
	'-rasize', type=float, required=True, default=30.,
	help='RA pixel size in arcsec.')

parser.add_argument(
	'-decsize', type=float, required=True, default=30.,
	help='DEC pixel size in arcsec.')


# optional arguments

parser.add_argument(
	'-box', nargs='+', type=float, default=[],
	help="""list of ra dec values specifying the 4 corners of the box.
	Example: -box ra1 ra2 dec1 dec2""")

parser.add_argument('-verbose', type=bool, required=False, default=True,
	help='Add a time constraint to the attitude file. Only the periods within the specified time would be considered.')

#
#parser.add_argument(
#	'-energy', nargs='+', type=float, default=[],
#	help="""Energy range of events to be used.
#	Multiple energy ranges are allowed, but the inputs must include
#	both the lower and upper bounds of each energy range.
#	Example: '-energy 5 10 12 30' would use events with energy between 5-30 keV
#	except those between 10-12 keV.""")
#

parser.add_argument(
	'-gtiname', type=str, required=False,
	help='Good Time Interval table file name.')

args = parser.parse_args()

verbose = args.verbose
vprint = verboseprint(verbose)


evtname = args.evtname
outname = args.outname
rasize = args.rasize
decsize = args.decsize
box = args.box

evttab = fits.getdata(evtname)

# Get RA/DEC range
if (len(box) == 4) | (len(box) == 0):
	if len(box) == 4:
		ra1, ra2, dec1, dec2 = args.box
	else:
		ra1 = np.min(evttab['RA'])
		ra2 = np.max(evttab['RA'])
		dec1 = np.min(evttab['DEC'])
		dec2 = np.max(evttab['DEC'])
	npixra = np.ceil((ra2 - ra1) / (rasize / 3600.)).astype(int)
	npixdec = np.ceil((dec2 - dec1) / (decsize / 3600.)).astype(int)
	# define WCS
	w = wcs.WCS(naxis=2)
	w.wcs.crpix = [npixra / 2., npixdec / 2.]
	w.wcs.cdelt = np.array([rasize / 3600., decsize / 3600.])
	w.wcs.crval = [(ra1 + ra2) / 2., (dec1 + dec2) / 2.]
	w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
	header_out = w.to_header()
else:
	parser.error(
		'Must provide four radec values (e.g., \
		-box ra1 ra2 dec1 dec2)')

vprint(
	'ra = [' + str(np.round(ra1, 4)) + ', ' +
	str(np.round(ra2, 4)) + ']')
vprint(
	'dec = [' + str(np.round(dec1, 4)) + ', ' +
	str(np.round(dec2, 4)) + ']')
vprint(
	'radecsize = [' + str(np.round(rasize, 2)) + ', ' +
	str(np.round(decsize, 2)) + '] (arcsec)')

ragrid = np.linspace(ra1, ra2, npixra)
decgrid = np.linspace(dec1, dec2, npixdec)
ras, decs = np.meshgrid(ragrid, decgrid)
ras_1d = ras.flatten()
decs_1d = decs.flatten()

h = plt.hist2d(
	np.asrray(evttab['RA']),
	np.asarray(evttab['DEC']), bins=[npixra, npixdec])

hdu = fits.PrimaryHDU(h[0], header=header_out)
hdu.writeto(outname, overwrite=True)
