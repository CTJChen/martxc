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
import astropy.coordinates as cd
import astropy.units as u
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
	'-evt', type=str, required=True,
	help='Name of the event file.')

parser.add_argument(
	'-out', type=str, required=True, default='expmap.fits',
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

parser.add_argument('-overwrite', type=bool, required=False, default=True,
	help='Overwrite if set as True.')

parser.add_argument(
	'-gti', type=str, required=False,
	help='Good Time Interval table file name.')

args = parser.parse_args()

verbose = args.verbose
vprint = verboseprint(verbose)

overwrite = args.overwrite

evt = args.evt
out = args.out
rasize = args.rasize
decsize = args.decsize
box = args.box

evttab = fits.getdata(evt)

# Get RA/DEC range
if (len(box) == 4) | (len(box) == 0):
	if len(box) == 4:
		ra1, ra2, dec1, dec2 = box
	else:
		ra1 = np.min(evttab['RA'])
		ra2 = np.max(evttab['RA'])
		dec1 = np.min(evttab['DEC'])
		dec2 = np.max(evttab['DEC'])
else:
	parser.error(
		'Must provide four radec values (e.g., \
		-box ra1 ra2 dec1 dec2)')

npixra = np.ceil((ra2 - ra1) / (rasize / 3600.)).astype(int)
npixdec = np.ceil((dec2 - dec1) / (decsize / 3600.)).astype(int)

rabins = np.linspace(0.5, float(npixra)+0.5, npixra+1, endpoint=True)
decbins = np.linspace(0.5, float(npixdec)+0.5, npixdec+1, endpoint=True)


# define WCS
w = wcs.WCS(naxis=2)
w.wcs.crpix = [npixra / 2., npixdec / 2.]
w.wcs.cdelt = np.array([rasize / 3600., decsize / 3600.])
w.wcs.crval = [(ra1 + ra2) / 2., (dec1 + dec2) / 2.]
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
w.wcs.cunit = ["deg", "deg"]
header_out = w.to_header()


vprint(
	'ra = [' + str(np.round(ra1, 4)) + ', ' +
	str(np.round(ra2, 4)) + ']')
vprint(
	'dec = [' + str(np.round(dec1, 4)) + ', ' +
	str(np.round(dec2, 4)) + ']')
vprint(
	'radecsize = [' + str(np.round(rasize, 2)) + ', ' +
	str(np.round(decsize, 2)) + '] (arcsec)')

c = cd.SkyCoord(evttab['RA'],evttab['DEC'],unit=(u.degree,u.degree))
rapix, decpix = w.world_to_pixel(c)
#img = np.zeros((npixdec,npixra))
rapix = np.floor(rapix).astype(int)
decpix = np.floor(decpix).astype(int)


# should include energy in the mask in the future.
mask = np.where((rapix <= npixra - 1) & (decpix <= npixdec -1))[0]
rapix = rapix[mask]
decpix = decpix[mask]

img, xx, yy = np.histogram2d(rapix, decpix, bins=[rabins, decbins])



hdu = fits.PrimaryHDU(img.transpose(), header=header_out)
hdu.writeto(out, overwrite=True)
