#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC analysis package. 
This script takes an image file and compute 
based on a given event file. \n
Currently calculate on a box with a given RA/DEC range that can be specified.\n
Since Eulerian distance is used instead of angular distance, the boxsize should
be limited to several deg.
"""

import sys
import numpy as np
import astropy.io.fits as fits
import argparse
from martxclib.martxcfun import *
from astropy import wcs


class HelpfulParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser = HelpfulParser(description=__doc__,
	epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)


# disabled for now. - CTC 20190620
# parser.add_argument('-xmlname', type=str, required=False,
# help='name of the telescope instrument specifications, an XML file.\n If not specified, the fov parameter should be specified.').
#

# Required arguments
parser.add_argument('-imgname', type=str, required=True,
	help='Name of the event file.')

parser.add_argument('-outname', type=str, required=True, default='output.fits',
	help='Name of the output file.')

parser.add_argument('-combine', type=bool, required=True, default=True, 
	help='Set to 0 would return the poisson background only. Otherwise a image combined with background would be saved.')

parser.add_argument('-nrate', type=bool, required=True, default=1e-3, 
	help='Noise count rate in count per second per pixel.')

# optional arguments
parser.add_argument('-exptime', type=float, required=False,
	help='Exposure value for calculating the background.')

parser.add_argument('-expname', type=str, required=False,
	help='Name of the exposure map.')

args = parser.parse_args()

print(args)
imgname = args.imgname
outname = args.outname
combine = args.combine
exptime = args.exptime
expname = args.expname
nrate = args.nrate


# noise function
def pnoise(exp, nrate=nrate):
	return np.random.poisson(nrate * exp)


# Reade image, some fits files have a blanck primary HDU, check this
ihdu = fits.open(imgname)
if len(ihdu) == 1:
	img = ihdu[0].data
else:
	img = ihdu[1].data
	if expname is not None:
		ehdu = fits.open(expname)
		if len(ehdu) == 1:
			expmap = ehdu[0].data
		else:
			expmap = ehdu[1].data
		exptime = np.max(expmap)
		# Now draw poisson noise based on expvalue on a pixel-by-pixel basis
		noise = pnoise(expmap)
		# if exposure value is zero, assume no backgrounds.
		inoexp = np.where(expmap == 0)
		noise[inoexp] = 0.
	else:
		# if no exp fits file is provided,
		# use image fits file and assume const. exp time.
		ishape = np.shape(img)
		noise = np.random.poisson(nrate * exptime, size=ishape)


if combine:
	img += noise
else:
	img = noise

if len(ihdu) == 1:
	ihdu[0].data = img
else:
	ihdu[1].data = img


ihdu.writeto(outname, overwrite=True)
