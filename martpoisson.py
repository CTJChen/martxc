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
from martxclib.martxcfun import *
from astropy import wcs



parser = Parser(description=__doc__,
	epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)


# disabled for now. - CTC 20190620
# parser.add_argument('-xmlname', type=str, required=False,
# help='name of the telescope instrument specifications, an XML file.\n If not specified, the fov parameter should be specified.').
#

# Required arguments
parser.add_argument('-img', type=str, required=True,
	help='Name of the event file.')

parser.add_argument('-out', type=str, required=True, default='output.fits',
	help='Name of the output file.')

parser.add_argument('-combine', type=bool, required=True, default=True, 
	help='Set to 0 would return the poisson background only. Otherwise a image combined with background would be saved.')

parser.add_argument('-nrate', type=float, required=True, default=1e-3, 
	help='Noise count rate in count per second per pixel.')

# optional arguments
parser.add_argument('-exptime', type=float, required=False,
	help='Exposure value for calculating the background.')

parser.add_argument('-exp', type=str, required=False,
	help='Name of the exposure map.')

parser.add_argument('-verbose', type=bool, required=False, default=True,
	help='Add a time constraint to the attitude file. Only the periods within the specified time would be considered.')

parser.add_argument('-overwrite', type=bool, required=False, default=True,
	help='Overwrite if set as True.')

args = parser.parse_args()

verbose = args.verbose
vprint = verboseprint(verbose)
vprint(args)

overwrite = args.overwrite

img = args.img
out = args.out
combine = args.combine
exptime = args.exptime
exp = args.exp
nrate = args.nrate


# noise function
def pnoise(exp, nrate=nrate):
	return np.random.poisson(nrate * exp)


# Reade image, some fits files have a blanck primary HDU, check this
ihdu = fits.open(img)
if len(ihdu) == 1:
	imgtab = ihdu[0].data.copy()
else:
	imgtab = ihdu[1].data.copy()

if exp is not None:
	vprint('Use exposuremap')
	ehdu = fits.open(exp)
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
	vprint('Use exposure time')
	ishape = np.shape(imgtab)
	noise = np.random.poisson(nrate * exptime, size=ishape)


if combine:
	imgtab += noise
else:
	imgtab = noise

if len(ihdu) == 1:
	ihdu[0].data = imgtab
else:
	ihdu[1].data = imgtab


ihdu.writeto(out, overwrite=overwrite)
