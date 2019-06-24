#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package. This script computes the exposuremap
for an arbitrary telescope with a circular field-of-view. \n
Currently calculate on a box with a given RA/DEC range that can be specified.\n
Since Eulerian distance is used instead of angular distance, the boxsize should
be limited to several deg.
"""

import sys
import numpy as np
import astropy.io.fits as fits
import argparse
from tqdm import tqdm
from scipy import spatial
from scipy.interpolate import griddata as gd
from scipy.interpolate import interp1d
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
# parser.add_argument('--xmlname', type=str, required=False,
# help='name of the telescope instrument specifications, an XML file.\n If not specified, the fov parameter should be specified.').
#

# Required arguments
parser.add_argument('--attname', type=str, required=True,
	help='Name of the attitude file.')

parser.add_argument('--outname', type=str, required=True, default='expmap.fits',
	help='Name of the output file.')

parser.add_argument('--fov', type=float, required=True,
	help='Field-of-view of the telescope to be simulated in arcmin. Supercedes the XML settings.')

parser.add_argument('--rasize', type=float, required=True, default=30.,
	help='RA pixel size in arcsec.')

parser.add_argument('--decsize', type=float, required=True, default=30.,
	help='DEC pixel size in arcsec.')


# optional arguments

parser.add_argument('--vigname', type=str, required=False,
	help='Name of the vignetting file. Only a fits file in compliance with the OGIP standard is supported (CAL/GEN/92-021). \
	If not provided, the script would only calculate the "raw" exposure map.')

parser.add_argument('--ra1', type=float, required=False,
	help='RA lower limit (in deg) for the area of interests.')

parser.add_argument('--ra2', type=float, required=False,
	help='RA upper limit (in deg) for the area of interests.')

parser.add_argument('--dec1', type=float, required=False,
	help='DEC lower limit (in deg) for the area of interests.')

parser.add_argument('--dec2', type=float, required=False,
	help='DEC upper limit (in deg) for the area of interests.')

parser.add_argument('--energy', type=float, required=False,
	help='The user can specify the energy from which the vignetting function the script should use. Must specify --vigname to do this. \
	If none is provided, the lowest energy vignetting function is used')

parser.add_argument('--attres', type=float, required=False,
	help='Time resolution in seconds. The input attitude file can be rebinned or interpolated at the supplied time resolution.\
	Re-binning/interpolating can speed up or slow down the computation time significantly.')

args = parser.parse_args()

print(args)
attname = args.attname
vigname = args.vigname
outname = args.outname
rasize = args.rasize
decsize = args.decsize
energy = args.energy
attres = args.attres
fov = args.fov
fovdeg = fov / 60.

# Get the vignetting function at the designated energy --
# only if a vignetting file name is provided
if vigname is not None:
	print('Vignetting function will be used.')
	vigtab = fits.getdata(vigname)
	vig_theta = vigtab['THETA'][0]
	if len(np.shape(vigtab['VIGNET'][0])) == 2:
		vig_vig = vigtab['VIGNET'][0]
	elif len(np.shape(vigtab['VIGNET'][0])) == 3:
		vig_vig = vigtab['VIGNET'][0][0]
	vig_elo = vigtab['ENERG_LO'][0]
	vig_ehi = vigtab['ENERG_HI'][0]
	vig_emed = (vig_elo + vig_ehi) / 2.

# interpolate the vignetting function at the input energy
	grid_e, grid_th = np.meshgrid(np.array([energy]), vig_theta)
	points = np.transpose(
	                      [np.tile(vig_emed, len(vig_theta)), \
	                       np.repeat(vig_theta, len(vig_emed))]
	                       )
	values = vig_vig.flatten()
	int_vig = gd(points, values, (grid_e, grid_th), method='nearest').flatten()

	vigfunc = interp1d(vig_theta*60., int_vig)

	def vig2d(ra, dec, pnt):
		'''
		For a given pointing position,
		calculate the corresponding Vignetting values at the input RA/DEC values.
		ra = x-axis
		dec = y-axis
		pnt = pointing position, with array([ra_pnt, dec_pnt]), in degrees
		'''
		pntra, pntdec = pnt
		radist = ra - pntra
		decdist = dec - pntdec
		dist = np.sqrt((radist)**2 + (decdist)**2) * 60.  # in arcmin
		vigout = np.zeros_like(dist)
		vigout = vigfunc(dist)
		return vigout
else:
	print('Vignetting function NOT found, calculate a raw exposure map.')


# The function for rebinning the attitude file
def rebinatt(atttime, attra, attdec, attres):
	'''
	#dist -- not used for now.
	Rebin the ATT array according to the new attres.
	Based on this scipy coockbook:
	https://scipy-cookbook.readthedocs.io/items/Rebinning.html
	'''
	a = np.vstack((atttime, attra, attdec))
	oldres = atttime[1] - atttime[0]
	if oldres == attres:
		return a
	else:
		k = int(oldres / attres)
		newlen = len(atttime) * k
		newshape = (3, newlen)
		assert len(a.shape) == len(newshape)
		slices = [slice(0, old, float(old) / new) for old, new in zip(a.shape, newshape)]
		coordinates = np.mgrid[slices]
		indices = coordinates.astype('i')  # choose the biggest smaller integer index
		newatttime, newattra, newattdec = a[tuple(indices)]
		newatttime = np.linspace(newatttime[0], newatttime[-1], len(newatttime))
		for raid in range(len(attra) - 1):
			newattra[raid * k:(raid + 1) * k] = \
			np.linspace(attra[raid], attra[raid + 1], k)
		for decid in range(len(attdec) - 1):
			newattdec[decid * k:(decid+1) * k] = \
			np.linspace(attdec[decid], attdec[decid + 1], k)
		return newatttime, newattra, newattdec


atttab = fits.getdata(attname)

# use the RA/DEC range in the attitude files
if args.ra1 is None:
	ra1= np.min(atttab['RA'])
if args.ra2 is None:
	ra2= np.max(atttab['RA'])
if args.dec1 is None:
	dec1= np.min(atttab['DEC'])
if args.dec1 is None:
	dec2= np.max(atttab['DEC'])

npixra = np.ceil((ra2 - ra1) / (rasize/3600.)).astype(int)
npixdec = np.ceil((dec2 - dec1) / (decsize/3600.)).astype(int)
ragrid = np.linspace(ra1, ra2, npixra)
decgrid = np.linspace(dec1, dec2, npixdec)
ras, decs = np.meshgrid(ragrid, decgrid)
ras_1d = ras.flatten()
decs_1d = decs.flatten()

atttab = atttab[(atttab['RA'] >= ra1) & (atttab['RA'] <= ra2) & \
                (atttab['DEC'] >= dec1) & (atttab['DEC'] <= dec2)]
atttime = atttab['TIME']
attra = atttab['RA']
attdec = atttab['DEC']

if attres is None:
    '''
    Assuming the time steps in the attitude file is a constant
    '''
    attres = atttime[1] - atttime[0]
    print('no attres parameter given, no attitude rebinning would be carried out')
else:
	print('re-bininng the attitude file using the new time step')
newatttime, newattra, newattdec = rebinatt(atttime, attra, attdec, attres)
del(atttime, attra, attdec)

# Work on the meshgrid -- use cKDTree to find the mesh points that's 
# within the FOV radius at a given pointing position
coordinates = np.c_[ras.ravel(), decs.ravel()]
tree = spatial.cKDTree(coordinates)

expmap = np.zeros(np.shape(coordinates)[0])

'''
Loop through all attitude entries
(rebinned according to the time resolution parameter, if exists)
The following tasks were done within this loop:
1. find the mesh points within the field-of-view radius
2. Compute the Vignetting values at each position -- if 
3. 
'''

if vigname is None:
    print('Raw exposure map (not corrected for Vignetting effects)')
    for aid in tqdm(range(len(newattra))):
        pnt = np.array([newattra[aid],newattdec[aid]])
        ix = tree.query_ball_point(pnt, fovdeg)
        expmap[ix] += attres
else:
    print('Vignetted exposure map')
    for aid in tqdm(range(len(newattra))):
        pnt = np.array([newattra[aid],newattdec[aid]])
        ix = tree.query_ball_point(pnt, fovdeg)
        expmap[ix] += vig2d(ras_1d[ix], decs_1d[ix], pnt) * attres

# define WCS 
w = wcs.WCS(naxis=2)
w.wcs.crpix = [npixra/2.,npixdec/2.]
w.wcs.cdelt = np.array([rasize/3600., decsize/3600.])
w.wcs.crval = [np.median(newattra), np.median(newattdec)]
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
header = w.to_header()

# header is an astropy.io.fits.Header object.  We can use it to create a new
# PrimaryHDU and write it to a file.
hdu = fits.PrimaryHDU(expmap.reshape(npixra,npixdec), header=header)
# Save to FITS file
hdu.writeto(outname, overwrite=True)
