#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package. This script computes the exposuremap
for an arbitrary telescope with a circular field-of-view. \n
"""

import sys
import numpy as np
import astropy.io.fits as fits
#from tqdm import tqdm
from scipy import spatial
from scipy.interpolate import griddata as gd
from scipy.interpolate import interp1d
from astropy import wcs
import astropy.coordinates as cd
import astropy.units as u
from martxclib.martxcfun import *
from martxclib.martxcfits import rebinatt



parser = Parser(description=__doc__,
	epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)


# disabled for now. - CTC 20190620
# parser.add_argument('-xmlname', type=str, required=False,
# help='name of the telescope instrument specifications, an XML file.\n If not specified, the fov parameter should be specified.').
#

# Required arguments
parser.add_argument('-att', type=str, required=True,
	help='Name of the attitude file.')

parser.add_argument('-out', type=str, required=True, default='expmap.fits',
	help='Name of the output file.')

parser.add_argument('-fov', type=float, required=True,
	help='Field-of-view of the telescope to be simulated in arcmin. Supercedes the XML settings.')


# optional arguments

parser.add_argument('-vig', type=str, required=False,
	help='Name of the vignetting file. Only a fits file in compliance with the OGIP standard is supported (CAL/GEN/92-021). \
	If not provided, the script would only calculate the "raw" exposure map.')

parser.add_argument('-img', type=str, required=False,
	help='Name of an image file, the pixel size, RA/DEC range, and WCS keywords will be used. \
	If not provided, must specify all the following RA/DEC pixel sizes.')

parser.add_argument('-rasize', type=float, required=False, default=30.,
	help='RA pixel size in arcsec.')

parser.add_argument('-decsize', type=float, required=False, default=30.,
	help='DEC pixel size in arcsec.')

parser.add_argument(
	'-box', nargs='+', type=float, default=[],
	help="""list of ra dec values specifying the 4 corners of the box.
	Example: -box ra1 ra2 dec1 dec2""")

parser.add_argument('-energy', type=float, required=False,
	help='The user can specify the energy from which the vignetting function the script should use. Must specify -vig to do this. \
	If none is provided, the lowest energy vignetting function is used')

parser.add_argument('-attres', type=float, required=False,
	help='Time resolution in seconds. The input attitude file can be rebinned or interpolated at the supplied time resolution.\
	Re-binning/interpolating can speed up or slow down the computation time significantly.')

parser.add_argument('-time', type=float, required=False,
	help='Add a time constraint to the attitude file. Only the periods within the specified time would be considered.')

parser.add_argument('-offaxis', type=bool, required=False, default=False,
	help='If set as True, an additional fits image storing the exposure-time averaged off-axis \
	angle values would be added to the exposure map.')

parser.add_argument('-verbose', type=bool, required=False, default=True,
    help='set for verbosity')

parser.add_argument('-overwrite', type=bool, required=False, default=True,
	help='Overwrite if set as True.')

args = parser.parse_args()

verbose = args.verbose

vprint = verboseprint(verbose)
overwrite = args.overwrite
offaxis = args.offaxis
att = args.att
vig = args.vig
img = args.img
out = args.out
rasize = args.rasize
decsize = args.decsize
energy = args.energy
attres = args.attres
timelim = args.time
fov = args.fov
fovdeg = fov / 60.

if verbose:
	vprint('arguments passed:')
	vprint(args)


# Read attitude file 
atttab = fits.getdata(att)


# Get the vignetting function at the designated energy --
# only if a vignetting file name is provided
if vig is not None:
	vprint('Vignetting function will be used.\n')
	vigtab = fits.getdata(vig)		
	vig_theta = vigtab['THETA'][0]
	if len(np.shape(vigtab['VIGNET'][0])) == 2:
		vig_vig = vigtab['VIGNET'][0]
	elif len(np.shape(vigtab['VIGNET'][0])) == 3:
		vig_vig = vigtab['VIGNET'][0][0]
	vig_elo = vigtab['ENERG_LO'][0]
	vig_ehi = vigtab['ENERG_HI'][0]
	vig_emed = (vig_elo + vig_ehi) / 2.
	if energy is None:
		energy = vig_emed[0]
	# interpolate the vignetting function at the input energy
	grid_e, grid_th = np.meshgrid(np.array([energy]), vig_theta)
	points = np.transpose(
						  [np.tile(vig_emed, len(vig_theta)), \
						   np.repeat(vig_theta, len(vig_emed))]
						   )
	values = vig_vig.flatten()
	int_vig = gd(points, values, (grid_e, grid_th), method='nearest').flatten()
	# define a vignetting function which takes a off-axis angle in arcmin
	# and returns a vignetting fraction
	vigfunc = interp1d(vig_theta*60., int_vig)
else:
	vprint('Vignetting function NOT found, calculate a raw exposure map.')



# Get RA/DEC range and pixel sizes
if img is not None:
	# This would override all radec related arguments
	ihdu = fits.open(img)
	if len(ihdu) == 1:
		img = ihdu[0].data
		hdr = ihdu[0].header
	else:
		img = ihdu[1].data
		hdr = ihdu[1].header
	npixra = hdr['NAXIS1']
	npixdec = hdr['NAXIS2']
	rasize = hdr['CDELT1'] * 3600.
	decsize = hdr['CDELT2'] * 3600.
	w = wcs.WCS(hdr)
	ra1 = w.pixel_to_world(0,0).ra.deg
	dec1 = w.pixel_to_world(0,0).dec.deg
	ra2 = w.pixel_to_world(npixra-1,npixdec-1).ra.deg
	dec2 = w.pixel_to_world(npixra-1,npixdec-1).dec.deg
	header_out = w.to_header()
elif (len(args.box) == 4) | (len(args.box) == 0):
	if len(args.box) == 4:
		ra1, ra2, dec1, dec2 = args.box
	else:
		box = None
		ra1 = np.min(atttab['RA'])
		ra2= np.max(atttab['RA'])
		dec1= np.min(atttab['DEC'])
		dec2= np.max(atttab['DEC'])
	npixra = np.ceil((ra2 - ra1) / (rasize/3600.)).astype(int)
	npixdec = np.ceil((dec2 - dec1) / (decsize/3600.)).astype(int)
	# define WCS 
	w = wcs.WCS(naxis=2)
	w.wcs.crpix = [npixra/2.,npixdec/2.]
	w.wcs.cdelt = np.array([rasize/3600., decsize/3600.])
	w.wcs.crval = [(ra2 + ra1) * 0.5, (dec2 + dec1) * 0.5]
	w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
	header_out = w.to_header()
else:
	parser.error('Must provide four radec values (e.g., python test_args.py -box ra1 ra2 dec1 dec2)')


# only crop out attitudes if a box is used.
atttab = atttab[(atttab['RA'] >= ra1) & (atttab['RA'] <= ra2) & \
		(atttab['DEC'] >= dec1) & (atttab['DEC'] <= dec2)]

vprint(
	'ra = [' + str(np.round(ra1, 4)) + ', ' +
	str(np.round(ra2, 4)) + ']')
vprint(
	'dec = [' + str(np.round(dec1, 4)) + ', ' +
	str(np.round(dec2, 4)) + ']')
vprint(
	'radecsize = [' + str(np.round(rasize, 2)) + ', ' +
	str(np.round(decsize, 2)) + '] (arcsec)')



# work only with attitude entries within the area of interests
if timelim is not None:
	atttab = atttab[atttab['TIME'] <= atttab['TIME'][0] + timelim]

atttime = atttab['TIME']
attra = atttab['RA']
attdec = atttab['DEC']

if attres is None:
	'''
	Assuming the time steps in the attitude file is a constant
	'''
	attres = atttime[1] - atttime[0]
	vprint('no attres parameter given, no attitude rebinning would be carried out')
else:
	vprint('re-bininng the attitude file using the new time step')
newatttime, newattra, newattdec = rebinatt(atttime, attra, attdec, attres)
del(atttime, attra, attdec)


# Work on the meshgrid -- use cKDTree to find the mesh points that's 
# within the FOV radius at a given pointing position
x = np.arange(npixra)
y = np.arange(npixdec)
X, Y = np.meshgrid(x, y)
ras, decs = w.wcs_pix2world(X, Y, 0)
ras_1d = ras.flatten()
decs_1d = decs.flatten()
coordinates = np.c_[ras.ravel(), decs.ravel()]
cat_grid = cd.SkyCoord(ras_1d,decs_1d,unit=(u.degree,u.degree))
expmap = np.zeros(npixra * npixdec)

'''
Loop through all attitude entries
(rebinned according to the time resolution parameter, if exists)
The following tasks were done within this loop:
1. find the mesh points within the field-of-view radius
2. Compute the Vignetting values at each position -- if 
3. 
'''

if vig is None:
	vprint('Raw exposure map (not corrected for Vignetting effects)')
	for aid in progressbar(range(len(newattra))):
		pnt = cd.SkyCoord(newattra[aid],newattdec[aid],unit=(u.degree,u.degree))
		dist = pnt.separation(cat_grid)
		ix = np.where(dist.deg <= fovdeg)[0]
		expmap[ix] += attres
else:
	vprint('Vignetted exposure map')
	for aid in progressbar(range(len(newattra))):
		pnt = cd.SkyCoord(newattra[aid],newattdec[aid],unit=(u.degree,u.degree))
		dist = pnt.separation(cat_grid)
		ix = np.where(dist.deg <= fovdeg)[0]
		expmap[ix] += vigfunc(dist[ix].arcmin) * attres



if offaxis is True:
	'''
	This would calculate exposure time weighted off-axis angle at each pixel
	'''
	offaxistab = np.zeros(npixra * npixdec)

	vprint('Make an additional table storing off-axis values')
	for aid in progressbar(range(len(newattra))):
		pnt = cd.SkyCoord(newattra[aid],newattdec[aid],unit=(u.degree,u.degree))
		dist = pnt.separation(cat_grid)
		ix = np.where(dist.deg <= fovdeg)[0]
		# Add exposure-time normalized off-axis values, incrementally.
		offaxistab[ix] += (dist[ix].arcmin * vigfunc(dist[ix].arcmin) * attres / expmap[ix])
	hdulist = [fits.PrimaryHDU(expmap.reshape(npixdec,npixra), header=header_out)]
	hdulist.append(fits.ImageHDU(offaxistab.reshape(npixdec,npixra), header=header_out))
	explist = fits.HDUList(hdulist)
	explist[0].name = 'expmap'
	explist[1].name = 'offaxis'
	explist[1].header.add_his
	explist.writeto(out,overwrite=overwrite)
else:
	# header is an astropy.io.fits.Header object.  We can use it to create a new
	# PrimaryHDU and write it to a file.
	hdu = fits.PrimaryHDU(expmap.reshape(npixdec,npixra), header=header_out)
	# Save to FITS file
	hdu.writeto(out, overwrite=overwrite)
