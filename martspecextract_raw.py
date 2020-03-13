#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package. 
This is a script for internal useage only.
Essentially, it searches for the location of the brightest object
in the FOV, then extract the spectrum at that position using a 5 pixel radius.
No off-axis ARF, PSF, and vignetting corrections are applied here. 
"""

import sys
import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import interp1d
from scipy.interpolate import griddata as gd
from astropy import wcs
import astropy.coordinates as cd
import astropy.units as u
from martxclib.martxcfun import *
from matplotlib.colors import LogNorm
from astropy.table import Table as tab
from matplotlib import pyplot as plt
#caldb = check_caldb()


# ART-XC has 48 by 48 pixels
# with a pixel angular size of 0.724 arcmin
npixx = 48
npixy = 48
pixscale = 0.724

parser = Parser(description=__doc__,
	epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required arguments
parser.add_argument('-evt', type=str, required=True, 
	help='Name of the event file.')

parser.add_argument('-out', type=str, required=True,
	help='Name of the output file.')

parser.add_argument('-rmf', type=str, required=True, 
	help='Name of the RMF file.')
#default='art-xc_v0.0.arf',

# optional arguments
parser.add_argument('-arf', type=str, required=False, 
	help='Name of the ARF file.')
#default='art-xc_v0.0.arf',

parser.add_argument('-switch', type=bool, required=False, default=False,
	help='Swith on for applying off-axis ARF correction')

parser.add_argument('-vig', type=str, required=False, 
	help='Name of the Vignetting file.')
#default='artxc_vignetting.fits',
parser.add_argument('-psf', type=str, required=False, 
	help='Name of the PSF file.')

parser.add_argument(
	'-box', nargs='+', type=float, default=[],
	help="""list of raw_x, raw_y specifying the 4 corners of the box.
	Example: -box x1 x2 y1 y2""")


parser.add_argument(
	'-bbox', nargs='+', type=float, default=[],
	help="""list of raw_x, raw_y specifying the 4 corners of the background box.
	Example: -box x1 x2 y1 y2""")

# Global optinal arguments
parser.add_argument('-verbose', type=bool, required=False, default=False,
	help='Seeting verbosity.')

parser.add_argument('-overwrite', type=bool, required=False, default=False,
	help='Overwrite if set as True.')

parser.add_argument('-savepdf', type=bool, required=False, default=False,
	help='savepdf if set as True.')

args = parser.parse_args()

verbose = args.verbose
vprint = verboseprint(verbose)
overwrite = args.overwrite
savepdf = args.savepdf

evt = args.evt
switch = args.switch
arf = args.arf
out = args.out
psf = args.psf
vig = args.vig
rmf = args.rmf
box = args.box
bbox = args.bbox
bkgout = out[:-4] + 'bkg.fits' 


evthdu = fits.open(evt)
evttab = tab(evthdu[1].data)
# kick out events with weird PI (should've been done with martfilter 
# already, but just in case)
evttab = evttab[(evttab['PI']>=0) & (evttab['PI'] <= 511)]

expvalue = evthdu[0].header['EXPOSURE']

# flattend index of the detector pixel coordinate of each event
detcoor_id = evttab['RAW_X'].quantity.value * npixx + evttab['RAW_Y'].quantity.value

# detector mask
dmask = create_circular_mask(48,48,radius=24)


img, xx, yy = np.histogram2d(evttab['RAW_X'], evttab['RAW_Y'], bins=[np.arange(npixx + 1), np.arange(npixy +1)])
if np.min(img) == 0:
    ctmin = 1
else:
    ctmin = np.min(img)
ctmax = np.max(img)
if ctmax / ctmin <= 100:
	'''
	if the highest count pixel has no more than 100 times more counts than the lowest cont pixel,
	assuming everything is a background and extrat region is set at the center of the FOV
	'''
	positionx = 23
	positiony = 23
else:
	'''
	'''
	positiony = np.where(img == ctmax)[0][0]
	positionx = np.where(img == ctmax)[1][0]
	# avoid extracting region being outside the detector
	if (positiony - 5 <= 0):
		positiony += 5 - positiony
	elif (47 - positiony <= 5):
		positiony -= 47 - positiony
	if (positionx - 5 <= 0):
	    positionx += 5 - positionx
	elif (47 - positionx <= 5):
	    positionx -= 47 - positionx



if len(box) == 4:
	imask = create_box(48,48,box) * dmask
else:
	imask = create_circular_mask(48,48,radius=5,center=[positionx,positiony]) * dmask

backscal_src = float(np.sum(imask)) / np.sum(dmask)

if len(bbox) == 4:
	bmask = create_box(48,48,bbox) * dmask
else:	
	bmask_inner = np.invert(create_circular_mask(48,48,radius=10,center=[positionx,positiony]))
	bmask_outer = create_circular_mask(48,48,radius=19)
	# background would be extrated away from the detector boundary and the source
	bmask = bmask_inner * bmask_outer
	bkg_len = len(np.where(bmask_inner * bmask_outer)[0])
	bkg_id = int(np.random.uniform(0, len(np.where(bmask_inner * bmask_outer)[0])))
	bkg_positionx = np.where(bmask_inner * bmask_outer)[0][bkg_id]
	bkg_positiony = np.where(bmask_inner * bmask_outer)[1][bkg_id]
	bmask = create_circular_mask(48,48,radius=5,center=[bkg_positionx,bkg_positiony])

backscal_bkg = float(np.sum(bmask))/ np.sum(dmask)

areascal_bkg = float(np.sum(imask)/np.sum(bmask))
fig = plt.figure()
ax = fig.add_subplot(111)
cs = ax.imshow(img,  norm=LogNorm())
ax.set_xlim(0,48)
ax.set_ylim(0,48)
plt.minorticks_on()
cbar = fig.colorbar(cs)
plt.savefig(out + '.pdf',format='pdf',bbox_inches='tight')


mimg = img.copy()
mimg[~(imask)] = 0
fig = plt.figure()
ax = fig.add_subplot(111)
cs = ax.imshow(mimg,  norm=LogNorm())
ax.set_xlim(0,48)
ax.set_ylim(0,48)
cbar = fig.colorbar(cs)
plt.minorticks_on()
plt.savefig(out + '.src.pdf',format='pdf',bbox_inches='tight')



mimg = img.copy()
mimg[~(bmask)] = 0
fig = plt.figure()
ax = fig.add_subplot(111)
cs = ax.imshow(mimg,  norm=LogNorm())
ax.set_xlim(0,48)
ax.set_ylim(0,48)
plt.minorticks_on()
cbar = fig.colorbar(cs)
plt.savefig(out + '.bkg.pdf',format='pdf',bbox_inches='tight')


# off-axis angle in arcmin, assuming the optical axis center is at the center of the detector.
offaxis = np.sqrt((positionx - (npixx-1) * 0.5) ** 2 + (positiony - (npixx-1) * 0.5) ** 2 ) * pixscale
regradius = 5. * pixscale

if verbose:
	vprint('arguments passed:')
	vprint(args)

# If switch is turned on, an off-axis ARF with vignetting and PSF aperture corretion would be created.
if switch:
	# Read on-axis ARF and PSF, Vignetting files
	arfhdu = fits.open(arf)
	# define the name of the new arf file
	arfout = out[:-4] + 'arf' 
	arftab = arfhdu[1].data.copy()
	arfen = 0.5 * (arftab['ENERG_LO'] + arftab['ENERG_HI'])

	psfhdu = fits.open(psf)
	psftab = psfhdu[1].data
	psftab = psfhdu[1].data
	# energy group used when measuring EEF. See Swartz et al. (2013)
	psf_egrp = np.array([4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20,25, 30],dtype=float)
	psfen = np.unique(psftab['ENERGY'])
	psfoffaxis = np.unique(psftab['THETA'])
	psfrad = psftab[0]['RADIUS']
	psfoffaxis = np.unique(psftab['THETA']) / 60.
	# get index of off-axis array for a given off-axis in arcmin
	oang_func = interp1d(psfoffaxis, np.arange(len(psfoffaxis)), fill_value = 'extrapolate')
	# Simpy scaled by the distances from the supplied offaxis angle
	# to the two nearest offaxis angle in the psfoffaxis array.
	theta0 = psfoffaxis[np.floor(oang_func(offaxis)).astype(int)]
	theta1 = psfoffaxis[np.ceil(oang_func(offaxis)).astype(int)]
	if theta0 == theta1:
	    coeff0 = 1.
	    coeff1 = 0.
	else:
	    coeff0 = (offaxis - theta0)/(theta1 - theta0)
	    coeff1 = (theta1 - offaxis)/(theta1 - theta0)
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
	vigfunc = interp1d(vig_emed, np.arange(len(vig_emed)),fill_value='extrapolate')


	# Aperture correction
	# starting at unity, compute if radius is not None
	arfcorr_aper = np.zeros(len(arftab),dtype=float) + 1. 
	if regradius is not None:
		radius_arcmin = regradius
		for i in range(len(psfen)):
		    '''
		    For a given off-axis angle, find the EEF
		    '''
		    e0 = psf_egrp[i]
		    e1 = psf_egrp[i+1]
		    ee = psfen[i]
		    id_en = np.where((arfen >= e0) & (arfen < e1))[0]
		    eeftab0 = psftab[(psftab['ENERGY'] == ee) & (psftab['THETA'] == theta0 * 60.)]
		    eeftab1 = psftab[(psftab['ENERGY'] == ee) & (psftab['THETA'] == theta1 * 60.)]
		    eeffunc0 = interp1d(psfrad, eeftab0[0]['EFF'],fill_value='extrapolate')
		    eeffunc1 = interp1d(psfrad, eeftab1[0]['EFF'],fill_value='extrapolate')
		    arfcorr_aper[id_en] = coeff0 * eeffunc0(radius_arcmin) + coeff1 * eeffunc1(radius_arcmin)



	# Vignetting correction 
	# starting at unity, compute if offaxis is not 0.
	arfcorr_vig = np.zeros(len(arftab),dtype=float) + 1. 
	if offaxis != 0.:
		for i in progressbar(range(len(arfen))):
		    energy = arfen[i]
		    if (energy >= vig_elo[0]) and (energy <= vig_ehi[-1]):
		        grid_e, grid_th = np.meshgrid(np.array([energy]), vig_theta)
		        points = np.transpose(
		                              [np.tile(vig_emed, len(vig_theta)), \
		                               np.repeat(vig_theta, len(vig_emed))]
		                               )
		        values = vig_vig.flatten()
		        int_vig = gd(points, values, (grid_e, grid_th), method='nearest').flatten()
		        # define a vignetting function which takes an off-axis angle in arcmin
		        # and returns a vignetting fraction
		        vigfunc = interp1d(vig_theta*60., int_vig)
		        arfcorr_vig[i] = coeff0 * vigfunc(theta0) + coeff1 * vigfunc(theta1)


	arfhdu[1].data['SPECRESP'] = arftab['SPECRESP'] * arfcorr_vig * arfcorr_aper
	arfhdu[0].header.add_history('Applied aperture and off-axis corrections using martmkarf.py')
	arfhdu.writeto(arfout,overwrite=overwrite)
else:
	print('Not calculating off-axis arf')
	arfout = arf
	arfhdu = fits.open(arf)
	arftab = arfhdu[1].data.copy()

# Now extrat the source and background spectra



dspec = evttab[np.where(np.in1d(detcoor_id, np.where(imask.flatten())[0]))[0]]
dbkg = evttab[np.where(np.in1d(detcoor_id, np.where(bmask.flatten())[0]))[0]]
# spec_grp = dspec.groupby('PI')
# bkg_grp = dbkg.groupby('PI')
# spec = spec_grp['RAW_X'].count().values
# bkg = bkg_grp['RAW_X'].count().values

spec_grp = dspec.group_by('PI')
bkg_grp = dbkg.group_by('PI')
spec = spec_grp.groups.aggregate(len)['RAW_X'].quantity.value
bkg = bkg_grp.groups.aggregate(len)['RAW_X'].quantity.value
pi_spec = spec_grp.groups.aggregate(len)['PI'].quantity.value
pi_bkg = bkg_grp.groups.aggregate(len)['PI'].quantity.value



dfct = tab()
dfct['CHANNEL'] = np.arange(len(arftab))
dfct['COUNTS'] = np.zeros(len(arftab))

dfct['COUNTS'][pi_spec] = spec
print('here')

dfbk = tab()
dfbk['CHANNEL'] = np.arange(len(arftab))
dfbk['COUNTS'] = np.zeros(len(arftab))
dfbk['COUNTS'][pi_bkg] = bkg




ctlist = [
    fits.PrimaryHDU(),
    fits.table_to_hdu(dfct)
]

bllist = [
    fits.PrimaryHDU(),
    fits.table_to_hdu(dfbk)
]



ctlist[0].header['TELESCOP'] = 'SRG'
ctlist[0].header['INSTRUME'] = 'ART-XC'
ctlist[0].header['DETNAM'] = 'M1'
ctlist[0].header['FILTER'] = 'NONE'

ctlist[1].header['TELESCOP'] = 'SRG'
ctlist[1].header['INSTRUME'] = 'ART-XC'
ctlist[1].header['DETNAM'] = 'M1'
ctlist[1].header['EXPOSURE'] = expvalue
ctlist[1].header['FILTER'] = 'NONE'
ctlist[1].header['AREASCAL'] = 1
ctlist[1].header['BACKSCAL'] = backscal_src
ctlist[1].header['BACKFILE'] = bkgout
ctlist[1].header['RESPFILE'] = rmf
ctlist[1].header['ANCRFILE'] = arfout
ctlist[1].header['CHANTYPE'] = 'PI'
ctlist[1].header['EXTNAME'] = 'SPECTRUM'
ctlist[1].header['HDUCLASS'] = 'OGIP'
ctlist[1].header['HDUCLAS1'] = 'SPECTRUM'
ctlist[1].header['HDUVERS'] = '1.2.1'
ctlist[1].header['DETCHANS'] = len(arftab)
ctlist[1].header['HDUCLAS4'] = 'TYPE:I'
ctlist[1].header['CORRFILE'] = 'none'
ctlist[1].header['CORRSCAL'] = 1
ctlist[1].header['POISSERR '] = True

bllist[0].header['TELESCOP'] = 'SRG'
bllist[0].header['INSTRUME'] = 'ART-XC'
bllist[0].header['DETNAM'] = 'M1'
bllist[0].header['FILTER'] = 'NONE'

bllist[1].header['TELESCOP'] = 'SRG'
bllist[1].header['INSTRUME'] = 'ART-XC'
bllist[1].header['DETNAM'] = 'M1'
bllist[1].header['EXPOSURE'] = expvalue
bllist[1].header['FILTER'] = 'NONE'
bllist[1].header['AREASCAL'] = areascal_bkg
bllist[1].header['BACKSCAL'] = backscal_bkg
bllist[1].header['BACKFILE'] = bkgout
bllist[1].header['RESPFILE'] = rmf
bllist[1].header['ANCRFILE'] = arfout
bllist[1].header['CHANTYPE'] = 'PI'
bllist[1].header['EXTNAME'] = 'SPECTRUM'
bllist[1].header['HDUCLASS'] = 'OGIP'
bllist[1].header['HDUCLAS1'] = 'SPECTRUM'
bllist[1].header['HDUVERS'] = '1.2.1'
bllist[1].header['DETCHANS'] = len(arftab)
bllist[1].header['HDUCLAS4'] = 'TYPE:I'
bllist[1].header['CORRFILE'] = 'none'
bllist[1].header['CORRSCAL'] = 1
bllist[1].header['POISSERR '] = True


fits.HDUList(ctlist).writeto(out ,overwrite=overwrite)
fits.HDUList(bllist).writeto(bkgout,overwrite=overwrite)    
print('saved ' + out + ', and ', bkgout)






