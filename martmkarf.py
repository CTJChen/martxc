#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package. This script applies
aperture and off-axis corrections based on the PSF and Vignetting functions
"""
import sys
import numpy as np
import astropy.io.fits as fits
import argparse
from scipy.interpolate import interp1d
from scipy.interpolate import griddata as gd
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
		sys.stderr.write('run \"python martexpmap.py -h\" for helps\n')
		sys.exit(2)



parser = HelpfulParser(description=__doc__,
	epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required arguments
parser.add_argument('-arf', type=str, required=True, default='art-xc_v0.0.arf',
	help='Name of the ARF file.')

parser.add_argument('-vig', type=str, required=True, default='artxc_vignetting.fits',
	help='Name of the Vignetting file.')

parser.add_argument('-psf', type=str, required=True, default='artxc_psf_eef.fits',
	help='Name of the PSF file.')

parser.add_argument('-out', type=str, required=True, default='example_output.arf',
	help='Name of the output ARF file.')

# optional arguments
parser.add_argument('-img', type=str, required=False, default=None,
	help='Image name, should include the on-axis pointing position. Works with region file')

parser.add_argument('-region', type=str, required=False,
	help='Name of the DS9 region file. Only supports one circular region with WCS values for now. \
	Takes precedence over all other spatial parameters.')

parser.add_argument('-offaxis', type=float, required=False, default=0., 
	help='Off-axis angle of the pointing position in arcmin. \
	Only works if -region and -imgname parameters are both None.')

parser.add_argument('-radius', type=float, required=False, default=None,
	help='Radius of a source extraction region in arcsec, does not work if a \
	DS9 region file is provided.')

parser.add_argument(
	'-box', nargs='+', type=float, default=[],
	help="""list of ra dec values specifying the 4 corners of a region box.
	Does not work with the -region parameter.""")

parser.add_argument('-verbose', type=bool, required=False, default=True,
	help='Seeting verbosity.')

parser.add_argument('-overwrite', type=bool, required=False, default=True,
	help='Overwrite if set as True.')


args = parser.parse_args()

verbose = args.verbose
vprint = verboseprint(verbose)
overwrite = args.overwrite

arf = args.arf
out = args.out
psf = args.psf
vig = args.vig
region = args.region
offaxis = args.offaxis
img = args.img
box = args.box

if verbose:
	vprint('arguments passed:')
	vprint(args)


# Read on-axis ARF and PSF, Vignetting files
arfhdu = fits.open(arf)
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





if (region is not None) and (img is not None):
	coordstr = parse_ds9(region)
	radec = cd.SkyCoord(coordstr[1][0] + ' ' + coordstr[1][1],frame = coordstr[0], unit=(u.hourangle, u.deg))  # in degree
	radius = float(re.findall(r'\d+\.\d+', coordstr[1][-1])[0]) # in radius
	'''
	Use the CRVAL keywords as the pointing position for now, 
	WILL be updated in the future. 
	'''
	imghdu = fits.open(img)
	if len(imghdu) == 1:
		ra_pnt = imghdu[0].header['CRVAL1']
		dec_pnt = imghdu[0].header['CRVAL2']
	else:
		if 'CRVAL1' in list(imghdu[0].header.keys()):
			ra_pnt = imghdu[0].header['CRVAL1']
			dec_pnt = imghdu[0].header['CRVAL2']
		else:
			ra_pnt = imghdu[1].header['CRVAL1']
			dec_pnt = imghdu[1].header['CRVAL2']
	radec_pnt = cd.SkyCoord(ra_pnt, dec_pnt, unit=(u.degree, u.degree))
	offaxis = radec_pnt.separation(radec).arcmin
elif (region is None) and (img is None):
	vprint('No DS9 region and image, continue. ')
else:
	sys.stderr.write('Must provide both -region and -img parameters, \n')
	sys.stderr.write('when making an ARF for a source-extraction region. ')
	sys.exit(2)




# Aperture correction
# starting at unity, compute if radius is not None
arfcorr_aper = np.zeros(len(arftab),dtype=float) + 1. 
if radius is not None:
	radius_arcmin = radius / 60.
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
arfhdu.writeto(out,overwrite=overwrite)








