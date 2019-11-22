#!/usr/bin/env python
# -*- coding: utf-8 -*-
__doc__ = """
Part of the MSFC ART-XC software package. 
This is a script for internal useage only.
It searches for the location of the brightest object in the FOV,
then find the range of azimuthal angles for the region
to be used for calculating the encircled energy fraction. 
Output - a set of encircled total counts for a "partial circle" 
centered at the peak pixel, as a function of 
radius = np.arange(20) + 1
Also saves a pdf for the partial circle image at radius = 20
"""
import numpy as np
import astropy.io.fits as fits
from astropy import wcs
from martxclib.martxcfun import *
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.weight'] = 'bold'#'bold','medium'
plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
plt.rcParams['legend.fontsize'] = 0.6*plt.rcParams['font.size']
plt.rcParams['xtick.labelsize'] = 0.9*plt.rcParams['font.size']
plt.rcParams['ytick.labelsize'] = 0.9*plt.rcParams['font.size']
#plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.width'] = 1
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.loc'] = 'center left'
plt.rcParams['axes.linewidth'] = 1.0
plt.gca().spines['right'].set_color('none')
plt.gca().spines['top'].set_color('none')
plt.gca().xaxis.set_ticks_position('bottom')
plt.gca().yaxis.set_ticks_position('left')
from matplotlib.colors import LogNorm
from astropy.table import Table as tab
from scipy.interpolate import interp1d

#caldb = check_caldb()

# ART-XC has 48 by 48 pixels
# with a pixel angular size of 0.724 arcmin
npixx = 48
npixy = 48

centx = int(npixx / 2)
centy = int(npixy / 2)

pixscale = 0.724

parser = Parser(description=__doc__,
    epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required arguments
parser.add_argument('-evt', type=str, required=True, 
    help='Name of the event file.')

parser.add_argument('-out', type=str, required=True,
    help='Name of the output file.')

# optional arguments
parser.add_argument(
    '-xy', nargs='+', type=float, required=False, default=[],
    help="""raw_x and raw_y coordinates""")

# Global optional arguments
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
out = args.out
xy = args.xy

evthdu = fits.open(evt, ignore_missing_end=True)
evttab = tab(evthdu[1].data)
evttab = evttab[(evttab['PI']>=0) & (evttab['PI'] <= 511)]

# flattend index of the detector pixels each event
detcoor_id = evttab['RAW_X'].quantity.value * npixx + evttab['RAW_Y'].quantity.value
dmask = create_circular_mask(48,48,radius=24)


img, xx, yy = np.histogram2d(evttab['RAW_X'], evttab['RAW_Y'], bins=[np.arange(npixx + 1), np.arange(npixy +1)])

if np.min(img) == 0:
    ctmin = 1
else:
    ctmin = np.min(img)
ctmax = np.max(img)

if len(xy) > 0:
    '''
    use raw_x raw_y from arguments if possible
    '''
    positionx = xy[0]
    positiony = xy[1]
elif ctmax / ctmin <= 100:
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

# 
bmask_inner = np.invert(create_circular_mask(48,48,radius=15,center=[positionx,positiony]))
# background value per pixel to be used
bmask = bmask_inner * dmask
bkg_median = np.median(img[bmask])

d=np.sqrt((positionx - centx)**2 + (positiony - centy)**2)
'''
if d + 20 >= 24:
    invert = False
    angle_min, angle_max = determin_angles(centx, centy, 24, positionx, positiony, 20)
    if angle_min > angle_max: angle_min, angle_max = angle_max, angle_min
    if (angle_min <= 90) & (angle_max <= 180): invert = True
    ctpix = [np.max(img)]
    ctsum = [0]
    for radius in np.arange(20)+1:
        pmask = create_partial_circle_mask(48,48,radius=radius,center=[positionx,positiony], 
            angle_min=angle_min, angle_max=angle_max, invert=invert)
        ctsum.append(np.sum(img[pmask]) - bkg_median * len(img[pmask]))
        ctpix.append((np.sum(img[pmask]) - bkg_median * len(img[pmask])) / len(img[pmask]))
else:
    ctsum = [0]
    ctpix = [np.max(img)]

    for radius in np.arange(20)+1:
        pmask = create_circular_mask(48,48,radius=radius,center=[positionx,positiony])
        ctsum.append(np.sum(img[pmask]) - bkg_median * len(img[pmask]))    
        ctpix.append((np.sum(img[pmask]) - bkg_median * len(img[pmask])) / len(img[pmask]))
'''
ctsum = [0]
ctpix = [np.max(img)]

for radius in np.arange(20)+1:
    pmask = create_circular_mask(48,48,radius=radius,center=[positionx,positiony])
    ctsum.append(np.sum(img[pmask]) - bkg_median * len(img[pmask]))    
    ctpix.append((np.sum(img[pmask]) - bkg_median * len(img[pmask])) / len(img[pmask]))

xdata = (np.arange(21))
ydata = np.array(ctsum) / max(ctsum)
#eef = pchip(xdata,ydata)
eef = interp1d(xdata, ydata,fill_value='extrapolate')


eefarr = eef(np.linspace(0,20,1000))
xarr =  np.linspace(0,20,1000)
hpd = np.round(xarr[np.argmin(np.abs(eefarr - 0.5))] * pixscale, 4)
print(evt + ' HPD = ' + str(hpd) + ' arcmin')

if savepdf:
    '''
    plt.cla()    
    plt.plot(xdata * pixscale, ydata, 'o',linestyle='-')
    plt.plot(xarr * pixscale, eef(xarr),color='red')
    plt.xlim(0.,15.)
    plt.ylim(0.45,1.0)
    plt.title(evt)
    plt.minorticks_on()
    plt.text(hpd + 0.5, 0.5 , 'HPD = ' + str(hpd))
    plt.xlabel('Radius (arcmin)')
    plt.ylabel('EEF')
    '''
    plt.cla()
    fig, axarr = plt.subplots(2, 1, sharex=True, figsize=(8, 8))
    fig.subplots_adjust(hspace=0.02)
    #ax.title('Including color prior reduce spurious rate')

    ax1 = axarr[0]
    ax1.plot(np.arange(21), ctpix, 'o',linestyle='-')
    ax1.set_yscale('log')
    ax1.set_xlim(0.,21)
    ax1.set_ylabel('Counts / Pixel')
    ax1.text(12.5, np.max(ctpix)/2, 'Center = (' + str(positionx) + ', ' + str(positiony) + ')')
    ax1.text(12.5, np.max(ctpix)/10, 'Median background \n(cts/pix) = ' + str(bkg_median))

    ax = axarr[1]
    ax.plot(np.arange(21), ydata, 'o',linestyle='-')
    ax.plot(np.linspace(0,25,1000), eef(np.linspace(0,25,1000)),color='green')
    ax.set_xlabel('Radious (detector)')
    ax.set_ylabel('EEF')
    ax.text(12.5, 0.5 , 'HPD = ' + str(np.round(hpd,2)) + ' arcmin')
    plt.minorticks_on()    
    plt.savefig('eef' + out + '.pdf',format='pdf',bbox_inches='tight')

if savepdf:
    fig, axarr = plt.subplots(1, 2, figsize=(9, 4))
    ax1 = axarr[0]
    ax2 = axarr[1]

    mimg = img.copy()
    mimg[~(pmask)] = 0
    im = ax1.imshow(img,  norm=LogNorm())
    im = ax2.imshow(mimg,  norm=LogNorm())
    ax1.set_xlim(0,48)
    ax1.set_ylim(0,48)
    ax2.set_xlim(0,48)
    ax2.set_ylim(0,48)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.01, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.title(evt)
    plt.savefig('img' + out + '.pdf',format='pdf',bbox_inches='tight')

np.savetxt(out + 'eef.txt', ydata)
np.savetxt(out + 'cts.txt', ctsum)


