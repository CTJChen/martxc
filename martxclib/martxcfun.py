"""
Functions for MSFC ART-XC
"""
from __future__ import print_function, division
import numpy as np
import sys
import re
import argparse
import os

def verboseprint(verbose):
    '''
    Only print if verbose == True
    '''
    if verbose:
        def vprint(*args):
            # Print each argument separately 
            for arg in args:
                print(arg)
    else:
        def vprint(*args):
            None
    return vprint


def progressbar(it, prefix="", size=60, file=sys.stdout):
    '''
    dependency-free progressbar by:
    https://stackoverflow.com/users/1207193/eusoubrasileiro
    see this stackoverflow thread:
    https://stackoverflow.com/questions/3160699/python-progress-bar
    '''
    count = len(it)

    def show(j):
        x = int(size * j / count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i + 1)
    file.write("\n")
    file.flush()

def parse_ds9(region):
    """
    Parse ds9 region
    only has limited capabiliy for now.
    Will include astropy region pacakge or pyregion in the future.
    """
    # Skip blanks
    coordinate_systems = ['fk5', 'fk4', 'icrs', 'galactic', 'wcs', 
    'physical', 'image', 'ecliptic', 'J2000']
    outputlist = []
    with open(region) as fp:
        line = fp.readline()
        for line in fp:
            if line.strip() in coordinate_systems:
                print('DS9 region system : ' + line.strip())
                outputlist.append(line.strip())
            if line.strip()[0:6] == 'circle':
                print('Reading region ' + line.strip())
                outputlist.append(line.strip()[7:].split(','))
        return outputlist

class Parser(argparse.ArgumentParser):
    '''
    Handling errors
    '''
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        sys.stderr.write('use the -h option\" for helps\n')
        sys.exit(2)


def create_circular_mask(npixy, npixx, center=None, radius=None):
    #  use the middle of the image
    if center is None:
        center = [int(npixx / 2), int(npixy / 2)]
    #  use the smallest distance between the center and image walls
    if radius is None: 
        radius = min(center[0], center[1], npixx - center[0], npixy - center[1])

    Y, X = np.ogrid[:npixy, :npixx]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    mask = dist_from_center <= radius
    return mask

def check_caldb():
    try:
        caldb = os.environ['CALDB']
        return caldb
    except:
        sys.exit('$CALDB not found, aborting...')

def create_partial_circle_mask(npixy, npixx, center=None, radius=None, angle_min=None, angle_max=None, invert=False):
    #  use the middle of the image
    if center is None:
        center = [int(npixx / 2), int(npixy / 2)]
    #  use the smallest distance between the center and image walls
    if radius is None: 
        radius = min(center[0], center[1], npixx - center[0], npixy - center[1])
    Y, X = np.ogrid[:npixy, :npixx]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    angle_from_center = np.arctan((Y-center[1]) / (X - center[0])) * 57.2957795
    xyminus = np.where(((X - center[0]) * (Y*0 + 1) < 0) & ((Y - center[1]) * (X*0 + 1) < 0))
    angle_from_center[xyminus] =  angle_from_center[xyminus] + 180
    xmyp = np.where(((X - center[0]) * (Y*0 + 1) < 0) & ((Y - center[1]) * (X*0 + 1) >= 0))
    angle_from_center[xmyp] =  angle_from_center[xmyp] + 180
    xpym = np.where(((X - center[0]) * (Y*0 + 1) >= 0) & ((Y - center[1]) * (X*0 + 1) <= 0))
    angle_from_center[xpym] =  angle_from_center[xpym] + 360
    if center[0] >= int(npixx / 2):
        mask = (dist_from_center <= radius) & (angle_from_center >= angle_min) & ((angle_from_center <= angle_max))
    else:
        mask = (dist_from_center <= radius) & ((angle_from_center <= angle_min) | ((angle_from_center >= angle_max)))
    if invert:
        mask = (dist_from_center <= radius) & ((angle_from_center <= angle_min) | ((angle_from_center >= angle_max)))       
    return mask

def determin_angles(x0, y0, r0, x1, y1, r1):    
    d=np.sqrt((x1-x0)**2 + (y1-y0)**2)
    a=(r0**2 - r1**2 + d**2)/(2*d)
    h=np.sqrt(r0**2 - a**2)
    x2=x0+a*(x1-x0)/d   
    y2=y0+a*(y1-y0)/d   
    x3=x2+h*(y1-y0)/d 
    y3=y2-h*(x1-x0)/d
    x4=x2-h*(y1-y0)/d
    y4=y2+h*(x1-x0)/d
    ang1 = np.arctan((y3-y1) / (x3-x1)) * 57.2957795
    if ((y3-y1) < 0) and ((x3-x1) <0):
        ang1 += 180
    if ((y3-y1) > 0) and ((x3-x1) <0):
        ang1 += 180
    if ((y3-y1) < 0) and ((x3-x1) >0):
        ang1 += 360
    ang2 = np.arctan((y4-y1) / (x4-x1)) * 57.2957795
    if ((y4-y1) < 0) and ((x4-x1) <0):
        ang2 += 180
    if ((y4-y1) > 0) and ((x4-x1) <0):
        ang2 += 180
    if ((y4-y1) < 0) and ((x4-x1) >0):
        ang2 += 360    
    return ang1, ang2

