"""
Functions for MSFC ART-XC
"""
from __future__ import print_function, division
import numpy as np
import sys
import re

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
    coordinate_systems = ['fk5', 'fk4', 'icrs', 'galactic', 'wcs', 'physical', 'image', 'ecliptic', 'J2000']
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
