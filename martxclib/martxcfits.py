"""
Functions for MSFC ART-XC
"""
from __future__ import print_function, division
import numpy as np


def rebinatt(atttime, attra, attdec, attres):
	'''
	The function to rebinning the attitude file
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
			newattdec[decid * k:(decid + 1) * k] = \
			np.linspace(attdec[decid], attdec[decid + 1], k)
		return newatttime, newattra, newattdec


def verboseprint(verbose):
	if verbose:
		def vprint(*args):
			# Print each argument separately so caller doesn't need to
			# stuff everything to be printed into a single string
			for arg in args:
				print(arg)
	else:
		def vprint(*args):
			None
	return vprint
