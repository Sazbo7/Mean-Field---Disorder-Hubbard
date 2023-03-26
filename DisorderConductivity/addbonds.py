#!/usr/bin/env python

"""Paramagnetic Susceptibility
  
This code computes the paramagnetic current-gauge susceptibility of a disordered Hubbard model. Using a self-consistent solution of a Hartree-Fock-Bogoliubov equation of a disordered Hubbard model, current response of every bond to a uniform (and time-varying) gauge field is computed.

"""


__author__ = 'Kyungmin Lee'
__copyright__ = ''
__credits__ = ['Tim M. McCormick']
__version__ = '1.0.0'
__license__ = ''
__maintainer__ = 'Kyungmin Lee'
__email__ = 'kyungmin.lee.42@gmail.com'

import re
from collections import OrderedDict
from typing import List
import sys, os
import argparse
import numpy as np
import scipy as sp
import astropy.io.fits as pyfits
import datetime


def append_bonds(fitsfilepath: str):
  """Check whether the FITS file contains HoppingBonds extension, and append one if it does not. 

  HoppingBonds is a BinTableHDU with the following columns
  
  1. RowSite : int32
  2. ColSite : int32
  3. Amplitude : complex128
  4. DisplacementX : float64
  5. DisplacementY : float64
  """

  with pyfits.open(fitsfilepath, mode='append') as fitsfile:
    bondpattern = re.compile(r'HoppingBonds', re.I)
    if any(bondpattern.match(hdu.name) for hdu in fitsfile): return
 
    bonds = OrderedDict()
    prihdr = fitsfile['PRIMARY'].header
    a1x = prihdr['LATTICE BASIS1 X ANGS']
    a1y = prihdr['LATTICE BASIS1 Y ANGS']
    a2x = prihdr['LATTICE BASIS2 X ANGS']
    a2y = prihdr['LATTICE BASIS2 Y ANGS']
    b1  = prihdr['LATTICE NSITES B1']
    b2  = prihdr['LATTICE NSITES B2']
    t   = prihdr['HUBBARD MODEL T']

    nnmap = fitsfile['LatticeNearestNeighborsMap'].data
    sitemap = fitsfile['LatticeSiteOrdering'].data
    for (idx, (num_nn, *nns)) in enumerate(nnmap):
      for nn in nns:
        i, j = (idx+1, nn)
        if i > j: (i,j) = (j,i)
        bonds[(i,j)] = None
    bonds = np.array([k for k in bonds.keys()], dtype=np.int32)

    displacements = []
    for (i,j) in bonds:
      xi = np.array(sitemap[i-1])
      xj = np.array(sitemap[j-1])
      
      rij = xj - xi
      rij[0] = (rij[0] + b1//2) % b1 - b1//2
      rij[1] = (rij[1] + b2//2) % b2 - b2//2
      displacements.append(rij)
    displacements = np.array(displacements)

    rowsites = bonds[:,0]
    colsites = bonds[:,1]
    amplitudes = t * np.ones_like(colsites, dtype=np.complex128)
    displacement_x = a1x * displacements[:,0] + a2x * displacements[:,1]
    displacement_y = a1y * displacements[:,0] + a2y * displacements[:,1]

    # Format strings can be found in 
    # https://heasarc.gsfc.nasa.gov/docs/software/fitsio/quick/node10.html

    col1 = pyfits.Column(name='RowSite', format='J', unit='unitless', array=rowsites)
    col2 = pyfits.Column(name='ColSite', format='J', unit='unitless', array=colsites)
    col3 = pyfits.Column(name='Amplitude', format='M', unit='t', array=amplitudes)
    col4 = pyfits.Column(name='DisplacementX', format='D', unit='??', array=displacement_x)
    col5 = pyfits.Column(name='DisplacementY', format='D', unit='??', array=displacement_y)

    tbhdu = pyfits.BinTableHDU.from_columns(pyfits.ColDefs([col1, col2, col3, col4, col5]))
    tbhdu.header['EXTNAME'] = 'HoppingBonds'
    fitsfile.append(tbhdu)
    fitsfile.flush()


def main():
  parser = argparse.ArgumentParser('add list of bondsy')
  parser.add_argument('fitsfilepath', type=str, help='Input FITS file')
  args = parser.parse_args()

  append_bonds(args.fitsfilepath)

if __name__ == '__main__':
  main()

