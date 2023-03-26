#!/usr/bin/env python

"""Diamagnetic Susceptibility

This code computes the diamagnetic current-gauge susceptibility of a disordered Hubbard model. Using a self-consistent solution of a Hartree-Fock-Bogoliubov equation of a disordered Hubbard model, current response of every bond to a uniform (and time-varying) gauge field is computed.

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
import git
import datetime
import cffi
import ctypes

from util import *
from addbonds import *
import c_susc


def compute_diamagnetic_susceptibility(
    temperature: float,
    #frequency: float, broadening: float,
    rowsites, colsites, amplitudes, displacements_x, displacements_y,
    eigenvalues, eigenvectors_up, eigenvectors_dn):
  """Wrapper function for computing diagmagnetic susceptibility
  """
  assert(temperature >= 0)
  #assert(frequency >= 0)
  #assert(broadening >= 0)
  assert(rowsites.shape == colsites.shape == amplitudes.shape == displacements_x.shape == displacements_y.shape)

  (num_bond, ) = rowsites.shape
  (num_eigen,) = eigenvalues.shape
  (num_eigen2, b1, b2) = eigenvectors_up.shape
  assert(num_eigen2 == num_eigen)
  assert((num_eigen2, b1, b2) == eigenvectors_dn.shape)
  # Jared's FITS file uses base-1
  rowsites        = np.ascontiguousarray(rowsites - 1, dtype=ctypes.c_int)
  colsites        = np.ascontiguousarray(colsites - 1, dtype=ctypes.c_int)
  amplitudes      = np.ascontiguousarray(amplitudes, dtype=np.complex128)
  displacements_x = np.ascontiguousarray(displacements_x, dtype=np.float64)
  displacements_y = np.ascontiguousarray(displacements_y, dtype=np.float64)
  eigenvalues     = np.ascontiguousarray(eigenvalues, dtype=np.float64)
  eigenvectors_up = np.ascontiguousarray(eigenvectors_up, dtype=np.complex128)
  eigenvectors_dn = np.ascontiguousarray(eigenvectors_dn, dtype=np.complex128)

  susceptibility_x_up = np.zeros(num_bond, dtype=np.complex128)
  susceptibility_x_dn = np.zeros(num_bond, dtype=np.complex128)
  susceptibility_y_up = np.zeros(num_bond, dtype=np.complex128)
  susceptibility_y_dn = np.zeros(num_bond, dtype=np.complex128)

  ffi = cffi.FFI()
  p_rowsites         = ffi.cast("int const *", rowsites.ctypes.data)
  p_colsites         = ffi.cast("int const *", colsites.ctypes.data)
  p_amplitudes       = ffi.cast("void const *", amplitudes.ctypes.data)
  p_displacements_x  = ffi.cast("double const *", displacements_x.ctypes.data)
  p_displacements_y  = ffi.cast("double const *", displacements_y.ctypes.data)
  p_eigenvalues      = ffi.cast("double const *", eigenvalues.ctypes.data)
  p_eigenvectors_up  = ffi.cast("void const *", eigenvectors_up.ctypes.data)
  p_eigenvectors_dn  = ffi.cast("void const *", eigenvectors_dn.ctypes.data)
  p_susceptibility_x_up = ffi.cast("void *", susceptibility_x_up.ctypes.data)
  p_susceptibility_x_dn = ffi.cast("void *", susceptibility_x_dn.ctypes.data)
  p_susceptibility_y_up = ffi.cast("void *", susceptibility_y_up.ctypes.data)
  p_susceptibility_y_dn = ffi.cast("void *", susceptibility_y_dn.ctypes.data)

  print("Computing diamagnetic susceptibility...")
  c_susc.lib.compute_diamagnetic_susceptibility(
      temperature, #frequency, broadening,
      num_bond,
      p_rowsites, p_colsites, p_amplitudes,
      p_displacements_x, p_displacements_y,
      num_eigen, b1, b2,
      p_eigenvalues,
      p_eigenvectors_up, p_eigenvectors_dn,
      p_susceptibility_x_up,
      p_susceptibility_x_dn,
      p_susceptibility_y_up,
      p_susceptibility_y_dn)
  print("Finished.")

  return (susceptibility_x_up, susceptibility_x_dn,
          susceptibility_y_up, susceptibility_y_dn)

def append_diamagnetic_susceptibility(fitsfilepath: str):
  """Compute the diamagnetic current-gauge susceptibility and append it to the FITS file.

  Diamagnetic susceptibility has no dependence on frequency.

  Parameters
  ----------
  fitsfilepath : str
      Path of FITS file.
  """
  with pyfits.open(fitsfilepath, mode='readonly') as fitsfile:
    diapattern = re.compile(r'DiamagneticSusceptibility', re.I)
    # Skip if result exists already
    for hdu in fitsfile:
      if diapattern.match(hdu.name): return 1

    # Read data
    hopping_bonds = np.ascontiguousarray(fitsfile['HoppingBonds'].data)
    rowsites = np.ascontiguousarray(hopping_bonds['RowSite'], dtype=ctypes.c_int)
    colsites = np.ascontiguousarray(hopping_bonds['ColSite'], dtype=ctypes.c_int)
    amplitudes = np.ascontiguousarray(hopping_bonds['Amplitude'], dtype=np.complex128)
    displacements_x = np.ascontiguousarray(hopping_bonds['DisplacementX'], dtype=np.float64)
    displacements_y = np.ascontiguousarray(hopping_bonds['DisplacementY'], dtype=np.float64)

    temperature = float(fitsfile['PRIMARY'].header['HUBBARD MODEL K_BT'])
    eigenvalues = np.ascontiguousarray(fitsfile['FinalEigenvalues'].data['Eigenvalues'], dtype=np.float64)
    eigenvectors_up = np.ascontiguousarray((fitsfile['FinalPsiUpEigenstates_Real'].data +
                                            1j * fitsfile['FinalPsiUpEigenstates_Imag'].data), dtype=np.complex128)
    eigenvectors_dn = np.ascontiguousarray((fitsfile['FinalPsiDownEigenstates_Real'].data +
                                            1j * fitsfile['FinalPsiDownEigenstates_Imag'].data), dtype=np.complex128)

  (diamagnetic_susceptibility_x_up,
   diamagnetic_susceptibility_x_dn,
   diamagnetic_susceptibility_y_up,
   diamagnetic_susceptibility_y_dn) = \
        compute_diamagnetic_susceptibility(
                        temperature,
                        #frequency, broadening,
                        rowsites, colsites, amplitudes, displacements_x, displacements_y,
                        eigenvalues,
                        eigenvectors_up, eigenvectors_dn)

  assert(diamagnetic_susceptibility_x_up.shape == (len(hopping_bonds), ) )
  assert(diamagnetic_susceptibility_x_dn.shape == (len(hopping_bonds), ) )
  assert(diamagnetic_susceptibility_y_up.shape == (len(hopping_bonds), ) )
  assert(diamagnetic_susceptibility_y_dn.shape == (len(hopping_bonds), ) )

  # Format strings can be found in
  # https://heasarc.gsfc.nasa.gov/docs/software/fitsio/quick/node10.html
  col1 = pyfits.Column(name='DiamagneticSusceptibilityXUp', format='M', array=diamagnetic_susceptibility_x_up)
  col2 = pyfits.Column(name='DiamagneticSusceptibilityXDn', format='M', array=diamagnetic_susceptibility_x_dn)
  col3 = pyfits.Column(name='DiamagneticSusceptibilityYUp', format='M', array=diamagnetic_susceptibility_y_up)
  col4 = pyfits.Column(name='DiamagneticSusceptibilityYDn', format='M', array=diamagnetic_susceptibility_y_dn)
  tbhdu = pyfits.BinTableHDU.from_columns( pyfits.ColDefs([col1, col2, col3, col4]) )

  timestamp = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M%p on %B %d, %Y")

  tbhdu.header.append(('HIERARCH PIPELINE CGS PROGRAM NAME', os.path.basename(sys.argv[0]), 'Name of simulation program' ))
  tbhdu.header.append(('HIERARCH PIPELINE CGS PROGRAM VERSION', __version__, 'Version of simulation program'))
  tbhdu.header.append(('HIERARCH PIPELINE CGS PROGRAM GIT VERSION', get_git_hash(), 'Git version hash of simulation program'))
  tbhdu.header.append(('HIERARCH PIPELINE CGS PROGRAM EXEC TIME', timestamp, 'Time (UTC)'))
  tbhdu.header['EXTNAME']= 'DiamagneticSusceptibility'

  pyfits.append(fitsfilepath, tbhdu.data, header=tbhdu.header)

  return 0

def main():
  parser = argparse.ArgumentParser('diamagnetic.py')
  parser.add_argument('fitsfilepath', type=str, help='Input FITS file')
  args = parser.parse_args()

  append_bonds(args.fitsfilepath)
  append_diamagnetic_susceptibility(args.fitsfilepath)

if __name__ == '__main__':
  main()
