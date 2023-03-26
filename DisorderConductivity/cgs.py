#!/usr/bin/env python

__author__ = 'Kyungmin Lee'
__copyright__ = ''
__credits__ = ['Tim M. McCormick']
__version__ = '1.0.0'
__license__ = ''
__maintainer__ = 'Kyungmin Lee'
__email__ = 'kyungmin.lee.42@gmail.com'

import argparse
from addbonds import *
from diamagnetic import *
from paramagnetic import *

def main():
  parser = argparse.ArgumentParser('cgs.py')
  parser.add_argument('fitsfilepath', type=str, help='Input FITS file')
  parser.add_argument('-f', '--frequency', type=float, nargs='+')
  parser.add_argument('-b', '--broadening', type=float, required=True)
  args = parser.parse_args()

  append_bonds(args.fitsfilepath)
  append_diamagnetic_susceptibility(args.fitsfilepath)
  for freq in args.frequency:
    append_paramagnetic_susceptibility(args.fitsfilepath, freq, args.broadening)

if __name__=='__main__':
  main()
