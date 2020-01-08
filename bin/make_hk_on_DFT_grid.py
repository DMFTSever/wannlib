#!/usr/bin/env python3

import sys
import os
sys.path.insert(0,sys.argv[0].replace('/wannlib/bin/make_hk_on_DFT_grid.py',''))
from wannlib.core import wannlib

import argparse 
import numpy as np

parser = argparse.ArgumentParser(description="This script creates H(k) from a given \
                                              wannier90_hr on a k grid read from a \
                                              <wannier90>.win file and prints it to a file \
                                              formated for w2dynamics.")
parser.add_argument('HkFile', help='Name of the resulting H(k) file.', type=str)
parser.add_argument('--WannierWinFile', default='wannier90.win', \
                    help='Path of  the <wannier90>.win file containing the kpoints.', type=str)
parser.add_argument('--HrFile', default='wannier90_hr.dat', \
                    help='Path of the input <wannier90_hr>.dat file.', type=str)
args = parser.parse_args()
wannlib.make_wannier90_hk_on_DFTmesh(fnameHk=args.HkFile, fnameWin=args.WannierWinFile, fnameHr=args.HrFile)
