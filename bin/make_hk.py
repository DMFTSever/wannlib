#!/usr/bin/env python3

import sys
import os
sys.path.insert(0,sys.argv[0].replace('/wannlib/bin/make_hk.py',''))
from wannlib.core import wannlib

import argparse 
import numpy as np

def check_dims(dim):
    try:
        dim = int(dim)
        if (dim < 1):
            raise RuntimeError('meshDims has to be list of integers > 0')
    except:
        raise
    return dim

parser = argparse.ArgumentParser(description="This script creates H(k) from a given \
                                              wannier90_hr and prints it to a file \
                                              formated for w2dynamics.")
parser.add_argument('HkFile', help='Name of the resulting H(k) file.', type=str)
parser.add_argument('meshDims', nargs='*', help='Number of points per dimension', type=check_dims)
parser.add_argument('--HrFile', default='wannier90_hr.dat', \
                    help='Path of the input <wannier90_hr>.dat file.', type=str)
args = parser.parse_args()

wannlib.make_wannier90_hk(fnameHk=args.HkFile, meshDims=np.array(args.meshDims,dtype=int), fnameHr=args.HrFile)
