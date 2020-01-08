#!/usr/bin/env python3

import sys
import os
sys.path.insert(0,sys.argv[0].replace('/wannlib/bin/make_new_wannier90_win_from_old.py',''))
from wannlib.core import wannlib
from wannlib.core import wannlib_core as core

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

parser = argparse.ArgumentParser(description="This script creates a new <winFile>.win file\
                                              from an old <wannier90>.win file replacing the kmesh\
                                              with a monkhorst pack mesh with the user specifed dimensions")
parser.add_argument('winFileNew', help='Name of the new <wannier90>.win file', type=str)
parser.add_argument('meshDims', nargs='*', help='Number of poitns per dimension', type=check_dims)
parser.add_argument('--winFileOld', default='wannier90.win', help='Name of the old <wannier90>.win file', type=str)
args = parser.parse_args()

mesh = core.generate_direct_coord_monkhorst_pack_kmesh(args.meshDims)
wannlib.make_new_wannier90_win_from_old(fnameNew=args.winFileNew, mesh=mesh, fnameOld=args.winFileOld)

