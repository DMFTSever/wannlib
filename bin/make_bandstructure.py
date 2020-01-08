#!/usr/bin/env python3

import sys
import os
sys.path.insert(0,sys.argv[0].replace('/wannlib/bin/make_bandstructure.py',''))
from wannlib.core import wannlib
from wannlib.core import wannlib_core as core

import argparse 
import numpy as np

def load_points(fname='kpath.dat'):
    points = np.loadtxt(fname, dtype=float)
    if points.shape[-1] != 3:
        raise RuntimeError('Kpoints have to have 3 coordinates (x,y,z)')
    return points

parser = argparse.ArgumentParser(description="This script creates the bandstructure from \
                                              a given wannier90_hr and prints it to a \
                                              file formated for gnuplot plotting.")
parser.add_argument('BandsFile', help='name of the resulting bands file.', type=str)
parser.add_argument('steps', help='Number of steps per interval', type=int)
parser.add_argument('--kpathFile', default='kpath.dat' , \
                    help="File containing the points of the kpath", type=str)
parser.add_argument('--fnameHr', default='wannier90_hr.dat', help="File containing the wannier90 \
                     hamitlonian in real space", type=str)
parser.add_argument('--fnameWin', default='wannier90.win', help="File containing the wannier90 \
                     input", type=str)
parser.add_argument('--kpathFromWannier', default=False, action='store_true', 
                    help="Use to read path from the wannier90.win file")
args = parser.parse_args()
if args.kpathFromWannier:
    points = core.read_wannier90_kpath(fname=args.fnameWin)
else:
    points = load_points(args.kpathFile)
wannlib.make_bandstructure_from_wannier90(fnameOut=args.BandsFile, steps=args.steps, points=points, fnameHr=args.fnameHr, fnameWin=args.fnameWin)
