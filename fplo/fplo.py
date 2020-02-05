#########NEEDS TO BE TESTED##########
#########NEEDS TO BE TESTED##########
#########NEEDS TO BE TESTED##########
#########NEEDS TO BE TESTED##########
#########NEEDS TO BE TESTED##########
#########NEEDS TO BE TESTED########## #########NEEDS TO BE TESTED##########
#########NEEDS TO BE TESTED##########
#########NOT WORKING##########

import sys
import os
sys.path.insert(0,sys.argv[0].replace('/wannlib/fplo/fplo.py',''))

import numpy as np
#from . import fplo_core as fplo
from wannlib.fplo import fplo_core as fplo
import wannlib.core.wannlib_core as core

def make_wannier90_hk(fnameHk, meshDims, fnameHr):
    '''
    Routine to read fplo H(r) and Fourier transfrom it to
    a monhorst pack kmesh with user specifed dimensions
    '''
    rbasis, hr = fplo.read_fplo(fname=fnameHr)
    mesh = core.generate_direct_coord_monkhorst_pack_kmesh(npoints=meshDims).reshape(meshDims.prod(),3)
    hk = fplo.fourier_transform_hr(hr=hr,kpoints=mesh)
    core.write_wannier90_hk(fname=fnameHk, hk=hk, kpoints=mesh)

def generate_bandstructure_from_fplo(steps, points, fnameHr):
    #TODO GO ON HERE fix it!!!
    '''
    Generating a Bandstructure form FPLO outpout.
    '''
    rbasis, hr = fplo.read_fplo(fname=fnameHr)
    print(rbasis.shape, hr.shape)
    print(rbasis)
    kbasis = 2.*np.pi*np.linalg.inv(rbasis)
    kpath, cartkpath = core.generate_k_path(steps=steps, points=points, kbasis=kbasis)
    diffs = np.linalg.norm(cartkpath[1:]-cartkpath[0:-1],axis=1)
    dist = np.array([np.sum(diffs[:i]) for i in range(0,len(cartkpath))])
    hk = fplo.fourier_transform_hr(hr=hr, kpoints=cartkpath)
    print(cartkpath.shape, hk.shape)
    ee, ev = np.linalg.eigh(hk)
    #IT ALREADY WRONG HERE!
    print(ee.shape)
    np.savetxt('test.dat', ee)

    return dist, ee, ev

def make_bandstructure_from_fplo(fnameOut, steps, points, fnameHr):
    distance, eigenenergies, eigenvectors = generate_bandstructure_from_fplo(steps=steps, points=points, fnameHr=fnameHr)
    nintervals=len(points)-1
    norb = eigenenergies.shape[-1]
    nenergyvec = np.repeat(np.arange(1,norb+1),steps*nintervals)
    eigenenergies = eigenenergies.flatten(order='F')
    eigenvectors = eigenvectors.transpose(1,0,2).reshape((eigenvectors.shape[0]*eigenvectors.shape[1],eigenvectors.shape[2]))
    projections = np.real(eigenvectors*eigenvectors.conjugate())
    distance = np.repeat(distance[:,None], repeats=norb, axis=1).flatten(order='F')
    bands = np.concatenate((nenergyvec[:,None], distance[:,None], eigenenergies[:,None], projections), axis=1)
    formatspec = ['%6d','%12.8f','%12.8f']
    for i in range(norb):
        formatspec += ['%12.8f']
    with open(fnameOut,'w') as outstream:
        outstream.write("#PATH\n#")
        for point in points[:-1]:
            outstream.write(" {} -".format(point))
        outstream.write(" {}".format(points[-1]))
        outstream.write("\n\n")
        outstream.write('{index:^6} {kdist:^12} {energy:^12}'.format(index='#index' ,kdist='kdist', energy='energy'))
        for i in range(norb):
            outstream.write(' {orb:^12}'.format(orb='orb'+str((i+1))))
        outstream.write('\n')
        for i in range(0,norb):
            for j in range(0,nintervals):
                np.savetxt(fname=outstream, X=bands[(i*nintervals+j)*steps:(i*nintervals+j+1)*steps], fmt=formatspec)
                outstream.write('\n')
            outstream.write('\n')


if __name__ == '__main__':
    import sys
    import os
    sys.path.insert(0,sys.argv[0].replace('/wannlib/fplo/fplo.py',''))
    from wannlib.core import wannlib
    from wannlib.core import wannlib_core as core
    make_bandstructure_from_fplo('bands.dat', 50, np.array([[0,0,0],[0.5,0.5,0]]),'hamdata')
