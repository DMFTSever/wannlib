###### NEEDS TO BE TESTED #########
###### NEEDS TO BE TESTED #########
###### NEEDS TO BE TESTED #########
###### NEEDS TO BE TESTED #########
###### NEEDS TO BE TESTED #########
###### NEEDS TO BE TESTED #########
import numpy as np
import sys
import wannlib.core as core

def read_fplo(fname):
    '''
    Routine to read in the rbasis and hr from and fplo calculation.
    Note that this produces output formated different than the read
    routines for wannier90 hamiltonians. In the output format of hr
    is as follows:
    [band1, band2, rpoints, 5]
    where in the last dimension the first 3 entries are the coordinates
    of the R vectors of the lattice, the 4 entry is the real and the 5
    entry is the imaginary part of Hij(R). Since R is real we store
    real and imaginary part seperatly.
    '''
    norb = 0
    nrpoints = 10
    hr = np.zeros((norb,norb,nrpoints,5))
    for i in range(0,20):
        try:
            with open(fname, 'r') as instream:
                read = False
                for line in instream:
                    splitline = line.split()
                    if ('Tij,' in splitline) and ('Hij:' in splitline) and ('end' in splitline):
                        read = False
                    elif ('Tij,' in splitline) and ('Hij:' in splitline):
                        read = True
                        i, j = instream.readline().split()
                        i = int(i)-1
                        j = int(j)-1
                        nr = 0
                    elif read == True:
                        hr[i,j,nr] = splitline
                        nr += 1
                    elif 'nwan:' in splitline:
                        norb = int(instream.readline())
                        hr = np.zeros((norb,norb,nrpoints,5), dtype=float)
                    elif 'lattice_vectors:' in splitline:
                        rbasis = np.loadtxt(instream, max_rows=3)
        except IndexError:
            nrpoints = 2*nrpoints
            continue
        break
    # NOTE UNSORTED HR (Needs to be resorted such that spin index is varyin slowest!)
    return rbasis, hr 

def fourier_transfrom_hr(hr, kpoints):
    '''
    hr file formated as it is outputed from read_fplo. Therefore kpoints
    have to be cartesian, since fplo prints the R points in cartesian coord.
    '''
    norb = hr.shape[0]
    nk = len(kpoints)
    hk = np.zeros((norb,norb,nk), dtype=complex)

    for i in range(0, norb):
        for j in range(0, norb):
            hk[i,j]=np.sum(np.exp(np.sum(kpoints[:,None,:]*hr[i,j,None,:,0:3], axis=2))\
                    *(hr[i,j,None,:,3]+1j*hr[i,j,None,:,4]),axis=1)

    # follwing is the real numpy way of writting the fourier transform. However in this
    # case the exponential has to be calculated for a lot more values effectively slowing
    # down the calculation so much that for 64 bands the upper routine is faster!
    #
    # hknew=np.sum(np.exp(np.sum(kpoints[None,None,:,None,:]*hr[:,:,None,:,0:3], axis=4))\
    #       *(hr[:,:,None,:,3]+1j*hr[:,:,None,:,4]),axis=3)
    
    return hk.transpose(2,0,1)
            

if __name__ == "__main__":
    sys.path.insert(0,sys.path[0].replace('fplo','core'))
    import wannlib_core as core

    rbasis, hr = read_fplo(fname='hamdata')
    kbasis = 2.*np.pi*np.linalg.inv(rbasis)
    path, cartkpath = core.generate_k_path(steps=100, points=np.array([[0,0,0],[0.5,0.5,0.5]]), kbasis=kbasis)
    hk = fourier_transfrom_hr(hr=hr, kpoints=cartkpath)
