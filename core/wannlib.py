from . import wannlib_core as core
import numpy as np

def make_new_wannier90_win_from_old(fnameNew, mesh, fnameOld="wannier90.win"):
    '''
    Reads in old wannier90.win file and creates a new file fnameNew where the
    kmesh of the old file is replaced with the supplied kmesh
    '''
    ndim = mesh.shape[-1]
    naxdims = mesh.shape[0:-1]
    mesh = mesh.reshape(np.prod(naxdims),ndim)
    with open(fnameNew, 'w') as outstream:
        with open(fnameOld, 'r') as instream:
            for line in instream:
                if not 'mp_grid' in line.split():
                    outstream.write(line)
                else:
                    outstream.write("mp_grid = ")
                    for element in naxdims:
                        outstream.write(" {:6d} ".format(element))
                    outstream.write("\n\n")
                    outstream.write("begin kpoints\n")
                    np.savetxt(outstream, mesh,fmt="%17.12f")
                    outstream.write("end kpoints")
                    break

def make_wannier90_hk(fnameHk, meshDims, fnameHr='wannier90_hr.dat'):
    '''
    Reads in the supplied wannier90_hr file to build H(r). Then it fourier transforms 
    it to a monkhorst pack kmesh with the dimensions supplied with meshDims (nxpoints,
    nypoints, ...). Eventually the fourier transformed H(k) is written to fnameHK.
    '''
    weights, rpoints, hr = core.read_wannier90_hr(fname=fnameHr)
    kpoints = core.generate_direct_coord_monkhorst_pack_kmesh(npoints=meshDims).reshape(meshDims.prod(),3)
    hk = core.fourier_transform_hr(hr=hr, kpoints=kpoints, rpoints=rpoints, weights=weights)
    core.write_wannier90_hk(fname=fnameHk, hk=hk, kpoints=kpoints)

def make_wannier90_hk_on_DFTmesh(fnameHk, fnameWin='wannier90.win', fnameHr='wannier90_hr.dat'):
    '''
    Reads the supplied wannier_hr.dat and wannier.win file to build H(r) and a kmesh. Then
    it fourier transforms H(r) on the kmesh and writes the H(k) to fnameHk.
    '''
    weights, rpoints, hr = core.read_wannier90_hr(fname=fnameHr)
    temp, kpoints = core.read_wannier90_win(fname=fnameWin)
    hk = core.fourier_transform_hr(hr=hr, kpoints=kpoints, rpoints=rpoints, weights=weights)
    core.write_wannier90_hk(fname=fnameHk, hk=hk, kpoints=kpoints)

def generate_bandstructure_from_wannier90(steps, points, fnameHr='wannier90_hr.dat', fnameWin='wannier90.win'):
    '''
    Reads the supplied wannier_hr.dat and wannier.win file to build H(r) and the rbasis. Then 
    the kbasis is calculated and used to calculate the kpath in cartesian coordinates. After
    that it fourier transforms H(r) on the kpath and returns the a distance vector, the eigenvalues
    along the distance and the eigenvectors along the distance.
    '''
    rbasis, kpoints = core.read_wannier90_win(fname=fnameWin)
    kbasis = 2.*np.pi*np.linalg.inv(rbasis)
    weights, rpoints, hr = core.read_wannier90_hr(fname=fnameHr)
    kpath, cartkpath = core.generate_k_path(steps=steps, points=points, kbasis=kbasis)
    diffs = np.linalg.norm(cartkpath[1:]-cartkpath[0:-1],axis=1)
    dist = np.array([np.sum(diffs[:i]) for i in range(0,len(cartkpath))])
    ee, ev = np.linalg.eigh(core.fourier_transform_hr(hr=hr, kpoints=kpath, rpoints=rpoints, weights=weights))
    
    return dist, ee, ev

def make_bandstructure_from_wannier90(fnameOut, steps, points, fnameHr='wannier90_hr.dat', fnameWin='wannier90.win'):
    '''
    Uses generate_bandstructure_from_wannier90 to generate the bandstructure and than writtes it to a file
    formated such that it can be plottet with gnuplot. First colum is the band index, second is the distance
    along the kpath, third the eigenenergy and after that come the entries of the eigenvectors, i.e. the
    projections on the initial orbitals.
    '''
    distance, eigenenergies, eigenvectors = generate_bandstructure_from_wannier90(steps=steps, points=points, fnameHr=fnameHr, fnameWin=fnameWin)
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


#make_new_wannier90_win_from_old(fnameNew='wannier90_60x60x1.win', mesh=core.generate_direct_coord_monkhorst_pack_kmesh(npoints=[60,60,1]), fnameOld="wannier90.win")
#make_wannier90_hk(fnameHk='myHk.dat', meshDims=np.array([60,60,1]))
#make_wannier90_hk_on_DFTmesh(fnameHk='myHkDFT.dat', fnameWin='wannier90_60x60x1.win', fnameHr='wannier90_hr.dat')
#make_wannier90_hk_on_DFTmesh(fnameHk='test2.dat')
#make_bandstructure_from_wannier90(fnameOut='bands.dat', steps=100, points=np.array([[0,0,0],[0.333333,0.333333,0],[0.0,0.5,0],[-0.333333,0.666666,0],[0,0,0]]))

