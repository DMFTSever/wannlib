import numpy as np
import datetime

def read_wannier90_win(fname='wannier90.win'):
    with open(fname) as instream:
        kpoints = []
        for line in instream:
            if ("unit_cell_cart" in line.split()) and ("begin" in line.split()):
                a = np.array(instream.readline().split(),dtype=float)
                b = np.array(instream.readline().split(),dtype=float) 
                c = np.array(instream.readline().split(),dtype=float)
                rbasis= np.array([a,b,c],dtype=float)
            if "mp_grid" in line.split():
                nk = int(np.array(line.split()[-3::],dtype=int).prod())
            if ("kpoints" in line.split()) and ("begin" in line.split()):
                kpoints = np.loadtxt(instream, max_rows=nk)
    return rbasis, kpoints

def read_wannier90_kpath(fname='wannier90.win'):
    with open(fname) as instream:
        read_kpath=False
        kpath = []
        for line in instream:
            splitline=line.split()
            if ("begin" in splitline) and ("kpoint_path" in splitline):
                read_kpath=True
            elif ("end" in splitline) and ("kpoint_path" in splitline):
                kpath.append(last)
                read_kpath=False
                break
            elif read_kpath==True:
                kpath.append(np.array(line.split()[1:4],dtype=float))
                last = np.array(line.split()[5:8],dtype=float)
    return kpath
            

def read_wannier90_hr(fname='wannier90_hr.dat'):
    with open(fname) as instream:
        comment_line = instream.readline()
        norb = int(instream.readline())
        nrpoints = int(instream.readline())
        weights = []
        for i in range(0,int(nrpoints/15)+1):
            weights += instream.readline().split()
        weights = np.array(weights, dtype=int)
        hr = np.loadtxt(instream)
    rpoints = np.array(hr[::norb*norb,0:3],dtype=float)
    bandsindices1= np.array(hr[:,3],dtype=int)-1
    bandsindices2= np.array(hr[:,4],dtype=int)-1
    rindices = np.repeat(np.arange(0,nrpoints),norb**2)
    hrMatrix = np.zeros((nrpoints,norb,norb),dtype=complex)
    hrMatrix[rindices,bandsindices1,bandsindices2]=hr[:,5]+1j*hr[:,6]

    #This is only for comparison to the straight forward poor mans approach (about 20-100 times slower)
    #hrnew2 = np.zeros((nrpoints,norb,norb),dtype=complex)
    #for i in range(0,nrpoints*norb**2):
    #    rint = int(i/(norb**2))
    #    hrnew2[rint,int(hr[i,3])-1,int(hr[i,4])-1]=hr[i,5]+1j*hr[i,6]
    #print(np.allclose(hrMatrix,hrnew2))
    
    return weights, rpoints, hrMatrix

def write_wannier90_hr(rpoints, weights, hrMatrix, fname='wannier90_hr.dat'):
    with open(fname, 'w') as outstream:
        outstream.write('written on '+str(datetime.datetime.now())+'\n') 
        for i in range(1,len(weights)+1):
            outstream.write('{:4d}'.format(weights[i-1])+' ')
            if i%15 == 0:
                outstream.write('\n')
            lasti = i
        if lasti%15 != 0:
            outstream.write('\n')
            
        for i in range(0, len(rpoints)):
            for b1 in range(0, hrMatrix.shape[1]):
                for b2 in range(0, hrMatrix.shape[2]):
                    outstream.write('{:4d}'.format(int(rpoints[i,0]))) 
                    outstream.write(' ')
                    outstream.write('{:4d}'.format(int(rpoints[i,1]))) 
                    outstream.write(' ')
                    outstream.write('{:4d}'.format(int(rpoints[i,2]))) 
                    outstream.write(' ')
                    outstream.write('{:4d}'.format(b1)) 
                    outstream.write(' ')
                    outstream.write('{:4d}'.format(b2)) 
                    outstream.write(' ')
                    outstream.write('{:12.8f}'.format(np.real(hrMatrix[i, b1, b2]))) 
                    outstream.write(' ')
                    outstream.write('{:12.8f}'.format(np.imag(hrMatrix[i, b1, b2]))) 
                    outstream.write('\n')

def write_wannier90_hk(fname, kpoints, hk):
    shape=hk.shape
    hk = np.concatenate((np.real(hk)[...,None], np.imag(hk)[...,None]), axis=-1).reshape(hk.shape[:-1]+(2*hk.shape[-1],))
    with open(fname,'w') as outstream:
        outstream.write('{a:8d}{b:8d}{c:8d}\n'.format(a=hk.shape[0],b=hk.shape[1],c=hk.shape[1]))
        for i, kpoint in enumerate(kpoints):
            outstream.write('{kx:13.10f} {ky:13.10f} {kz:13.10f}\n'.format(kx=kpoint[0], ky=kpoint[1], kz=kpoint[2]))
            np.savetxt(outstream, hk[i], fmt='%15.10f')

def generate_direct_coord_monkhorst_pack_kmesh(npoints):
    '''written such that the mesh is of monkhorst pack type as used in standard meshing in VASP'''
    ndim = len(npoints)
    meshvecs = []
    for n in npoints:
        if n == 1:
            meshvecs.append(np.array([0]))
        else:
            meshvecs.append(-np.linspace(start=-0.5,stop=0.5,num=n,endpoint=False)[::-1])
    return np.moveaxis(np.array(np.meshgrid(*meshvecs, indexing='ij')),source=0,destination=-1)

def generate_k_path(steps, points, kbasis=None):
        points = np.array(points)
        path = np.zeros((steps*(len(points)-1),3))
        cartpath = np.zeros((steps*(len(points)-1),3))
        if np.any(kbasis)==None:
            for point in range(0,len(points)-1):
                path[point*steps:(point+1)*steps,:] = np.linspace(points[point,:],points[point+1,:],steps)
            return path
        else:
            for point in range(0,len(points)-1):
                path[point*steps:(point+1)*steps,:] = np.linspace(points[point,:],points[point+1,:],steps)
                cartpath[point*steps:(point+1)*steps,:] = np.dot(np.linspace(points[point,:],\
                                                                             points[point+1,:],steps) \
                                                                 ,kbasis.transpose())
            return path, cartpath

def fourier_transform_hr(hr, kpoints, rpoints, weights=None):
    nr = len(rpoints)
    nk = len(kpoints)
    norb = hr.shape[-1]
    
    if nr*nk*norb**2*16 > 1000000000:
        print('To many kpoints for fast routine, switching to kpoint looping')
        hk = np.zeros((nk,norb,norb), dtype=complex)
        if np.any(weights) == None:
            for i, kpoint in enumerate(kpoints): 
                hk[i] = np.sum(hr[:,:,:] * np.exp(2*np.pi*1j*np.sum(rpoints[:,:] * \
                              kpoint[None,:],axis=1))[:,None,None], axis=0)
        else:
            for i, kpoint in enumerate(kpoints):
                hk[i] = np.sum(hr[:,:,:] * (np.exp(2*np.pi*1j*np.sum(rpoints[:,:] * \
                          kpoint[None,:],axis=1))/weights[:])[:,None,None], axis=0)
        return hk
        
    else:
        if np.any(weights) == None:
            return np.sum(hr[None,:,:,:] * np.exp(2*np.pi*1j*np.sum(rpoints[None,:,:] * \
                          kpoints[:,None,:],axis=2))[:,:,None,None], axis=1)
        else:
            return np.sum(hr[None,:,:,:] * (np.exp(2*np.pi*1j*np.sum(rpoints[None,:,:] * \
                          kpoints[:,None,:],axis=2))/weights[None,:])[:,:,None,None], axis=1)
