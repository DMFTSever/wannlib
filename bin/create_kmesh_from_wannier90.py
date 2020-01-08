import numpy as np
import argparse

def write_new_wannier90_win(fnameNew, mesh, fnameOld="wannier90.win"):
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

def generate_equidistant_kmesh(npoints):
    '''written such that the mesh is of monkhorst pack type as used in standard meshing in VASP'''
    ndim = len(npoints)
    meshvecs = []
    for n in npoints:
        if n == 1:
            meshvecs.append(np.array([0]))
        else:
            meshvecs.append(-np.linspace(start=-0.5,stop=0.5,num=n,endpoint=False)[::-1])
    return np.moveaxis(np.array(np.meshgrid(*meshvecs, indexing='ij')),source=0,destination=-1)

def check_naxdims(dim):
    try:
        dim = int(dim)
        if dim == 0:
            raise RuntimeError("dimension hast to be integer > 0!")
    except:
        raise
    
    return dim

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="This script creates a new wannier90 win file \
                                                  with a user specified k grid such that HrtoHk\
                                                  uses the new grid")
    parser.add_argument('newFileName', help='name of new .win file', type=str)
    parser.add_argument('naxdims', help='List of points along each axis', nargs='*', type=check_naxdims)
    parser.add_argument('--oldFileName', default="wannier90.win", help="Name of old wannier90 win \
                                                                        file. Defaults to wannier90.win",\
                        type=str)
    args = parser.parse_args()
    print(args.points)
    mesh = generate_equidistant_kmesh(args.naxdims)
    write_new_wannier90_win(fnameNew=args.newFileName, mesh=mesh)
    
