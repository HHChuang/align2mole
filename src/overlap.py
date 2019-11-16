#!/usr/bin/env python
#################################################################################################
#                                                                                               #
#   Program: Overlap two rigid bodies as much as possible.                                      #
#                                                                                               #
#   Input:                                                                                      #
#       1. structure to be varied; $1, extension: .xyz                                          #
#       2. reference structure; $2, extension: .xyz                                             #
#                                                                                               #
#   Output:                                                                                     #
#       1. Std-out: initial RMSD, final RMSD                                                    #
#       2. modified structure; shifted.$1                                                       #
#                                                                                               #
# History:                                                                                      #
# 2019/11/15, Grace, Thanks for the help form Dr. L.P. Wang, UCD.                               #
#   reference code: https://github.com/leeping/forcebalance/blob/master/src/molecule.py#L692    #
#                                                                                               #
#################################################################################################

import numpy as np
import re
import sys

def main():
    # 1. Import two structures from input arguments
    coord_var_atoms, coord_var = get_coord(str(sys.argv[1]))
    coord_ref_atoms, coord_ref = get_coord(str(sys.argv[2]))

    # 2. Calculate the initial RMSD
    rmsd1 = rmsd(coord_var, coord_ref)

    # 3. Move two structures to their centroid; translation
    trans = centroid(coord_ref)
    coord_var_cen = coord_var - centroid(coord_var)
    coord_ref_cen = coord_ref - centroid(coord_ref)
    # 4. Generate rotation matrix by Kabsch algorithm
    R = kabsch(coord_var_cen, coord_ref_cen)
    

    # 5. Rotate and translate
    coord_var_shifted = np.dot(coord_var_cen,R) + trans

    # 6. Export new structure into file named shifted.$1 (*.xyz)
    output_struc(coord_var_atoms,coord_var_shifted)

    # 7. Std-out initial and final RMSD
    rmsd2 = rmsd(coord_var_shifted, coord_ref)
    print("% 7.4f % 7.4f" % (rmsd1, rmsd2))

def get_coord(filename):
    f = open(filename, 'r')
    V = list()
    atoms = list()
    n_atoms = int(f.readline())

    f.readline() # comment line

    for lines_read, line in enumerate(f):

        if lines_read == n_atoms:
            break

        atom = re.findall(r'[a-zA-Z]+', line)[0]
        # atom = re.findall(r'[0-9]', line)[0]
        atom = atom.upper()

        numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', line)
        numbers = [float(number) for number in numbers]

        if len(numbers) == 3:
            V.append(np.array(numbers))
            atoms.append(atom)
        else:
            exit("Problem with reading input files".format(
                lines_read + 2))

    f.close()
    atoms = np.array(atoms)
    V = np.array(V)
    return atoms, V


def rmsd(A, B):
    Coord = len(A[0])
    NAtom = len(A)
    cum = 0.0
    for i in range(NAtom):
        for j in range(Coord):
            cum += (A[i][j] - B[i][j])**2.0
    return np.sqrt(cum / NAtom)


def centroid(A):
    A = A.mean(axis=0)
    return A


def kabsch(coord_var, coord_ref):

    # covariance matrix
    covar = np.dot(coord_var.T, coord_ref)
    print(covar)

    # SVD  http://en.wikipedia.org/wiki/Kabsch_algorithm
    v, s, wt = np.linalg.svd(covar)

    # Transposition of v,wt
    d = (np.linalg.det(v) * np.linalg.det(wt)) < 0.0
    # right-hand coord
    if d:
        s[-1] = -s[-1]
        v[:, -1] = -v[:, -1]

    # Create Rotation matrix R
    R = np.dot(v, wt)

    return R

def output_struc(atoms,coord):
    NAtoms = len(atoms)
    filename = 'shifted.' + str(sys.argv[1])
    shifted = open(filename,'w')
    shifted.write(str(NAtoms) + '\n')
    shifted.write(filename + '\n')
    for j in range(NAtoms):
            line = str(atoms[j]) + ' ' + str(coord[j][0]) +\
                ' ' + str(coord[j][1]) + ' ' + str(coord[j][2]) + '\n'
            shifted.write(line)
    shifted.close()

main()