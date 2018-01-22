#!/usr/bin/env python

import numpy as np
from niggli import niggli_reduce


def get_niggli(lattice, pos):
    rlattice = np.transpose(niggli_reduce(np.transpose(lattice)))
    tranmat = np.matmul(rlattice, np.linalg.inv(lattice))
    # print(tranmat)

    rpos_tmp = np.dot(pos, np.linalg.inv(tranmat))
    rpos = rpos_tmp % 1.0
    # print(rpos)
    return (rlattice, rpos)


def refine(inposcar, outposcar):
    buff = []
    with open(inposcar) as f:
        for line in f:
            buff.append(line.split())

    lattice = np.array(buff[2:5], float)
    typt = np.array(buff[5], int)
    nat = sum(typt)
    pos = np.array(buff[7: 7 + nat], float)
    rlattice, rpos = get_niggli(lattice, pos)
    f = open('POSCAR_niggli.vasp', 'w')
    f.write(' '.join(buff[0]) + '\n')
    f.write('1.0\n')
    for x in rlattice:
        f.write("%15.9f %15.9f %15.9f\n" % tuple(x))
    f.write(' '.join(buff[5]) + '\n')
    f.write('Direct\n')
    for x in rpos:
        f.write("%15.9f %15.9f %15.9f\n" % tuple(x))
    f.close()


if __name__ == '__main__':
    refine('POSCAR', 'POSCAR_niggli.vasp')
