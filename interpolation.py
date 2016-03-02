#!/usr/bin/python -u

import sys
import numpy as np
from copy import deepcopy as cp
from itdbase import Cell, lat2lcons, lcons2lat


def interpolation(P1, P2):
    lat1 = P1.get_lattice()
    lat2 = P2.get_lattice()
    cons1 = lat2lcons(lat1)
    cons2 = lat2lcons(lat2)
    pos1 = P1.get_positions()
    pos2 = P2.get_positions()

    cons3 = (cons1 + cons2) / 2
    lat3 = lcons2lat(cons3)
    pos3 = (pos1 + pos2) / 2
    print pos1
    print pos2
    print pos3

    P3 = cp(P1)

    P3.set_lattice(lat3)
    P3.set_positions(pos3)

    return P3


def utest(vasp1, vasp2):
    P1 = Cell()
    P2 = Cell()
    buff = []
    with open(vasp1) as f:
        for line in f:
            buff.append(line.split())

    lat1 = np.array(buff[2:5], float)
    typt1 = np.array(buff[5], int)
    pos1 = np.array(buff[7:], float)
    P1.set_name('P1')
    P1.set_lattice(lat1)
    P1.set_positions(pos1)
    P1.set_typt(typt1)

    buff = []
    with open(vasp2) as f:
        for line in f:
            buff.append(line.split())

    lat2 = np.array(buff[2:5], float)
    typt2 = np.array(buff[5], int)
    pos2 = np.array(buff[7:], float)
    P2.set_name('P2')
    P2.set_lattice(lat2)
    P2.set_positions(pos2)
    P2.set_typt(typt2)

    P3 = interpolation(P1, P2)

    lat3 = P3.get_lattice()
    pos3 = P3.get_positions()

    f = open('P3.vasp', 'w')
    f.write('P3\n')
    f.write('1.0\n')
    for x in lat3:
        f.write("%15.9f %15.9f %15.9f\n" % tuple(x))
    for ix in typt2:
        f.write(str(ix) + ' ')
    f.write('\n')
    f.write('Direct\n')
    for x in pos3:
        f.write("%15.9f %15.9f %15.9f\n" % tuple(x))


if __name__ == '__main__':
    args = sys.argv
    vasp1 = args[1]
    vasp2 = args[2]
    utest(vasp1, vasp2)






