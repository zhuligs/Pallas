#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from zfunc import set_cell_from_vasp
from copy import deepcopy as cp


def normcell(cell1, cell2):
    lat1 = cell1.get_lattice()
    lat2 = cell2.get_lattice()
    pos1 = cell1.get_positions()
    pos2 = cell2.get_positions()
    nat = len(pos1)
    for i in range(nat):
        for j in range(3):
            if pos1[i][j] >= 1.0:
                pos1[i][j] -= int(pos1[i][j])
            if pos1[i][j] < 0:
                pos1[i][j] += (int(pos1[i][j]) + 1)
            if pos2[i][j] >= 1.0:
                pos2[i][j] -= int(pos2[i][j])
            if pos2[i][j] < 0:
                pos2[i][j] += (int(pos2[i][j]) + 1)

    normpos2 = []
    for i in range(nat):
        p1 = cp(pos1[i])
        p1xyz = np.dot(p1, lat1)
        p2 = cp(pos2[i])
        p2xyz = np.dot(p2, lat2)
        pdd = p2xyz - p1xyz
        print np.sqrt(np.dot(pdd, pdd))
        pp = np.zeros(3)
        ddd = []
        for ix in range(-1, 2):
            for iy in range(-1, 2):
                for iz in range(-1, 2):
                    pp[0] = p2[0] + float(ix)
                    pp[1] = p2[1] + float(iy)
                    pp[2] = p2[2] + float(iz)
                    ppxyz = np.dot(pp, lat2)
                    p12 = ppxyz - p1xyz
                    dist = np.sqrt(np.dot(p12, p12))
                    print dist, ix, iy, iz
                    ddd.append([dist, ppxyz])
        sddd = sorted(ddd, key=lambda x: x[0])
        print sddd[0]
        normpos2.append(sddd[0][1])

    normparray = np.array(normpos2)

    print normparray
    print '++++++++++++++'
    print np.dot(pos1, lat1)


def normcell2(cell1, cell2):
    pos1 = cell1.get_positions()
    pos2 = cell2.get_positions()
    nat = len(pos1)
    for i in range(nat):
        for j in range(3):
            if pos1[i][j] >= 1.0:
                pos1[i][j] -= int(pos1[i][j])
            if pos1[i][j] < 0:
                pos1[i][j] += (int(pos1[i][j]) + 1)
            if pos2[i][j] >= 1.0:
                pos2[i][j] -= int(pos2[i][j])
            if pos2[i][j] < 0:
                pos2[i][j] += (int(pos2[i][j]) + 1)

    normpos2 = []
    for i in range(nat):
        p1 = cp(pos1[i])
        # p1xyz = np.dot(p1, lat1)
        p2 = cp(pos2[i])
        # p2xyz = np.dot(p2, lat2)
        pdd = p2 - p1
        print np.sqrt(np.dot(pdd, pdd))
        pp = np.zeros(3)
        ddd = []
        for ix in range(-1, 2):
            for iy in range(-1, 2):
                for iz in range(-1, 2):
                    pp[0] = p2[0] + float(ix)
                    pp[1] = p2[1] + float(iy)
                    pp[2] = p2[2] + float(iz)
                    # ppxyz = np.dot(pp, lat2)
                    p12 = pp - p1
                    dist = np.sqrt(np.dot(p12, p12))
                    print dist, ix, iy, iz
                    ddd.append([dist, cp(pp)])
        sddd = sorted(ddd, key=lambda x: x[0])
        print sddd[0]
        normpos2.append(sddd[0][1])
    normparray = np.array(normpos2)
    print normparray


def main():
    x1 = set_cell_from_vasp('x1.vasp')
    x2 = set_cell_from_vasp('x2.vasp')
    normcell2(x1, x2)


if __name__ == "__main__":
    main()
