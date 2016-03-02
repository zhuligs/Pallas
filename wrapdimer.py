#!/usr/bin/python

import itin
import fppy
from zfunc import *
from copy import deepcopy as cp
import numpy as np


def checkident(cell1, cell2):
    cell1.cal_fp(itin.fpcut, itin.lmax)
    cell2.cal_fp(itin.fpcut, itin.lmax)
    lfp1 = cell1.get_lfp()
    lfp2 = cell2.get_lfp()
    types = cell1.get_types()
    fpdist = fppy.fp_dist(itin.ntyp, types, lfp1, lfp2)
    if fpdist < 0.01:
        return True
    else:
        return False


def connect(reac, prod):
    Rpool = []
    Ppool = []
    for i in range(itin.ndimMax):
        print "## R DIM, ", i
        Rpool.append(rundim(reac))
        print "## P DIM, ", i
        Ppool.append(rundim(prod))

    DCOMPT = []
    for i in range(itin.ndimMax):
        xreac = Rpool[i]
        print 'get_lfp'
        fpi = np.array(xreac.get_lfp())
        print 'get_lfp done'
        for j in range(itin.ndimMax):
            xprod = Ppool[i]
            types = xprod.get_types()
            print 'get_lfp j'
            fpj = np.array(xprod.get_lfp())
            print 'get_lfp j done'
            dist = fppy.fp_dist(itin.ntyp, types, fpi, fpj)
            print 'get dist done'
            DCOMPT.append([dist, [i, j]])
    DCOMP = sorted(DCOMPT, key=lambda x: x[0])
    (ix, iy) = DCOMP[0][1]
    xsp = Rpool[ix]
    ysp = Ppool[iy]
    # mindist = DCOMP[0][0]

    # xspl : saddle point to local minima
    xspl = gopt(xsp)
    yspl = gopt(ysp)

    rbe = xsp.get_e() - reac.get_e()
    pbe = ysp.get_e() - prod.get_e()

    fpxs = np.array(xspl.get_lfp())
    fpys = np.array(yspl.get_lfp())
    d = fppy.fp_dist(itin.ntyp, types, fpxs, fpys)
    print "## D", d

    return (d, rbe, pbe, xsp, ysp, xspl, yspl)


def rconnect(xreac, xprod):
    dmax = 0.01
    d = 1.0
    RC = []
    ist = 0
    while d > dmax:
        ist += 1
        if ist > 100: break
        (d, rbe, pbe, xsp, ysp, xspl, yspl) = connect(xreac, xprod)
        xreac = cp(xspl)
        xprod = cp(yspl)
        RC.append([d, rbe, pbe, xsp, ysp, xspl, yspl])
        print "## IST", ist
        print '## DRP: ', d, rbe, pbe

    return RC


def main():
    reac = set_cell_from_vasp('R.vasp')
    xreac = gopt(reac)
    prod = set_cell_from_vasp('P.vasp')
    xprod = gopt(prod)
    RC = rconnect(xreac, xprod)
    i = 0
    for x in RC:
        i += 1
        print "no, d, rbe, pbe", i, x[0], x[1], x[2]
        write_cell_to_vasp(x[3], 'xsp' + str(i) + '.vasp')
        write_cell_to_vasp(x[4], 'yxp' + str(i) + '.vasp')
        write_cell_to_vasp(x[5], 'xspl' + str(i) + '.vasp')
        write_cell_to_vasp(x[6], 'yspl' + str(i) + '.vasp')


if __name__ == '__main__':
    main()














