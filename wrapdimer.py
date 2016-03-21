#!/usr/bin/python

import itin
import fppy
from zfunc import *
from copy import deepcopy as cp
import numpy as np
from tsase.neb.util import vunit, vrand

# def checkident(cell1, cell2):
#     cell1.cal_fp(itin.fpcut, itin.lmax)
#     cell2.cal_fp(itin.fpcut, itin.lmax)
#     lfp1 = cell1.get_lfp()
#     lfp2 = cell2.get_lfp()
#     types = cell1.get_types()
#     fpdist = fppy.fp_dist(itin.ntyp, types, lfp1, lfp2)
#     if fpdist < 0.01:
#         return True
#     else:
#         return False


def get_mode(cell1, cell2):
    types = cell1.get_types()
    fp1 = np.array(cell1.get_lfp())
    fp2 = np.array(cell2.get_lfp())
    (d, m) = fppy.fp_dist(itin.ntyp, types, fp1, fp2)
    mode = np.zeros((itin.nat + 3, 3))
    mode[-3:] = cell2.get_lattice() - cell1.get_lattice()
    ilat = np.linalg.inv(cell1.get_lattice())
    vol = cell1.get_volume()
    jacob = (vol/itin.nat)**(1.0/3.0) * itin.nat**0.5
    mode[-3:] = np.dot(ilat, mode[-3:]) * jacob
    pos1 = cell1.get_cart_positions()
    pos2 = cell2.get_cart_positions()
    for i in range(itin.nat):
        mode[i] = pos2[m[i]] - pos1[i]
    # mode[:-3] = cell2.get_cart_positions() - cell1.get_cart_positions()
    mode = vunit(mode)
    return mode


def get_rmode():
    mode = np.zeros((itin.nat + 3, 3))
    mode = vrand(mode)
    mode = vunit(mode)
    return mode


def get_0mode():
    return np.zeros((itin.nat + 3, 3))


def con2(reac, prod):
    mode = get_mode(reac, prod)
    sdd = rundim(reac, mode)
    mode = -1.0 * get_mode(sdd, reac)
    newmini = gopt(sdd, mode)
    types = reac.get_types()
    fp0 = np.array(reac.get_lfp())
    fp1 = np.array(newmini.get_lfp())
    fp2 = np.array(prod.get_lfp())
    (d0, mt) = fppy.fp_dist(itin.ntyp, types, fp1, fp0)
    print "ZZ d0", d0
    (dist, mt) = fppy.fp_dist(itin.ntyp, types, fp1, fp2)
    return (dist, newmini)


def rcon2(xreac, xprod):
    dmax = 0.01
    d = 1.0
    ist = 0
    while d > dmax:
        ist += 1
        if ist > 100:
            break
        (d, newmini) = con2(xreac, xprod)
        xreac = cp(newmini)
        print 'ZZ# ist d', ist, d

    return 0


# def tc(xreac, xprod):
#     (d, newmini) = con2(xreac, xprod)
#     print 'ZZ# DD', d


def connect(reac, prod):

    Rpool = []
    Ppool = []
    for i in range(itin.ndimMax):
        print "ZZ# R DIM, ", i
        e = 151206.0
        while abs(e-151206) < 1.0:
            mode = get_rmode()
            tcc = cp(rundim(reac, mode))
            e = tcc.get_e()
            print 'ZZ# R DIM E', e
        Rpool.append(tcc)

        print "ZZ# P DIM, ", i
        e = 151206.0
        while abs(e-151206) < 1.0:
            mode = get_rmode()
            tcc = cp(rundim(prod, mode))
            e = tcc.get_e()
            print 'ZZ# P DIM E', e
        Ppool.append(tcc)

    DCOMPT = []
    for i in range(itin.ndimMax):
        xreac = cp(Rpool[i])
        # print 'ZZ# get_lfp'
        fpi = np.array(xreac.get_lfp())
        # print 'ZZ# get_lfp done'
        for j in range(itin.ndimMax):
            xprod = cp(Ppool[j])
            types = xprod.get_types()
            # print 'ZZ# get_lfp j'
            fpj = np.array(xprod.get_lfp())
            # print 'ZZ# get_lfp j done'
            (dist, mp) = fppy.fp_dist(itin.ntyp, types, fpi, fpj)
            print 'ZZ# get dist done', i, j, dist
            DCOMPT.append([dist, [i, j]])
    DCOMP = sorted(DCOMPT, key=lambda x: x[0])
    print 'ZZ# shortest D', DCOMP[0][0]
    (ix, iy) = DCOMP[0][1]
    print 'ZZ# ix iy', (ix, iy)
    xsp = cp(Rpool[ix])
    ysp = cp(Ppool[iy])
    fp1 = np.array(xsp.get_lfp())
    fp2 = np.array(ysp.get_lfp())
    (d1, m2) = fppy.fp_dist(itin.ntyp, types, fp1, fp2)
    print 'ZZ# CONF D', d1
    # # mindist = DCOMP[0][0]

    # xspl : saddle point to local minima
    modex = -1*get_mode(xsp, reac)
    modey = -1*get_mode(ysp, prod)
    # modex = get_rmode()
    modex = np.zeros((itin.nat + 3, 3))
    xspl = gopt(xsp, modex)
    ex = xspl.get_e()

    while abs(ex-151206) < 1.0:
        modex = get_rmode()
        xspl = gopt(xsp, modex)
        ex = xspl.get_e()
    # modey = get_rmode()
    modey = np.zeros((itin.nat + 3, 3))
    yspl = gopt(ysp, modey)
    ey = xspl.get_e()
    while abs(ey-151206) < 1.0:
        modey = get_rmode()
        yspl = gopt(ysp, modey)
        ey = yspl.get_e()

    rbe = xsp.get_e() - reac.get_e()
    pbe = ysp.get_e() - prod.get_e()

    fpxs = np.array(xspl.get_lfp())
    fpys = np.array(yspl.get_lfp())
    (d, mm) = fppy.fp_dist(itin.ntyp, types, fpxs, fpys)
    print "ZZ# D", d

    return (d, rbe, pbe, xsp, ysp, xspl, yspl)


def rconnect(xreac, xprod):
    dmax = itin.dist
    dd = 1.0
    RC = []
    ist = 0
    xfp0 = np.array(xreac.get_lfp())
    yfp0 = np.array(xprod.get_lfp())
    while dd > dmax:
        ist += 1
        if ist > 100:
            break
        (d, rbe, pbe, xsp, ysp, xspl, yspl) = connect(xreac, xprod)
        xreac = cp(xspl)
        xprod = cp(yspl)
        RC.append([d, rbe, pbe, xsp, ysp, xspl, yspl])

        dtt = []
        for i in range(len(RC)):
            xxl = RC[i][5]
            fpxxl = np.array(xxl.get_lfp())
            types = xxl.get_types()
            (dx, mm) = fppy.fp_dist(itin.ntyp, types, fpxxl, yfp0)
            dtt.append(dx)
            print 'ZZ I to PROD', i, dx
            for j in range(len(RC)):
                yyl = RC[j][6]
                fpyyl = np.array(yyl.get_lfp())
                (dy, mm) = fppy.fp_dist(itin.ntyp, types, fpyyl, xfp0)
                dtt.append(dy)
                print 'ZZ J to REAC', j, dy
                (dt, mm) = fppy.fp_dist(itin.ntyp, types, fpxxl, fpyyl)
                print "ZZ CON I J", i, j, dt
                dtt.append(dt)
        dd = min(dtt)
        print "ZZ# IST:", ist
        print "ZZ# DRP:", d, dd, rbe, pbe
        print "ZZ# X-S-E, X-L-E, Y-S-E, Y-L-E:", xsp.get_e(),\
            xspl.get_e(), ysp.get_e(), yspl.get_e()

    return RC


def main():
    mode = get_0mode()
    reac = set_cell_from_vasp('R.vasp')
    xreac = gopt(reac, mode)
    prod = set_cell_from_vasp('P.vasp')
    xprod = gopt(prod, mode)
    types = xreac.get_types()
    fp1 = np.array(xreac.get_lfp())
    fp2 = np.array(xprod.get_lfp())
    (d, m) = fppy.fp_dist(itin.ntyp, types, fp1, fp2)
    print 'ZZ INIT D', d
    print 'ZZ REAC-E', xreac.get_e()
    print 'ZZ PROD-E', xprod.get_e()
    # rcon2(xreac, xprod)
    RC = rconnect(xreac, xprod)
    i = 0
    for x in RC:
        i += 1
        print "ZZ# no, d, rbe, pbe", i, x[0], x[1], x[2]
        write_cell_to_vasp(x[3], 'xsp' + str(i) + '.vasp')
        write_cell_to_vasp(x[4], 'yxp' + str(i) + '.vasp')
        write_cell_to_vasp(x[5], 'xspl' + str(i) + '.vasp')
        write_cell_to_vasp(x[6], 'yspl' + str(i) + '.vasp')


if __name__ == '__main__':
    main()














