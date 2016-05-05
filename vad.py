#!/usr/bin/env python

# import numpy as np
import itin
import sdata
import fppy
from copy import deepcopy as cp
from wrapdimer import get_rmode, get_0mode, get_mode
from zfunc import set_cell_from_vasp, write_cell_to_vasp
from vfunc import runvdim, goptv

# def con(reac, prod):
# 	mode = get_mode(reac, prod)
# 	sdd = runvdim(reac, mode)
# 	mode = -1.0 * get_mode(sdd, reac)
# 	newmin = goptv(sdd, mode)
# 	types = sdata.types()
# 	fp0 = reac.get_lfp()


def con(reac, prod):
    rPool = []
    pPool = []
    for i in range(itin.ndimMax):
        print "ZLOG: R DIM", i
        mode = get_rmode()
        tcc = runvdim(reac, mode)
        rPool.append(tcc)
        print "ZLOG: P DIM", i
        mode = get_rmode()
        tcc = runvdim(prod, mode)
        pPool.append(tcc)

    dcompt = []
    for i in range(itin.ndimMax):
        xreac = cp(rPool[i])
        fpi = xreac.get_lfp()
        for j in range(itin.ndimMax):
            xprod = cp(pPool[j])
            fpj = xprod.get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpi, fpj)
            print "ZLOG: I %d J %d dist %8.6E" % (i, j, dist)
            dcompt.append([dist, [i, j]])
    dcomp = sorted(dcompt, key=lambda x: x[0])
    print "ZLOG: shortest dim D %8.6E" % (dcomp[0][0])
    (ix, iy) = dcomp[0][1]
    print "ZLOG: ix iy", ix, iy
    xsp = cp(rPool[ix])
    ysp = cp(pPool[iy])
    fp1 = xsp.get_lfp()
    fp2 = ysp.get_lfp()
    (d1, m1) = fppy.fp_dist(itin.ntyp, sdata.types, fp1, fp2)
    print "ZLOG: CONF D: %8.6E" % (d1)

    # modex = -1*get_mode(xsp, reac)
    # modey = -1*get_mode(ysp, prod)
    modex = get_0mode()
    xspl = goptv(xsp, modex)
    yspl = goptv(ysp, modex)

    rbe = xsp.get_e() - reac.get_e()
    pbe = ysp.get_e() - reac.get_e()
    fpxs = xspl.get_lfp()
    fpys = yspl.get_lfp()
    (d, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpxs, fpys)
    print "ZLOG: DD: ", d

    return(d, rbe, pbe, xsp, ysp, xspl, yspl)


def rcon(xreac, xprod):
    dmax = itin.dist
    dd = 1.0
    rc = []
    ist = 0
    xfp0 = xreac.get_lfp()
    yfp0 = xprod.get_lfp()
    while dd > dmax:
        ist += 1
        if ist > 200:
            break
        (d, rbe, pbe, xsp, ysp, xspl, yspl) = con(xreac, xprod)
        xreac = cp(xspl)
        xprod = cp(yspl)
        rc.append([d, rbe, pbe, xsp, ysp, xspl, yspl])

        dtt = []
        for i in range(len(rc)):
            xxl = rc[i][5]
            fpxxl = xxl.get_lfp()
            (dx, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpxxl, yfp0)
            dtt.append(dx)
            print "ZLOG: I %d to PROD dist %8.6E" % (i, dx)
            for j in range(len(rc)):
                yyl = rc[j][6]
                fpyyl = yyl.get_lfp()
                (dy, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpyyl, xfp0)
                dtt.append(dy)
                print "ZLOG: J %d to PROD dist %8.6E" % (j, dy)
                (dt, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpxxl, fpyyl)
                print "ZLOG: CONT I %2d J %2d dist %8.6E" % (i, j, dt)
                dtt.append(dt)
        dd = min(dtt)
        print "ZLOG: IST:", ist
        print "ZLOG: DRP:", d, dd, rbe, pbe
        print "ZLOG: X-S-E, X-L-E, Y-S-E, Y-L-E:", \
              xsp.get_e(), xspl.get_e(), ysp.get_e(), yspl.get_e()

    return rc


def initrun():
    reac0 = set_cell_from_vasp('R.vasp')
    prod0 = set_cell_from_vasp('P.vasp')
    mode = get_0mode()
    reac = goptv(reac0, mode)
    prod = goptv(prod0, mode)
    write_cell_to_vasp(reac, 'ROPT.vasp')
    write_cell_to_vasp(prod, 'POPT.vasp')
    sdata.types = reac.get_types()
    fpr = reac.get_lfp()
    fpp = prod.get_lfp()
    (d, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpr, fpp)
    print 'ZLOG: INIT DIST', d
    print 'ZLOG: REAC ENERGY', reac.get_e()
    print 'ZLOG: PROD ENERGY', prod.get_e()
    return (reac, prod)


def main():
    (reac, prod) = initrun()
    rc = rcon(reac, prod)
    i = 0
    for x in rc:
        i += 1
        print "ZZ# no, d, rbe, pbe", i, x[0], x[1], x[2]
        write_cell_to_vasp(x[3], 'xsp' + str(i) + '.vasp')
        write_cell_to_vasp(x[4], 'ysp' + str(i) + '.vasp')
        write_cell_to_vasp(x[5], 'xspl' + str(i) + '.vasp')
        write_cell_to_vasp(x[6], 'yspl' + str(i) + '.vasp')


if __name__ == '__main__':
    main()























