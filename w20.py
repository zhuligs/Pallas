#!/usr/bin/python -u
# encoding: utf-8

import sys
from copy import deepcopy as cp
import numpy as np

import itin
import sdata
import fppy
from wrapdimer import get_mode, get_0mode, get_rmode
from zfunc import gopt, rundim, set_cell_from_vasp, write_cell_to_vasp, getx
from w2 import gen_rsaddle, initrun


def wo(reac, prod):
    reace = reac.get_e()
    print 'ZLOG, reace', reace
    Xloc = []
    Yloc = []
    Xsad = []
    Ysad = []
    Xen = []
    Yen = []
    sdata.xllist.append(reac)
    sdata.yllist.append(prod)
    sdata.xslist.append(reac)
    sdata.yslist.append(prod)
    # DATABASE
    sdata.xmdb.append([[],[]]) # reac
    sdata.ymdb.append([[],[]]) # prod
    for ip in range(itin.npop):
        (xsad, vx) = gen_rsaddle(reac)
        (ysad, vy) = gen_rsaddle(prod)
        # update saddle point identical number
        xid = update_iden(sdata.xslist, xsad)
        yid = update_iden(sdata.yslist, ysad)
        xsad.set_iden(xid)
        ysad.set_iden(yid)
        xsad.set_sm('S')
        ysad.set_sm('S')
        xsad.add_nbor(reac.get_iden())
        ysad.add_nbor(prod.get_iden())
        gmod = get_0mode()
        xloc = gopt(xsad, gmod)
        yloc = gopt(ysad, gmod)
        xid = update_iden(sdata.xllist, xloc)
        yid = update_iden(sdata.yllist, yloc)
        xloc.set_iden(xid)
        yloc.set_iden(yid)
        xloc.set_sm('M')
        yloc.set_sm('M')
        xloc.add_nbor(xsad.get_iden())
        yloc.add_nbor(ysad.get_iden())
        print "ZLOG: INIT STEP, IP %4d X SAD EN: %8.7E, X LOC EN: %8.7E" %\
              (ip, xsad.get_e(), xloc.get_e())
        print "ZLOG: INIT STEP, IP %4d Y SAD EN: %8.7E, Y LOC EN: %8.7E" %\
              (ip, ysad.get_e(), yloc.get_e())
        Xsad.append(xsad)
        Ysad.append(ysad)
        sdata.vx.append(vx)
        sdata.vy.append(vy)
        Xloc.append(xloc)
        Yloc.append(yloc)
        Xen.append(xsad.get_e())
        Yen.append(ysad.get_e())
        sdata.xllist.append(xloc)
        sdata.yllist.append(yloc)
        sdata.xslist.append(xsad)
        sdata.yslist.append(ysad)
        # the init pbest
        sdata.pbestx.append(xloc)
        sdata.pbesty.append(yloc)
    sdata.xlocs.append(Xloc)
    sdata.ylocs.append(Yloc)


    xydist = []
    for ix in range(itin.npop):
        fpx = Xloc[ix].get_lfp()
        ex = Xsad[ix].get_e() - reace
        for iy in range(itin.npop):
            fpy = Yloc[iy].get_lfp()
            ey = Ysad[iy].get_e() - reace
            ee = max(ex, ey)
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            # mdist for multiobjective opt
            if dist < 1e-4:
                xdist = 1e-4
            else:
                xdist = dist
            mdist = np.log10(xdist) + ee
            print 'ZLOG: mdist, dist, log(dist), ee', mdist, dist, np.log(dist), ee
            xydist.append((mdist, (ix, iy), dist, ee))
    xydistSort = sorted(xydist,  key=lambda x: x[0])
    sdata.bestdist = xydistSort[0][2]
    sdata.bestmdist = xydistSort[0][0]
    sdata.gbestx = cp(Xloc[xydistSort[0][1][0]])
    sdata.gbesty = cp(Yloc[xydistSort[0][1][1]])
    print "ZLOG: INIT STEP, bestDist: %8.7E, bestmD: %8.7E, X-Y: %4d %4d" %\
          (xydistSort[0][2], xydistSort[0][0],
           xydistSort[0][1][0], xydistSort[0][1][1])

    # update pdist x
    for ix in range(itin.npop):
        xytdist = []
        fpx = Xloc[ix].get_lfp()
        ex = Xsad[ix].get_e() - reace
        for iy in range(len(sdata.yllist)):
            fpy = sdata.yllist[iy].get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            ey = sdata.yslist[iy].get_e() - reace
            ee = max(ex, ey)
            if dist < 1e-4:
                xdist = 1e-4
            else:
                xdist = dist
            mdist = np.log10(xdist) + ee
            print 'ZLOG: mdist, dist, log(dist), ee', mdist, dist, np.log(dist), ee
            xytdist.append(mdist)
        xytdistb = sorted(xytdist)[0]
        sdata.pdistx.append(xytdistb)

    # update pdist y
    for iy in range(itin.npop):
        yxtdist = []
        fpy = Yloc[iy].get_lfp()
        ey = Ysad[iy].get_e() - reace
        for ix in range(len(sdata.xllist)):
            fpx = sdata.xllist[ix].get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            ex = sdata.xslist[ix].get_e() - reace
            ee = max(ex, ey)
            if dist < 1e-4:
                xdist = 1e-4
            else:
                xdist = dist
            mdist = np.log10(xdist) + ee
            print 'ZLOG: mdist, dist, log(dist), ee',\
                  mdist, dist, np.log(dist), ee
            yxtdist.append(mdist)
        yxtdistb = sorted(yxtdist)[0]
        sdata.pdisty.append(yxtdistb)

    # iratio = int(itin.psoration * itin.npop)
    # if iratio >= itin.npop:
    #     iratio = itin.npop - 1
    for ip in range(itin.npop):
        sdata.ifpsox[ip] = True
        sdata.ifpsoy[ip] = True


def woo(reac, prod):
    # init step
    reace = reac.get_e()
    wo(reac, prod)
    # pso step
    for istep in range(itin.instep):
        Xloc = []
        Yloc = []
        Xsad = []
        Ysad = []
        xlocs = sdata.xlocs[istep]
        ylocs = sdata.ylocs[istep]
        for ip in range(itin.npop):
            xloc = xlocs[ip]
            yloc = ylocs[ip]
            if sdata.ifpsox[ip]:
                (xsad, vx) = gen_psaddle('x', xloc, istep, ip)
            else:
                (xsad, vx) = gen_rsaddle(xloc)
            if sdata.ifpsoy[ip]:
                (ysad, vy) = gen_psaddle('y', yloc, istep, ip)
            else:
                (ysad, vy) = gen_rsaddle(yloc)
            sdata.vx[ip] = cp(vx)
            sdata.vy[ip] = cp(vy)
            Xsad.append(xsad)
            Ysad.append(ysad)
            gmod = get_0mode()
            xxloc = gopt(xsad, gmod)
            yyloc = gopt(ysad, gmod)
            Xloc.append(xxloc)
            Yloc.append(yyloc)
            sdata.xllist.append(xxloc)
            sdata.yllist.append(yyloc)
            sdata.xslist.append(xsad)
            sdata.yslist.append(ysad)
            print "ZLOG: STEP %4d, IP %4d X SAD EN: %8.7E, X LOC EN: %8.7E" %\
                  (istep, ip, xsad.get_e(), xxloc.get_e())
            print "ZLOG: STEP %4d, IP %4d Y SAD EN: %8.7E, Y LOC EN: %8.7E" %\
                  (istep, ip, ysad.get_e(), yyloc.get_e())

        # update gbest
        xyldist = []
        for ix in range(len(sdata.xllist)):
            fpx = sdata.xllist[ix].get_lfp()
            ex = sdata.xslist[ix].get_e() - reace
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_lfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                ey = sdata.yslist[iy].get_e() - reace
                ee = max(ex, ey)
                if dist < 1e-4:
                    xdist = 1e-4
                else:
                    xdist = dist
                mdist = np.log10(xdist) + ee
                print 'ZLOG: mdist, dist, log(dist), ee', mdist, dist, np.log(dist), ee
                xyldist.append([mdist, [ix, iy], dist, ee])
        xyldistSort = sorted(xyldist, key=lambda x: x[0])
        ix = xyldistSort[0][1][0]
        iy = xyldistSort[0][1][1]
        sdata.gbestx = cp(sdata.xllist[ix])
        sdata.gbesty = cp(sdata.yllist[iy])
        bestdist = xyldistSort[0][0]
        # print 'ZLOG: STEP %4d, bestDist: %8.7E' % (istep, bestdist)
        print "ZLOG: STEP %4d, bestDist: %8.7E, X-Y: %d %d" %\
              (istep, bestdist, ix, iy)
        print "ZLOG: X %d SAD-E: %8.7E LOC-E: %8.7E" % \
              (ix, sdata.xslist[ix].get_e(), sdata.xllist[ix].get_e())
        print "ZLOG: Y %d SAD-E: %8.7E LOC-E: %8.7E" % \
              (iy, sdata.yslist[iy].get_e(), sdata.yllist[iy].get_e())

        write_de(xyldist, reac.get_e())

        # updata pbest
        xydist = []
        for ix in range(itin.npop):
            fpx = Xloc[ix].get_lfp()
            ex = Xsad[ix].get_e() - reace
            xytdist = []
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_lfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                ey = sdata.yslist[iy].get_e() - reace
                ee = max(ex, ey)
                if dist < 1e-4:
                    xdist = 1e-4
                else:
                    xdist = dist
                mdist = np.log10(xdist) + ee
                print 'ZLOG: mdist, dist, log(dist), ee', mdist, dist, np.log(dist), ee
                xytdist.append([mdist, iy])
            xytdistSort = sorted(xytdist, key=lambda x: x[0])
            xytbestdist = xytdistSort[0][0]
            # iy = xytbestdist[0][1]
            # get the best dist for each particle in the pop
            xydist.append(xytbestdist)
            if xytbestdist < sdata.pdistx[ix]:
                sdata.pbestx[ix] = cp(Xloc[ix])
                sdata.pdistx[ix] = xytbestdist

        yxdist = []
        for iy in range(itin.npop):
            fpy = Yloc[iy].get_lfp()
            ey = Ysad[iy].get_e() - reace
            yxtdist = []
            for ix in range(len(sdata.xllist)):
                fpx = sdata.xllist[ix].get_lfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                ex = sdata.xslist[ix].get_e() - reace
                ee = max(ex, ey)
                if dist < 1e-4:
                    xdist = 1e-4
                else:
                    xdist = dist
                mdist = np.log10(xdist) + ee
                print 'ZLOG: mdist, dist, log(dist), ee', mdist, dist, np.log(dist), ee
                yxtdist.append([mdist, ix])
            yxtdistSort = sorted(yxtdist, key=lambda x: x[0])
            yxtbestdist = yxtdistSort[0][0]
            yxdist.append(yxtbestdist)
            if yxtbestdist < sdata.pdisty[iy]:
                sdata.pbesty[iy] = cp(Yloc[iy])
                sdata.pdisty[iy] = yxtbestdist

        sdata.xlocs.append(Xloc)
        sdata.ylocs.append(Yloc)

        #if abs(bestdist) < itin.dist:
        #    print "ZLOG: CONVERGED!"
        #    break


def update_iden(xlist, cell):
    fpc = cell.get_lfp()
    oldids = []
    for x in xlist:
        fpx = x.get_lfp()
        idx = x.get_iden()
        (d, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpc)
        if d < itin.dist:
            idc = idx
            return idc
        oldids.append(idx)
    idc = max(oldids) + 1
    return idc


def gen_psaddle(xy, xcell, istep, ip):
    if xy is 'x':
        pbest = cp(sdata.pbestx[ip])
        gbest = cp(sdata.gbestx)
        v0 = sdata.vx[ip]
    elif xy is 'y':
        pbest = cp(sdata.pbesty[ip])
        gbest = cp(sdata.gbesty)
        v0 = sdata.vy[ip]
    else:
        print 'ERROR: xy'
        sys.exit(0)
    c1 = 2.0
    c2 = 2.0
    w = 0.9 - 0.5 * (istep + 1) / itin.instep
    (r1, r2) = np.random.rand(2)
    v = v0 * w + c1 * r1 * getx(pbest, xcell) + c2 * r2 * getx(gbest, xcell)
    scell = rundim(xcell, v)
    return (scell, v)


def update_path():
    # input: saddlelist, minimalist, npop, istep
    # reactant/product minimalist[0]





# def getpbest(xy, istep, ip):
#     if xy is 'x':
#         pbest = cp(sdata.pbestx[])
#     elif xy is 'y':
#     else:
#         print "ERROR: xy"
#         sys.exit(1)


def write_de(xyldist, e0):
    k = 0
    # for ix in range(len(sdata.xllist)):
    #     for iy in range(len(sdata.yllist)):
    #         e1 = sdata.xslist[ix].get
    #         de = min()
    data = []
    for xs in sdata.xslist:
        for ys in sdata.yslist:
            ex = xs.get_e() - e0
            ey = ys.get_e() - e0
            ee = max(ex, ey)
            md = xyldist[k][0]
            d = xyldist[k][2]
            k += 1
            data.append([d, md, ee])

    f = open('de.dat', 'w')
    for x in data:
        f.write("%8.7E  %8.7E  %8.7E\n" % tuple(x))
    f.close()


def main():
    (reac, prod) = initrun()
    woo(reac, prod)
    # outputw()


if __name__ == "__main__":
    main()

















