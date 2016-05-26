#!/usr/bin/python -u
# encoding: utf-8


from copy import deepcopy as cp
import itin
import sdata
import fppy
from wrapdimer import get_mode, get_0mode, get_rmode
from zfunc import gopt, rundim, set_cell_from_vasp, write_cell_to_vasp, getx
from w2 import gen_rsaddle


def wo(reac, prod):
    Xloc = []
    Yloc = []
    Xsad = []
    Ysad = []
    Xmod = []
    Ymod = []
    Xen = []
    Yen = []
    sdata.xllist.append(reac)
    sdata.yllist.append(prod)
    for ip in range(itin.npop):
        (xsad, xmod) = gen_rsaddle(reac)
        (ysad, ymod) = gen_rsaddle(prod)
        gmod = get_0mode()
        xloc = gopt(xsad, gmod)
        yloc = gopt(ysad, gmod)
        Xsad.append(xsad)
        Ysad.append(ysad)
        Xmod.append(xmod)
        Ymod.append(ymod)
        Xloc.append(xloc)
        Yloc.append(yloc)
        Xen.append(xsad.get_e())
        Yen.append(ysad.get_e())
        sdata.xllist.append(xloc)
        sdata.yllist.append(yloc)
        # the init pbest
        sdata.pbestx.append(xloc)
        sdata.pbesty.append(yloc)
    sdata.xlocs.append(Xloc)
    sdata.ylocs.append(Yloc)

    xydist = []
    for ix in range(itin.npop):
        fpx = Xloc[ix].get_lfp()
        for iy in range(itin.npop):
            fpy = Yloc[iy].get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            xydist.append((dist, (ix, iy)))
    xydistSort = sorted(xydist,  key=lambda x: x[0])
    sdata.bestdist = xydistSort[0][0]
    sdata.gbestx = cp(Xloc[xydistSort[0][1][0]])
    sdata.gbesty = cp(Yloc[xydistSort[0][1][1]])

    # update pdist x
    for ix in range(itin.npop):
        xytdist = []
        fpx = Xloc[ix].get_lfp()
        for iy in range(len(sdata.yllist)):
            fpy = sdata.yllist[iy].get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            xytdist.append(dist)
        xytdistb = sorted(xytdist)[0]
        sdata.pdistx.append(xytdistb)

    # update pdist y
    for iy in range(itin.npop):
        yxtdist = []
        fpy = Yloc[iy].get_lfp()
        for ix in range(len(sdata.xllist)):
            fpx = sdata.xllist[ix].get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            yxtdist.append(dist)
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
            if sdata.ifposy[ip]:
                (ysad, vy) = gen_psaddle('y', yloc, istep, ip)
            else:
                (ysad, vy) = gen_rsaddle(yloc)
            Xsad.append(xsad)
            Ysad.append(ysad)
            gmod = get_0mode()
            xxloc = gopt(xsad, gmod)
            yyloc = gopt(ysad, gmod)
            Xloc.append(xxloc)
            Yloc.append(yyloc)
            sdata.xllist.append(xxloc)
            sdata.yllist.append(yyloc)

        # update gbest
        xyldist = []
        for ix in range(len(sdata.xllist)):
            fpx = sdata.xllist[ix].get_lfp()
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_lfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                xyldist.append([dist, [ix, iy]])
        xyldistSort = sorted(xyldist, key=lambda x: x[0])
        ix = xyldistSort[0][1][0]
        iy = xyldistSort[0][1][1]
        sdata.gbestx = cp(sdata.xllist[ix])
        sdata.gbesty = cp(sdata.yllist[iy])
        bestdist = xyldistSort[0][0]
        print 'ZLOG: STEP %4d, bestDist: %8.7E' % (istep, bestdist)

        # updata pbest
        xydist = []
        for ix in range(itin.npop):
            fpx = Xloc[ix].get_lfp()
            xytdist = []
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_lfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                xytdist.append([dist, iy])
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
            yxtdist = []
            for ix in range(len(sdata.xllist)):
                fpx = sdata.xllist[ix].get_lfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                yxtdist.append([dist, ix])
            yxtdistSort = sorted(yxtdist, key=lambda x: x[0])
            yxtbestdist = yxtdistSort[0][0]
            yxdist.append(yxtbestdist)
            if yxtbestdist < sdata.pdisty[iy]:
                sdata.pbesty[iy] = cp(Yloc[iy])
                sdata.pdisty[iy] = yxtbestdist

        sdata.xlocs.append(Xloc)
        sdata.ylocs.append(Yloc)

        if bestdist < itin.dist:
            print "ZLOG: CONVERGED!"
            break


def gen_psaddle(xy, xcell, istep, ip):
    pass























