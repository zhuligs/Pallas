#!/usr/bin/python -u
# encoding: utf-8

import itin
import sdata
import fppy
from copy import deepcopy as cp
# import cPickle as pick
from wrapdimer import get_mode, get_0mode, get_rmode
from zfunc import gopt, rundim, set_cell_from_vasp, write_cell_to_vasp
import numpy as np
import os
from tsase.neb.util import vunit

from w2 import initrun


def goo(reac, prod):
    fpp = prod.get_lfp()
    xlocs = []
    xmods = []
    dist2p = []
    x0s = []
    for ip in range(itin.npop):
        (xloc, xloc0, xmod) = get_rloc(reac)
        xlocs.append(xloc)
        xmods.append(xmod)
        sdata.localhistory.append(xloc)
        fpx = xloc.get_lfp()
        (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpp)
        dist2p.append([ip, dist])
        x0s.append(xloc0)
        print "ZLOG: STEP %d NEW LOC %d DIST: %8.7E" % (-1, ip, dist)
    sortdist2p = sorted(dist2p, key=lambda x: x[1])
    iratio = int(itin.psoratio * itin.npop)
    if iratio >= itin.npop:
        iratio = itin.npop - 1
    if iratio < 0:
        iratio = 0
    for ip in range(itin.npop):
        if dist2p[ip][1] > sortdist2p[iratio][1]:
            sdata.ifpso[ip] = False
        else:
            sdata.ifpso[ip] = True
    ipbest = sortdist2p[0][0]
    sdata.bestdist = sortdist2p[0][1]
    gbest = cp(xlocs[ipbest])
    sdata.wlocs.append(xlocs)
    sdata.wmodes.append(xmods)
    sdata.gbests.append(gbest)
    sdata.wdis.append(dist2p)
    sdata.x0s.append(x0s)


def goopso(reac, prod):
    goo(reac, prod)
    fpp = prod.get_lfp()
    for istep in range(itin.instep):
        xlocs = []
        x0s = []
        xmods = []
        dist2p = []
        for ip in range(itin.npop):
            if sdata.ifpso[ip]:
                (ploc, xloc0, xmod) = get_ploc(istep, ip)
            else:
                # (ploc, xloc0, xmod) = get_rloc(sdata.wlocs[istep][ip])
                (ploc, xloc0, xmod) = get_rloc(reac)
            xlocs.append(ploc)
            x0s.append(xloc0)
            xmods.append(xmod)
            fpx = ploc.get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpp)
            dist2p.append([ip, dist])
            print "ZLOG: STEP %d NEW LOC %d DIST: %8.7E" % (istep, ip, dist)
        sortdist2p = sorted(dist2p, key=lambda x: x[1])
        iratio = int(itin.psoratio * itin.npop)
        if iratio >= itin.npop:
            iratio = itin.npop - 1
        if iratio < 0:
            iratio = 0
        for ip in range(itin.npop):
            if dist2p[ip][1] > sortdist2p[iratio][1]:
                sdata.ifpso[ip] = False
            else:
                sdata.ifpso[ip] = True
        ipbest = sortdist2p[0][0]
        if sdata.bestdist > sortdist2p[0][1]:
            sdata.bestdist = sortdist2p[0][1]
            gbest = cp(xlocs[ipbest])
        else:
            gbest = cp(sdata.gbests[-1])
        sdata.gbests.append(gbest)
        sdata.wlocs.append(xlocs)
        sdata.wmodes.append(xmods)
        sdata.wdis.append(dist2p)
        sdata.x0s.append(x0s)

        if sortdist2p[0][1] < itin.dist:
            print "ZLOG: converged", sortdist2p[0][1]
            break


def get_ploc(istep, ip):
    x0 = sdata.x0s[istep][ip]
    pbest = sdata.wlocs[istep][ip]
    gbest = sdata.gbests[istep]
    v0 = sdata.wmodes[istep][ip]
    w = 0.9 - 0.5 * (istep + 1) / itin.instep
    c1 = 2.0
    c2 = 2.0
    (r1, r2) = np.random.rand(2)
    v = v0 * w + c1 * r1 * getx(pbest, x0) + c2 * r2 * getx(gbest, x0)
    # v = v0 * w + c1 * r1 * getx(x0, pbest) + c2 * r2 * getx(x0, gbest)
    # v = vunit(v)
    xopt = applymode(pbest, v)
    mode0 = get_0mode()
    ploc = gopt(xopt, mode0)
    return (ploc, xopt, v)


def get_rloc(reac):
    succ = False
    while not succ:
        e = 151206.0
        while abs(e - 151206) < 1.0:
            rmode = get_rmode()
            xopt = applymode(reac, rmode)
            mode0 = get_0mode()
            xloc = gopt(xopt, mode0)
            e = xloc.get_e()
        succ = checkloc(xloc)
    return (xloc, xopt, rmode)


def applymode(xcell, mode):
    vol = xcell.get_volume()
    jacob = (vol / itin.nat)**(1.0/3.0) * itin.nat**0.5
    lat = xcell.get_lattice() + np.dot(xcell.get_lattice(), mode[-3:]/jacob)
    xcell.set_lattice(lat)
    xcell.set_cart_positions(xcell.get_cart_positions() + mode[:-3])
    return xcell


def getx(cell1, cell2):
    mode = np.zeros((itin.nat + 3, 3))
    mode[-3:] = cell1.get_lattice() - cell2.get_lattice()
    ilat = np.linalg.inv(cell1.get_lattice())
    vol = cell1.get_volume()
    jacob = (vol / itin.nat)**(1.0 / 3.0) * itin.nat**0.5
    mode[-3:] = np.dot(ilat, mode[-3:]) * jacob
    pos1 = cell1.get_cart_positions()
    pos2 = cell2.get_cart_positions()
    for i in range(itin.nat):
        mode[i] = pos1[i] - pos2[i]
    mode = vunit(mode)
    return mode


def checkloc(xloc):
    if len(sdata.localhistory) == 0:
        return True
    xfp = xloc.get_lfp()
    for i, xstru in enumerate(sdata.localhistory):
        fp = xstru.get_lfp()
        (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, xfp, fp)
        if dist < itin.dist:
            print "ZLOG# same as loc %d, dist is %8.7E" % (i, dist)
            return False
    return True


def outputw2():
    f = open('ZOUT', 'w')
    n = len(sdata.wdis)
    for i in range(n):
        f.write("STEP: %d\n" % i)
        pdir = 'data' + str(i)
        os.system('mkdir ' + pdir)
        for j in range(itin.npop):
            xloc = sdata.wlocs[i][j]
            dis = sdata.wdis[i][j][1]
            write_cell_to_vasp(xloc, pdir + '/xloc_' + str(j) + '.vasp')
            f.write("%4d  %8.7E  %8.7E\n" %
                    (j, xloc.get_e(), dis))
    f.close()


def main():
    (reac, prod) = initrun()
    goopso(reac, prod)
    outputw2()


if __name__ == "__main__":
    main()

