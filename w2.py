#!/usr/bin/python -u

import itin
import sdata
import fppy
from copy import deepcopy as cp
# import cPickle as pick
from wrapdimer import get_mode, get_0mode, get_rmode
from zfunc import gopt, rundim, set_cell_from_vasp, write_cell_to_vasp
import numpy as np


def ww(reac, prod):
    XMode = []
    XSaddle = []
    XLoc = []
    dist2p = []
    fpp = prod.get_lfp()
    vt = []
    for ip in range(itin.npop):
        # xmode = get_rmode()
        #
        # xsad = rundim(reac, xmode)
        vt.append(get_0mode())
        (xsad, xmode) = gen_rsaddle(reac)
        XMode.append(xmode)
        XSaddle.append(xsad)
        gmode = -1 * get_mode(xsad, prod)
        xloc = gopt(xsad, gmode)
        if abs(xloc.get_e() - 151206) < 1.0:
            gmode = get_0mode()
            xloc = gopt(xsad, gmode)
        XLoc.append(xloc)
        sdata.saddlehistory.append(xsad)
        sdata.localhistory.append(xloc)
        fp = xloc.get_lfp()
        (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fp, fpp)
        dist2p.append([ip, dist])
        print "ZLOG# NEW LOC %d DIST: %g" % (ip, dist)
    sortdist2p = sorted(dist2p, key=lambda x: x[1])
    iratio = int(itin.psoratio * itin.npop) + 1
    for ip in range(itin.npop):
        if dist2p[ip][1] > sortdist2p[iratio][1]:
            sdata.ifpso[ip] = False
        else:
            sdata.ifpso[ip] = True
    ipbest = sortdist2p[0][0]
    pmode = cp(XMode[ipbest])
    sdata.wvs.append(vt)
    sdata.pmodes.append(pmode)
    sdata.gmodes.append(pmode)
    sdata.wmodes.append(XMode)
    sdata.bestdist = sortdist2p[0][1]
    sdata.wsads.append(XSaddle)
    sdata.wlocs.append(XLoc)


def psov(xmode, v0, w, pbest, gbest):
    c1 = 2.0
    c2 = 2.0
    (r1, r2) = np.random.rand(2)
    v = v0 * w + c1 * r1 * (pbest - xmode) + c2 * r2 * (gbest - xmode)
    return v


def w3(reac, prod):
    # step 0
    ww(reac, prod)
    fpp = prod.get_lfp()
    for istep in range(itin.instep):
        xlocs = sdata.wlocs[istep]
        vt = []
        XMode = []
        XSad = []
        XLoc = []
        dist2p = []
        for ip in range(itin.npop):
            xloc = xlocs[ip]
            if sdata.ifpso[ip]:
                (xsad, xmode, v) = gen_psallde(xloc, istep, ip)
            else:
                (xsad, xmode) = gen_rsaddle(xloc)
                v = get_0mode()
                # not real mode, just generate the ZERO velocity
            vt.append(v)
            XSad.append(xsad)
            XMode.append(xmode)
            gmode = -1 * get_mode(xsad, prod)
            nloc = gopt(xsad, gmode)
            if abs(nloc.get_e() - 151206) < 1.0:
                gmode = get_0mode()
                nloc = gopt(xsad, gmode)
            XLoc.append(nloc)
            fp = nloc.get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fp, fpp)
            dist2p.append([ip, dist])
            print "ZLOG# STEP %d NEW LOC %d DIST: %g" % (istep, ip, dist)
        sortdist2p = sorted(dist2p, key=lambda x: x[1])
        iratio = itin.psoratio * itin.npop
        for ip in range(itin.npop):
            if dist2p[ip][1] > dist2p[iratio][1]:
                sdata.ifpso[ip] = False
            else:
                sdata.ifpso[ip] = True
        ipbest = sortdist2p[0][0]
        pbest = cp(XMode[ipbest])
        gbest = sdata.gmodes[istep]
        if sortdist2p[0][1] < sdata.bestdist:
            sdata.bestdist = sortdist2p[0][1]
            sdata.gmodes.append(pbest)
        else:
            sdata.gmodes.append(gbest)
        sdata.wvs.append(vt)
        sdata.pmodes.append(pbest)
        sdata.wmodes.append(XMode)
        sdata.wsads.append(XSad)
        sdata.wlocs.append(XLoc)

        if sortdist2p[0][1] < itin.dist:
            print "ZLOG # converged", sortdist2p[0][1]
            break


def gen_psallde(reac, istep, ip):
    xmode = cp(sdata.wmodes[istep][ip])
    pbest = sdata.pmodes[istep]
    gbest = sdata.gmodes[istep]
    v0 = sdata.wvs[istep][ip]
    w = 0.9 - 0.5 * (istep + 1) / itin.instep
    v = psov(xmode, v0, w, pbest, gbest)
    nmode = xmode + v
    xsad = rundim(reac, xmode)
    if abs(xsad.get_e() - 151206) < 1.0:
        (xsad, nmode) = gen_rsaddle(reac)
        v = get_0mode()
    return (xsad, nmode, v)


def gen_rsaddle(reac):
    succ = False
    while not succ:
        e = 151206.0
        while abs(e - 151206) < 1.0:
            xmode = get_rmode()
            xsad = rundim(reac, xmode)
            e = xsad.get_e()
            print 'ZLOG# SAD E:', e
        succ = checkident(xsad)
    return (xsad, xmode)


def checkident(xcell):
    # if the xcell is new, then return True
    if len(sdata.localhistory) == 0 or len(sdata.saddlehistory) == 0:
        return True
    xfp = xcell.get_lfp()
    for i, xstru in enumerate(sdata.saddlehistory):
        fp = xstru.get_lfp()
        (dist, mt) = fppy.fp_dist(itin.ntyp, sdata.types, xfp, fp)
        if dist < itin.dist:
            print 'ZLOG# same as sad %d, dist is %g' % (i, dist)
            return False
    for i, xstru in enumerate(sdata.localhistory):
        fp = xstru.get_lfp()
        (dist, mt) = fppy.fp_dist(itin.ntyp, sdata.types, xfp, fp)
        if dist < itin.dist:
            print 'ZLOG# same as loc %d, dist is %g' % (i, dist)
            return False
    return True


def initrun():
    sdata.ifpso = []
    for i in range(itin.npop):
        sdata.ifpso.append(True)
    reac0 = set_cell_from_vasp('R.vasp')
    prod0 = set_cell_from_vasp('P.vasp')
    mode = get_0mode()
    reac = gopt(reac0, mode)
    prod = gopt(prod0, mode)
    write_cell_to_vasp(reac, 'ROPT.vasp')
    write_cell_to_vasp(prod, 'POPT.vasp')
    sdata.types = reac.get_types()
    fpr = reac.get_lfp()
    fpp = prod.get_lfp()
    (d, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpr, fpp)
    print 'ZLOG# INIT DIST', d
    print 'ZLOG# REAC ENERGY', reac.get_e()
    print 'ZLOG# PROD ENERGY', prod.get_e()
    return (reac, prod)


def main():
    (reac, prod) = initrun()
    w3(reac, prod)


if __name__ == "__main__":
    main()
