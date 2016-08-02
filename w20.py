#!/usr/bin/python -u
# encoding: utf-8

import sys
from copy import deepcopy as cp
import numpy as np
import cPickle as pick

from treelib import Node, Tree

import itin
import sdata
import fppy
from wrapdimer import get_mode, get_0mode, get_rmode
from zfunc import gopt, rundim, set_cell_from_vasp, write_cell_to_vasp, getx
from w2 import gen_rsaddle, initrun

sys.setrecursionlimit(100000)

def wo(reac, prod):
    reace = reac.get_e()
    print 'ZLOG, reace', reace
    Xloc = []
    Yloc = []
    Xsad = []
    Ysad = []
    Xen = []
    Yen = []
    reac.add_left(-1)
    prod.add_left(-1)
    sdata.xllist.append(reac)
    sdata.yllist.append(prod)
    sdata.xslist.append(reac)
    sdata.yslist.append(prod)
    # DATABASE
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
        xsad.add_left(reac.get_iden())
        ysad.add_left(prod.get_iden())
        reac.add_right(xid)
        prod.add_right(yid)
        gmod = get_0mode()
        xloc = gopt(xsad, gmod)
        yloc = gopt(ysad, gmod)
        xid = update_iden(sdata.xllist, xloc)
        yid = update_iden(sdata.yllist, yloc)
        xloc.set_iden(xid)
        yloc.set_iden(yid)
        xloc.set_sm('M')
        yloc.set_sm('M')
        xloc.add_left(xsad.get_iden())
        yloc.add_left(ysad.get_iden())
        xsad.add_right(xid)
        ysad.add_right(yid)
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
        # ex = Xsad[ix].get_e() - reace
        ex = get_barrier(sdata.xllist, sdata.xslist, reac, Xloc[ix])
        for iy in range(itin.npop):
            fpy = Yloc[iy].get_lfp()
            # ey = Ysad[iy].get_e() - reace
            ey = get_barrier(sdata.yllist, sdata.yslist, prod, Yloc[iy])
            ee = max(ex, ey)
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            # mdist for multiobjective opt
            if dist < 1e-4:
                xdist = 1e-4
            else:
                xdist = dist
            mdist = np.log10(xdist) + ee
            print 'ZLOG: mdist, dist, log(dist), ee',\
                  mdist, dist, np.log(dist), ee
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
        # ex = Xsad[ix].get_e() - reace
        ex = get_barrier(sdata.xllist, sdata.xslist, reac, Xloc[ix])
        for iy in range(len(sdata.yllist)):
            fpy = sdata.yllist[iy].get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            # ey = sdata.yslist[iy].get_e() - reace
            ey = get_barrier(sdata.yllist, sdata.yslist, prod,
                             sdata.yllist[iy])
            ee = max(ex, ey)
            if dist < 1e-4:
                xdist = 1e-4
            else:
                xdist = dist
            mdist = np.log10(xdist) + ee
            print 'ZLOG: mdist, dist, log(dist), ee',\
                  mdist, dist, np.log(dist), ee
            xytdist.append(mdist)
        xytdistb = sorted(xytdist)[0]
        sdata.pdistx.append(xytdistb)

    # update pdist y
    for iy in range(itin.npop):
        yxtdist = []
        fpy = Yloc[iy].get_lfp()
        # ey = Ysad[iy].get_e() - reace
        ey = get_barrier(sdata.yllist, sdata.yslist, prod, Yloc[iy])
        for ix in range(len(sdata.xllist)):
            fpx = sdata.xllist[ix].get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            # ex = sdata.xslist[ix].get_e() - reace
            ex = get_barrier(sdata.xllist, sdata.xslist, reac,
                             sdata.xllist[ix])
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
            xid = update_iden(sdata.xslist, xsad)
            yid = update_iden(sdata.yslist, ysad)
            xsad.set_iden(xid)
            ysad.set_iden(yid)
            xloc.add_right(xid)
            yloc.add_right(yid)
            xsad.set_sm('S')
            ysad.set_sm('S')
            xsad.add_left(xloc.get_iden())
            ysad.add_left(yloc.get_iden())
            sdata.vx[ip] = cp(vx)
            sdata.vy[ip] = cp(vy)
            Xsad.append(xsad)
            Ysad.append(ysad)
            gmod = get_0mode()
            xxloc = gopt(xsad, gmod)
            yyloc = gopt(ysad, gmod)
            xid = update_iden(sdata.xllist, xxloc)
            yid = update_iden(sdata.yllist, yyloc)
            xxloc.set_iden(xid)
            yyloc.set_iden(yid)
            xxloc.set_sm('M')
            yyloc.set_sm('M')
            xxloc.add_left(xsad.get_iden())
            yyloc.add_left(ysad.get_iden())
            xsad.add_right(xid)
            ysad.add_right(yid)
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
            # ex = sdata.xslist[ix].get_e() - reace
            ex = get_barrier(sdata.xllist, sdata.xslist, reac,
                             sdata.xllist[ix])
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_lfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                # ey = sdata.yslist[iy].get_e() - reace
                ey = get_barrier(sdata.yllist, sdata.yslist, prod,
                                 sdata.yllist[iy])
                ee = max(ex, ey)
                if dist < 1e-4:
                    xdist = 1e-4
                else:
                    xdist = dist
                mdist = np.log10(xdist) + ee
                print 'ZLOG: mdist, dist, log(dist), ee',\
                      mdist, dist, np.log(dist), ee
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
            # ex = Xsad[ix].get_e() - reace
            ex = get_barrier(sdata.xllist, sdata.xslist, reac, Xloc[ix])
            xytdist = []
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_lfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                # ey = sdata.yslist[iy].get_e() - reace
                ey = get_barrier(sdata.yllist, sdata.yslist, prod,
                                 sdata.yllist[iy])
                ee = max(ex, ey)
                if dist < 1e-4:
                    xdist = 1e-4
                else:
                    xdist = dist
                mdist = np.log10(xdist) + ee
                print 'ZLOG: mdist, dist, log(dist), ee',\
                      mdist, dist, np.log(dist), ee
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
            # ey = Ysad[iy].get_e() - reace
            ey = get_barrier(sdata.yllist, sdata.yslist, prod, Yloc[iy])
            yxtdist = []
            for ix in range(len(sdata.xllist)):
                fpx = sdata.xllist[ix].get_lfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                # ex = sdata.xslist[ix].get_e() - reace
                ex = get_barrier(sdata.xllist, sdata.xslist, reac,
                                 sdata.xllist[ix])
                ee = max(ex, ey)
                if dist < 1e-4:
                    xdist = 1e-4
                else:
                    xdist = dist
                mdist = np.log10(xdist) + ee
                print 'ZLOG: mdist, dist, log(dist), ee',\
                      mdist, dist, np.log(dist), ee
                yxtdist.append([mdist, ix])
            yxtdistSort = sorted(yxtdist, key=lambda x: x[0])
            yxtbestdist = yxtdistSort[0][0]
            yxdist.append(yxtbestdist)
            if yxtbestdist < sdata.pdisty[iy]:
                sdata.pbesty[iy] = cp(Yloc[iy])
                sdata.pdisty[iy] = yxtbestdist

        sdata.xlocs.append(Xloc)
        sdata.ylocs.append(Yloc)

        # if abs(bestdist) < itin.dist:
        #     print "ZLOG: CONVERGED!"
        #     break


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


def connect_path(ine, mlisted, slisted, xm, xend, fatherids, xpath):
    # input: saddlelist, minimalist, npop, istep
    # reactant/product minimalist[0]
    # xend : the end point, either product or reactant
    # find the first neighbor saddle for xend

    # MERGE mlist and slist
    # mlisted = mergelr(mlist)
    # slisted = mergelr(mlist)

    # snode0 = []
    # snode0id = []
    # for xs in slist:
    #    if 0 in xs.get_nbor():
    #        snode0.append(xs)
    #        snode0id.append(xs.get_iden())

    # print 'm left', xm.get_left()

    for sp_id in xm.get_left():
        if sp_id not in fatherids and sp_id > -1:
            fatherids.append(sp_id)
            sp = getx_fromid(sp_id, slisted)
            e = sp.get_e() - ine
            xpath.create_node('Saddle'+str(sp.get_iden())+'E'+str(e),
                              sp.get_nid(), parent=xm.get_nid(), data=sp)
            for m_id in sp.get_left():
                if m_id == 0:
                    # connect the xend
                    xend.set_nid(xend.get_nid()-1)
                    xpath.create_node('END', xend.get_nid(),
                                      parent=sp.get_nid(), data=xend)
                    # print '# ZLOG: CONNECTED MID', m_id
                else:
                    # print '# ZLOG: SON ID', m_id
                    mp = getx_fromid(m_id, mlisted)
                    e = mp.get_e() - ine
                    xpath.create_node('Minima' + str(mp.get_iden())+'E'+str(e),
                                      mp.get_nid(), parent=sp.get_nid(), data=mp)
                    connect_path(ine, mlisted, slisted, mp, xend, fatherids, xpath)
    return 0


def getx_fromid(xid, listed):
    for xterm in listed:
        if xterm.get_iden() == xid:
            sdata.nidp += 1
            xterm.set_nid(sdata.nidp)
            return cp(xterm)
    print 'ERROR: getx_fromid', xid
    exit(1)


def mergelist(xlist):
    xlisted = []
    xid = []
    for xc in xlist:
        xid.append(xc.get_iden())

    print 'set(xid)', set(xid)
    for idt in set(xid):
        simit = []
        lt = []
        rt = []
        for xc in xlist:
            # print xc.get_iden()
            if xc.get_iden() == idt:
                simit.append(xc)
                lt += xc.get_left()
                rt += xc.get_right()
                xt = cp(xc)
        ltt = list(set(lt))
        rtt = list(set(rt))
        xt.set_left(ltt)
        xt.set_right(rtt)
        xlisted.append(xt)
    return xlisted


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


def get_barrier(mlist, slist, startp, endp):
    # get barrier energy for startp (reactant/product)
    # to endp (one local minima)
    mlisted = mergelist(mlist)
    slisted = mergelist(slist)
    endp.set_nid(-1)
    ine = endp.get_e()
    if endp.get_iden() > 0:
        xpath = Tree()
        fatherids = []
        endp.set_nid(0)
        xpath.create_node('M' + str(endp.get_iden()), 0, data=endp)
        sdata.nidp = 0
        connect_path(ine, mlisted, slisted, endp, startp, fatherids, xpath)

        dd = []
        for xx in xpath.all_nodes():
            if xx.identifier < 0:
                d = []
                d.append(xx.data)
                nid = xx.bpointer
                nxx = xpath.get_node(nid)
                while True:
                    d.append(nxx.data)
                    if nxx.is_root():
                        break
                    nid = nxx.bpointer
                    nxx = xpath.get_node(nid)
                dd.append(d)
        xe = []
        for x in dd:
            ee = []
            for dx in x:
                exx = dx.get_e() - sdata.reace
                ee.append(exx)
            xe.append(max(ee))
        mxe = min(xe)
    else:
        mxe = 0.0

    return mxe




def main():
    (reac, prod) = initrun()
    woo(reac, prod)
    f = open('xm.dat', 'w')
    pick.dump(sdata.xllist, f)
    f.close()
    f = open('xs.dat', 'w')
    pick.dump(sdata.xslist, f)
    f.close()
    f = open('ym.dat', 'w')
    pick.dump(sdata.yllist, f)
    f.close()
    f = open('ys.dat', 'w')
    pick.dump(sdata.yslist, f)
    f.close()
    # outputw()


def utest1():
    f = open('xm.dat')
    xmlist = pick.load(f)
    f.close()
    f = open('xs.dat')
    xslist = pick.load(f)
    f.close()
    xmlisted = mergelist(xmlist)
    xslisted = mergelist(xslist)
    print 'nxmlist, nxmlisted', len(xmlist), len(xmlisted)
    print 'nxslist, nsmlisted', len(xslist), len(xslisted)
    xend = cp(xmlist[0])
    xm = cp(xmlist[-3])
    print xm.get_iden()
    fatherids = []
    xpath = Tree()
    xm.set_nid(0)
    xpath.create_node('Minima' + str(xm.get_iden()), 0, data=xm)
    xpath.show()
    sdata.nidp = 0
    xend.set_nid(-1)
    connect_path(xmlisted, xslisted, xm, xend, fatherids, xpath)
    xpath.show()

    dd = []
    for xxx in xpath.all_nodes():
        if xxx.identifier < 0:
            print 'WA'
            d = []
            d.append(xxx.data)
            nid = xxx.bpointer
            xx = xpath.get_node(nid)
            while True:
                d.append(xx.data)
                if xx.is_root():
                    break
                nid = xx.bpointer
                xx = xpath.get_node(nid)
            dd.append(d)

    print len(dd)

    for x in dd:
        print 'EE',
        ee = []
        for xx in x:
            print xx.get_e() - xend.get_e(),
            ee.append(xx.get_e() - xend.get_e())
        print
        print max(ee)


def utest2():
    f = open('xm.dat')
    xmlist = pick.load(f)
    f.close()
    f = open('xs.dat')
    xslist = pick.load(f)
    f.close()
    f = open('ym.dat')
    ymlist = pick.load(f)
    f.close()
    f = open('ys.dat')
    yslist = pick.load(f)
    f.close()

    xmlisted = mergelist(xmlist)
    xslisted = mergelist(xslist)
    ymlisted = mergelist(ymlist)
    yslisted = mergelist(yslist)

    xend = cp(xmlist[0])
    # print 'xend', xend.get_e()
    yend = cp(ymlist[0])

    sdata.types = xend.get_types()

    goodlist = []
    for xx in xmlist:
        fpx = xx.get_lfp()
        for yy in ymlist:
            fpy = yy.get_lfp()
            (d, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            if d < itin.dist:
                print 'dd', d
                goodlist.append([xx, yy])
                print xx.get_iden(), yy.get_iden()

    print 'len', len(goodlist)

    ine = xend.get_e()

    xend.set_nid(-1)
    yend.set_nid(-1)
    kk = 0
    mxy = []
    for xyxy in goodlist:
        print 'kk', kk
        kk += 1
        xx = xyxy[0]
        yy = xyxy[1]
        print 'xx, yy', xx.get_iden(), yy.get_iden()
        if xx.get_iden() > 0:
            xpath = Tree()
            fatherids = []
            xx.set_nid(0)
            xpath.create_node('Minima'+str(xx.get_iden()), 0, data=xx)
            sdata.nidp = 0
            connect_path(ine, xmlisted, xslisted, xx, xend, fatherids, xpath)

        if yy.get_iden() > 0:
            ypath = Tree()
            fatherids = []
            yy.set_nid(0)
            ypath.create_node('Minima'+str(yy.get_iden()), 0, data=yy)
            sdata.nidp = 0
            connect_path(ine, ymlisted, yslisted, yy, yend, fatherids, ypath)

        if xx.get_iden() > 0:
            dd = []
            for xxx in xpath.all_nodes():
                if xxx.identifier < 0:
                    d = []
                    d.append(xxx.data)
                    nid = xxx.bpointer
                    nxx = xpath.get_node(nid)
                    while True:
                        d.append(nxx.data)
                        if nxx.is_root():
                            break
                        nid = nxx.bpointer
                        nxx = xpath.get_node(nid)
                    dd.append(d)
            xe = []
            for x in dd:
                ee = []
                for dx in x:
                    exx = dx.get_e() - xend.get_e()
                    ee.append(exx)
                xe.append(max(ee))
            mxe = min(xe)
        else:
            mxe = 0.0

        if yy.get_iden() > 0:
            dd = []
            for yyy in ypath.all_nodes():
                if yyy.identifier < 0:
                    d = []
                    d.append(yyy.data)
                    nid = yyy.bpointer
                    nyy = ypath.get_node(nid)
                    while True:
                        d.append(nyy.data)
                        if nyy.is_root():
                            break
                        nid = nyy.bpointer
                        nyy = ypath.get_node(nid)
                    dd.append(d)
            ye = []
            for y in dd:
                ee = []
                for dy in y:
                    exx = dy.get_e() - xend.get_e()
                    ee.append(exx)
                ye.append(max(ee))
            mye = min(ye)
        else:
            mye = 0.0

        print 'BAR', max(mxe, mye), mxe, mye
        mxy.append([max(mxe, mye), xx, yy])

    mxysort = sorted(mxy, key=lambda x: x[0])
    print mxysort[0][0]
    print mxysort[0][1].get_iden()
    print mxysort[0][2].get_iden()

    goodx = cp(mxysort[0][1])
    idenx = goodx.get_iden()
    fatherids = []
    xpath = Tree()
    goodx.set_nid(0)
    xpath.create_node('Minima' + str(goodx.get_iden()), 0, data=goodx)
    xpath.show()
    sdata.nidp = 0
    xend.set_nid(-1)
    connect_path(ine, xmlisted, xslisted, goodx, xend, fatherids, xpath)
    xpath.show()
    xpath.save2file('xp.dat')






if __name__ == "__main__":
    # main()
    utest2()

















