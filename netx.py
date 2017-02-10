#!/usr/bin/python -u

import numpy
import networkx as nx
import matplotlib.pyplot as plt
import cPickle as pick
import itin
import fppy
from copy import deepcopy as cp


def foo():
    G = nx.Graph()
    # G.add_nodes_from([0, 1, 2, 3, 4, 5])
    G.add_node(0, energy=1.1)
    G.add_node(1, energy=1.3)
    G.add_node(2, energy=2.0)
    G.add_node(3, energy=3.0)
    G.add_node(4, energy=1.5)
    G.add_node(5, energy=2.1)
    G.add_edge(0, 1)
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    G.add_edge(3, 4)
    G.add_edge(4, 5)
    G.add_edge(5, 1)
    nx.write_gexf(G, "test1.gexf")
    # nx.draw(G)
    # plt.show()
    print G.nodes()

    # G.nodes(data=True)
    print G.nodes(data=True)


def netx():
    f = open('xllist.bin')
    xl = pick.load(f)
    f.close()
    f = open('xslist.bin')
    xs = pick.load(f)
    f.close()
    f = open('yllist.bin')
    yl = pick.load(f)
    f.close()
    f = open('yslist.bin')
    ys = pick.load(f)
    f.close()

    xled = mergelist(xl)
    xsed = mergelist(xs)
    yled = mergelist(yl)
    ysed = mergelist(ys)

    e0 = xl[0].get_e()
    ees = []
    vvs = []
    G = nx.Graph()
    for xll in xled:
        iden = xll.get_iden()
        nm = 'xl' + str(iden)
        e = xll.get_e() - e0
        if e > 2: e = 2.
        v = float(xll.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v, Polygon=0)
        # G.add_node(nm, energy=e)
        ees.append(e)
    for xss in xsed:
        iden = xss.get_iden()
        nm = 'xs' + str(iden)
        e = xss.get_e() - e0
        if e > 2: e = 2.
        v = float(xss.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v, Polygon=4)
        # G.add_node(nm, energy=e)
        ees.append(e)
    for yll in yled:
        iden = yll.get_iden()
        nm = 'yl' + str(iden)
        e = yll.get_e() - e0
        if e > 2: e = 2.
        v = float(yll.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v, Polygon=0)
        # G.add_node(nm, energy=e)
        ees.append(e)
    for yss in ysed:
        iden = yss.get_iden()
        nm = 'ys' + str(iden)
        e = yss.get_e() - e0
        if e > 2: e = 2.
        v = float(yss.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v, Polygon=4)
        # G.add_node(nm, energy=e)
        ees.append(e)

    # print sorted(ees)
    # print sorted(vvs)

    types = xled[0].get_types()
    for xll in xled:
        fpx = xll.get_lfp()
        ix = xll.get_iden()
        nmx = 'xl' + str(ix)
        for yll in yled:
            fpy = yll.get_lfp()
            iy = yll.get_iden()
            nmy = 'yl' + str(iy)
            (d, m) = fppy.fp_dist(itin.ntyp, types, fpx, fpy)
            if d < itin.dist:
                G.add_edge(nmx, nmy, weight=d)

    for xll in xled:
        fpx = xll.get_lfp()
        ix = xll.get_iden()
        nmx = 'xl' + str(ix)
        lr = xll.get_left() + xll.get_right()
        for iid in lr:
            nsx = 'xs' + str(iid)
            xsl = getx_fromid(iid, xsed)
            fpxx = xsl.get_lfp()
            (d, m) = fppy.fp_dist(itin.ntyp, types, fpx, fpxx)
            G.add_edge(nmx, nsx, weight=d)

    for yll in yled:
        fpy = yll.get_lfp()
        iy = yll.get_iden()
        nmy = 'yl' + str(iy)
        lr = yll.get_left() + yll.get_right()
        for iid in lr:
            nsy = 'ys' + str(iid)
            ysl = getx_fromid(iid, ysed)
            fpyy = ysl.get_lfp()
            (d, m) = fppy.fp_dist(itin,ntyp, types, fpy, fpyy)
            G.add_edge(nmy, nsy, weight=d)


    nx.write_gexf(G, "test.gexf")
    nx.draw(G)
    # plt.show()


def mergelist(xlist):
    # xlist = []
    # for xx in xlist0:
    #     if xx.get_e() < 1000:
    #         xlist.append(xx)
    xlisted = []
    xid = []
    for xc in xlist:
        xid.append(xc.get_iden())

    # print 'set(xid)', set(xid)
    for idt in set(xid):
        simit = []
        lt = []
        rt = []
        for xc in xlist:
            # print xc.get_iden()
            if xc.get_iden() == idt:
                simit.append(xc)
                # print 'lt'
                # print lt
                # print xc.get_left()
                # print 'rt'
                # print rt
                # print xc.get_right()
                lt += xc.get_left()
                rt += xc.get_right()
                xt = cp(xc)
        ltt = list(set(lt))
        rtt = list(set(rt))
        xt.set_left(ltt)
        xt.set_right(rtt)
        xlisted.append(xt)
    return xlisted


def getx_fromid(xid, listed):
    for xterm in listed:
        if xterm.get_iden() == xid:
            sdata.nidp += 1
            xterm.set_nid(sdata.nidp)
            return cp(xterm)
    print 'ERROR: getx_fromid', xid
    exit(1)


if __name__ == '__main__':
    netx()
