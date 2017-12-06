#!/usr/bin/python -u
####
# import numpy
import networkx as nx
# import matplotlib.pyplot as plt
import cPickle as pick
import itin
import fppy
from copy import deepcopy as cp
import sdata

from fplib2 import fpdist


# from w20 import get_barrier


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

    types = xl[0].get_types()
    print types
    return

    dd = []
    for x in xl:
        fpx = x.get_lfp()
        for y in yl:
            fpy = y.get_lfp()
            (d1, m) = fppy.fp_dist(itin.ntyp, types, fpx, fpy)
            d2 = fpdist(itin.ntyp, types, fpx, fpy)
            dd.append([d1, d2])

    f = open('t28.out', 'w')
    for x in dd:
        f.write("%g  %g\n" % tuple(x))
    f.close()
    return 0


    # for i in range(len(xl) - 1):
    #     fp1 = xl[i].get_lfp()
    #     for j in range(i + 1, len(xl)):
    #         fp2 = xl[j].get_lfp()
    #         (d, m) = fppy.fp_dist(itin.ntyp, types, fp1, fp2)
    #         if d < 0.001:
    #             print xl[i].get_iden(), xl[j].get_iden()

    # return 0

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
        if e > 2:
            e = 2.
        v = float(xll.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v, Polygon=0)
        # G.add_node(xll, energy=e, volume=v, Polygon=0)
        # G.add_node(nm, energy=e)
        ees.append(e)
    for xss in xsed:
        iden = xss.get_iden()
        nm = 'xs' + str(iden)
        e = xss.get_e() - e0
        if e > 2:
            e = 2.
        v = float(xss.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v, Polygon=4)
        # G.add_node(xss, energy=e, volume=v, Polygon=4)
        # G.add_node(nm, energy=e)
        ees.append(e)
    for yll in yled:
        iden = yll.get_iden()
        nm = 'yl' + str(iden)
        e = yll.get_e() - e0
        if e > 2:
            e = 2.
        v = float(yll.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v, Polygon=0)
        # G.add_node(yll, energy=e, volume=v, Polygon=0)
        # G.add_node(nm, energy=e)
        ees.append(e)
    for yss in ysed:
        iden = yss.get_iden()
        nm = 'ys' + str(iden)
        e = yss.get_e() - e0
        if e > 2:
            e = 2.
        v = float(yss.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v, Polygon=4)
        # G.add_node(yss, energy=e, volume=v, Polygon=4)
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
                if d < 0.0001:
                    d = 0.0001
                d = 1. / d
                G.add_edge(nmx, nmy, weight=d)
                # G.add_edge(nmx, nmy)

    for xll in xled:
        fpx = xll.get_lfp()
        ix = xll.get_iden()
        nmx = 'xl' + str(ix)
        lr = xll.get_left() + xll.get_right()
        for iid in lr:
            if iid > -1:
                nsx = 'xs' + str(iid)
                xsl = getx_fromid(iid, xsed)
                fpxx = xsl.get_lfp()
                (d, m) = fppy.fp_dist(itin.ntyp, types, fpx, fpxx)
                if d < 0.0001:
                    d = 0.0001
                d = 1. / d
                G.add_edge(nmx, nsx, weight=d)
                # G.add_edge(nmx, nsx)

    for yll in yled:
        fpy = yll.get_lfp()
        iy = yll.get_iden()
        nmy = 'yl' + str(iy)
        lr = yll.get_left() + yll.get_right()
        for iid in lr:
            if iid > -1:
                nsy = 'ys' + str(iid)
                ysl = getx_fromid(iid, ysed)
                fpyy = ysl.get_lfp()
                (d, m) = fppy.fp_dist(itin.ntyp, types, fpy, fpyy)
                if d < 0.0001:
                    d = 0.0001
                d = 1. / d
                G.add_edge(nmy, nsy, weight=d)
                # G.add_edge(nmy, nsy)

    # nx.write_gexf(G, "test.gexf")
    # nx.draw(G)
    # plt.show()

    for xss in xsed:
        fpx = xss.get_lfp()
        ix = xss.get_iden()
        nmx = 'xs' + str(ix)
        ll = xss.get_left() + xss.get_right()
        for iid in ll:
            if iid > -1:
                nmm = 'xl' + str(iid)
                xll = getx_fromid(iid, xled)
                fppx = xll.get_lfp()
                (d, m) = fppy.fp_dist(itin.ntyp, types, fpx, fppx)
                if d < 0.0001:
                    d = 0.0001
                d = 1. / d
                G.add_edge(nmx, nmm, weight=d)
                #G.add_edge(nmx, nmm)

    # for yy in yled:
        # print yy.get_iden()

    for yss in ysed:
        fpy = yss.get_lfp()
        iy = yss.get_iden()
        nmy = 'ys' + str(iy)
        ll = yss.get_left() + yss.get_right()
        for iid in ll:
            if iid > -1:
                nmm = 'yl' + str(iid)
                try:
                    yll = getx_fromid(iid, yled)
                except:
                    continue
                # fppy =
                fpyy = yll.get_lfp()
                (d, m) = fppy.fp_dist(itin.ntyp, types, fpy, fpyy)
                if d < 0.0001:
                    d = 0.0001
                d = 1. / d
                G.add_edge(nmy, nmm, weight=d)

    # G.add_edge('xl0', 'xl10', color='red')
    # G.add_edge('xl10', 'yl0', color='red')

    # G.add_edge('xl0', 'xl186')
    #sdata.reace = xl[0].get_e()
    # print sdata.reace
    # glist = []
    # for xxl in xl:
    #     fpx = xxl.get_lfp()
    #     # eex = get_barrier(xl, xs, xl[0], xxl)
    #     for yyl in yl:
    #         fpy = yyl.get_lfp()
    #         (d, m) = fppy.fp_dist(itin.ntyp, types, fpx, fpy)
    #         if d < itin.dist:
    #             glist.append((xxl, yyl))
    #             print xxl.get_iden(), yyl.get_iden(), xxl.get_e(), yyl.get_e()
    #
    # xx = getx_fromid(3, xl)
    # print xx.get_left()
    # print xx.get_right()
    # xx1 = getx_fromid(11, xs)
    # print xx1.get_left()
    # print xx1.get_right()

    #G.add_edge('xl0', 'xs11', color='red')
    #G.add_edge('xs11', 'xl3', color='red')
    #G.add_edge('xl3', 'ys2', color='red')
    #G.add_edge('ys2', 'yl0', color='red')
    # nx.write_gexf(G, "test.gexf")
    # nx.write_gml(G, 'network.gml')

    # for path in nx.all_simple_paths(G, source='xl0', target='yl0'):
    #    print path
    paths = nx.all_simple_paths(G, source='xl0', target='yl0', cutoff=10)
    # print len(list(paths))
    for path in paths:
        print path
        ee = []
        for node in path:
            e = G.node[node]['energy']
            ee.append(e)
        be = max(ee)
        print 'energy', be, ee




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


# def mergl(xlist):
#     n = len(xlist)
#     types = xlist[0].get_types()
#     for i in range(n - 1):
#         xi = xlist[i]
#         x = cp(xi)
#         fpi = xi.get_lfp()
#         xleft = xi.get_left()
#         xright = xi.get_right()
#         for j in range(i + 1, n):
#             xj = xlist[j]
#             fpj = xj.get_lfp()
#             (d, m) = fppy.fp_dist(itin.ntyp, types, fpi, fpj)
#             if d < itin.dist:
#                 left = fpj.get_left()       
#                 right = fpj.get_right()
#                 xleft += left
#                 xright += right
#         x.set_left(xleft)
#         x.set_righ()


def getx_fromid(xid, listed):
    for xterm in listed:
        if xterm.get_iden() == xid:
            sdata.nidp += 1
            xterm.set_nid(sdata.nidp)
            return cp(xterm)
    print 'ERROR: getx_fromid', xid
    return xterm


if __name__ == '__main__':
    netx()
