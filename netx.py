#!/usr/bin/python -u

import numpy
import networkx as nx
import matplotlib.pyplot as plt
import cPickle as pick
import itin
import fppy

from w20 import mergelist


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
        G.add_node(nm, energy=e, volume=v)
        # G.add_node(nm, energy=e)
        ees.append(e)
    for xss in xsed:
        iden = xss.get_iden()
        nm = 'xs' + str(iden)
        e = xss.get_e() - e0
        if e > 2: e = 2.
        v = float(xss.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v)
        # G.add_node(nm, energy=e)
        ees.append(e)
    for yll in yled:
        iden = yll.get_iden()
        nm = 'yl' + str(iden)
        e = yll.get_e() - e0
        if e > 2: e = 2.
        v = float(yll.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v)
        # G.add_node(nm, energy=e)
        ees.append(e)
    for yss in ysed:
        iden = yss.get_iden()
        nm = 'ys' + str(iden)
        e = yss.get_e() - e0
        if e > 2: e = 2.
        v = float(yss.get_volume() / itin.nat)
        vvs.append(v)
        G.add_node(nm, energy=e, volume=v)
        # G.add_node(nm, energy=e)
        ees.append(e)

    print sorted(ees)
    print sorted(vvs)

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
                G.add_edge(nmx, nmy)

    nx.write_gexf(G, "test.gexf")
    nx.draw(G)
    # plt.show()





if __name__ == '__main__':
    netx()
