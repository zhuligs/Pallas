#!/usr/bin/python -u

from pele.storage import Database
from pele.landscape import database2graph
from pele.utils.disconnectivity_graph import DisconnectivityGraph
import matplotlib.pyplot as plt
import sdata
from copy import deepcopy as cp
import cPickle as pick


def discon():
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

    db = Database(db="pathway.sqlite", accuracy=0.001)

    # allmin = []
    # for x in xl:
    #     ex = x.get_e()
    #     coorx = x.get_positions()
    #     minx = db.addMinimum(ex, coorx)
    #     allmin.append(minx)
    for x in xs:
        # esadx = x.get_e()
        # cosx = x.get_positions()
        left = x.get_left()
        right = x.get_right()
        xleft = getx_fromid(left[0], xl)
        xright = getx_fromid(right[0], xl)
        min1 = db.addMinimum(xleft.get_e(), xleft.get_positions())
        min2 = db.addMinimum(xright.get_e(), xright.get_positions())
        db.addTransitionState(x.get_e(), x.get_positions(), min1, min2)

    for x in ys:
        # esadx = x.get_e()
        # cosx = x.get_positions()
        left = x.get_left()
        right = x.get_right()
        xleft = getx_fromid(left[0], yl)
        xright = getx_fromid(right[0], yl)
        min1 = db.addMinimum(xleft.get_e(), xleft.get_positions())
        min2 = db.addMinimum(xright.get_e(), xright.get_positions())
        db.addTransitionState(x.get_e(), x.get_positions(), min1, min2)

    graph = database2graph(db)
    dg = DisconnectivityGraph(graph)
    dg.calculate()
    dg.plot()
    plt.show()
    plt.savefig("tree.pdf")


def getx_fromid(xid, listed):
    for xterm in listed:
        if xterm.get_iden() == xid:
            sdata.nidp += 1
            xterm.set_nid(sdata.nidp)
            return cp(xterm)
    print 'ERROR: getx_fromid', xid
    return xterm


if __name__ == '__main__':
    discon()
