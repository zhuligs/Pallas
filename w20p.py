#!/usr/bin/python -u

import cPickle as pick
import fppy
import sdata
import itin


def w20p():
    f = open('xllist.bin')
    xllist = pick.load(f)
    f.close()
    f = open('yllist.bin')
    yllist = pick.load(f)
    f.close()

    sdata.types = xllist[0].get_types()

    datt = []
    for i in range(len(xllist)):
        fpx = xllist[i].get_lfp()
        for j in range(len(yllist)):
            fpy = yllist[j].get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            datt.append([dist, [i, j]])

    sdatt = sorted(datt, key=lambda x: x[0])
    print sdatt[0]
    print len(xllist)


if __name__ == '__main__':
    w20p()
