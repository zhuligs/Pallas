#!/usr/bin/python -u

import itin
from zfunc import set_cell_from_vasp
import cPickle as pick
import os


def ovjob():
    try:
        pcell = set_cell_from_vasp('CONTCAR')
        try:
            e = float(os.popen("awk '/free  energy/{print $5}' OUTCAR|tail -1").read())
            h = itin.press * pcell.get_volume() / 1602.2 + e
        except:
            h = 31118.
        pcell.set_e(h)
    except:
        pcell = set_cell_from_vasp('POSCAR')
        h = 31118.
        pcell.set_e(h)
    f = open('pcell.bin', 'w')
    pick.dump(pcell, f)
    f.close()

    # print h


if __name__ == '__main__':
    ovjob()
