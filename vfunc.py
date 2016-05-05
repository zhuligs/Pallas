#!/usr/bin/env python

import itin
import numpy as np
from zfunc import write_cell_to_vasp, set_cell_from_vasp
import os


def goptv(xcell, mode):
    lat = xcell.get_lattice()
    vol = xcell.get_volume()
    jacob = (vol / itin.nat)**(1.0/3.0) * itin.nat**0.5
    latt = lat + np.dot(lat, mode[-3:]/jacob)
    xcell.set_lattice(latt)
    newpos = xcell.get_positions() + mode[:-3]
    xcell.set_positions(newpos)
    write_cell_to_vasp(xcell, "POSCAR")
    os.system("cp INCAR_OPT INCAR")
    os.system("sh runvasp.sh")
    e = float(os.popen("awk '/free  energy/{print $5}' OUTCAR|tail -1").read())
    pcell = set_cell_from_vasp("CONTCAR")
    pcell.set_e(e)

    return pcell


def runvdim(xcell, mode):
    lat = xcell.get_lattice()
    vol = xcell.get_volume()
    jacob = (vol/itin.nat)**(1.0/3.0) * itin.nat**0.5
    latt = lat + np.dot(lat, mode[-3:]/jacob)
    xcell.set_lattice(latt)
    f = open('MODECAR', 'w')
    for x in mode[:-3]:
        f.write("%15.9f  %15.9f  %15.9f\n" % tuple(x))
    f.close()
    write_cell_to_vasp(xcell, "POSCAR")
    os.system("cp INCAR_DIM INCAR")
    os.system("sh runvasp.sh")
    e = float(os.popen("awk '/free  energy/{print $5}' OUTCAR|tail -1").read())
    pcell = set_cell_from_vasp("CONTCAR")
    pcell.set_e(e)
    return pcell


