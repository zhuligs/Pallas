#!/usr/bin/env python

import itin
import numpy as np
from zfunc import write_cell_to_vasp, set_cell_from_vasp
import os
import glob
import sdata


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
    h = itin.press * pcell.get_volume() / 1602.2 + e
    pcell.set_e(h)
    gdirs = glob.glob('Gdir*')
    gdir = 'Gdir' + str(len(gdirs))
    os.system('mkdir -p ' + gdir)
    os.system('cp POSCAR OUTCAR CONTCAR XDATCAR ' + gdir)
    sdata.gdir = gdir

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
    h = itin.press * pcell.get_volume() / 1602.2 + e
    pcell.set_e(h)
    ddirs = glob.glob('Ddir*')
    ddir = 'Ddir' + str(len(ddirs))
    os.system('mkdir -p ' + ddir)
    os.system('cp POSCAR MODECAR OUTCAR XDATCAR DIMCAR ' + ddir)
    sdata.ddir = ddir

    return pcell


