import numpy as np
from itdbase import Cell
import itin
from copy import deepcopy as cp

# from tsase.optimize import MDMin
from ase.optimize.fire import FIRE
# from ase.optimize import BFGS
from ase import *
from ase.io import read, write
import os
import sys
import numpy as np
from tsase.mushybox import mushybox
# from tsase.calculators.vasp_ext import Vasp
# from tsase.calculators.lammps_ext import LAMMPS
from ase.calculators.lammpsrun import LAMMPS
from tsase.dimer import ssdimer
from tsase.dimer import lanczos
from tsase.neb.util import vunit, vrand


def gopt(xcell, mode):
    write_cell_to_vasp(xcell, 'POSCAR')
    p1 = read('POSCAR', format='vasp')
    # tags = [a.symbol == 'Si' for a in p1]

    #parameters = {'mass': ['1 1.0'], 'pair_style': 'lj/sf 2.5',
    #              'pair_coeff': ['1 1  1.0  1.0 2.5'],
    #              'pair_modify': 'shift yes'}
    parameters = itin.parameters
    calc = LAMMPS(parameters=parameters)

    p1.set_calculator(calc)

    natom = len(p1)
    vol = p1.get_volume()
    jacob = (vol/natom)**(1.0/3.0) * natom**0.5
    # mode = np.zeros((len(p)+3, 3))
    # mode = vrand(mode)
    # try:
    #     mode = vunit(mode)
    # except:
    #     pass
    cellt = p1.get_cell() + np.dot(p1.get_cell(), mode[-3:]/jacob)
    p1.set_cell(cellt, scale_atoms=True)
    p1.set_positions(p1.get_positions() + mode[:-3])

    pstress = p1.get_cell()*0.0

    p1box = mushybox(p1, pstress)
    # print len(p1box)
    # print p1box.jacobian
    # print p1box.get_potential_energy()
    try:
        dyn = FIRE(p1box, dt=0.1, maxmove=0.2, dtmax=0.2)
        dyn.run(fmax=0.01, steps=2000)

        io.write("CONTCAR", p1, format='vasp')
        pcell = set_cell_from_vasp("CONTCAR")
        e = p1box.get_potential_energy()
        pcell.set_e(e)
    except:
        pcell = cp(xcell)
        pcell.set_e(151206)

    return pcell


def rundim(xcell, mode):
    write_cell_to_vasp(xcell, 'DCAR')
    p = read('DCAR', format='vasp')
    parameters = itin.parameters
    calc = LAMMPS(parameters=parameters)
    p.set_calculator(calc)
    # E0 = p.get_potential_energy()
    natom = len(p)
    vol = p.get_volume()
    jacob = (vol/natom)**(1.0/3.0) * natom**0.5
    # mode = np.zeros((len(p)+3, 3))
    # mode = vrand(mode)
    # try:
    #     mode = vunit(mode)
    # except:
    #     pass
    cellt = p.get_cell() + np.dot(p.get_cell(), mode[-3:]/jacob)
    p.set_cell(cellt, scale_atoms=True)
    p.set_positions(p.get_positions() + mode[:-3])
    d = lanczos.lanczos_atoms(p, mode=mode, rotationMax=4,
                              ss=True, phi_tol=15)
    dyn = FIRE(d, dt=0.1, maxmove=0.2, dtmax=0.2)
    try:
        dyn.run(fmax=0.05, steps=2000)
        E1 = p.get_potential_energy()
        write("CDCAR", d.R0, format='vasp', direct=True)
        pcell = set_cell_from_vasp("CDCAR")
        pcell.set_e(E1)
    except:
        pcell = cp(xcell)
        pcell.set_e(151206)
    return pcell


def set_cell_from_vasp(pcar):
    xcell = Cell()
    buff = []
    with open(pcar) as f:
        for line in f:
            buff.append(line.split())

    lat = np.array(buff[2:5], float)
    typt = np.array(buff[5], int)
    pos = np.array(buff[7:], float)
    xcell.set_name(itin.sname)
    xcell.set_lattice(lat)
    if buff[6][0].strip()[0] == 'D':
        xcell.set_positions(pos)
    else:
        xcell.set_cart_positions(pos)
    xcell.set_typt(typt)
    xcell.set_znucl(itin.znucl)
    xcell.set_types()
    xcell.cal_fp(itin.fpcut, itin.lmax)
    return xcell


def write_cell_to_vasp(xcell, pcar):
    lat = xcell.get_lattice()
    typt = xcell.get_typt()
    pos = xcell.get_positions()

    f = open(pcar, 'w')
    f.write(itin.sname + '\n')
    f.write('1.0\n')
    for x in lat:
        f.write("%15.9f %15.9f %15.9f\n" % tuple(x))
    for ix in typt:
        f.write(str(ix) + '  ')
    f.write('\n')
    f.write('Direct\n')
    for x in pos:
        f.write("%15.9f %15.9f %15.9f\n" % tuple(x))
    f.close()


def getx(cell1, cell2):
    mode = np.zeros((itin.nat + 3, 3))
    mode[-3:] = cell1.get_lattice() - cell2.get_lattice()
    ilat = np.linalg.inv(cell1.get_lattice())
    vol = cell1.get_volume()
    jacob = (vol / itin.nat)**(1.0 / 3.0) * itin.nat**0.5
    mode[-3:] = np.dot(ilat, mode[-3:]) * jacob
    pos1 = cell1.get_cart_positions()
    pos2 = cell2.get_cart_positions()
    for i in range(itin.nat):
        mode[i] = pos1[i] - pos2[i]
    mode = vunit(mode)
    return mode
