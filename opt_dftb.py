#!/usr/bin/python -u 

import itin
from zfunc import write_cell_to_vasp, set_cell_from_vasp, getx


import numpy as np
from copy import deepcopy as cp
import cPickle as pick

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
# from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.dftb import Dftb
from tsase.dimer import ssdimer
from tsase.dimer import lanczos
from tsase.neb.util import vunit, vrand

def gopt_dftb(amode):
    p1 = read('POSCAR', format='vasp')
    # tags = [a.symbol == 'Si' for a in p1]

    #parameters = {'mass': ['1 1.0'], 'pair_style': 'lj/sf 2.5',
    #              'pair_coeff': ['1 1  1.0  1.0 2.5'],
    #              'pair_modify': 'shift yes'}
    # parameters = itin.parameters
    # calc = LAMMPS(parameters=parameters)
    calc = Dftb(label='SiO2',
            atoms=p1,
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_Si='"p"',
            Hamiltonian_MaxAngularMomentum_O='"p"',
            kpts=(5,5,5),)

    p1.set_calculator(calc)

    natom = len(p1)
    vol = p1.get_volume()
    jacob = (vol/natom)**(1.0/3.0) * natom**0.5
    mode = np.zeros((len(p1)+3, 3))
    mode = vrand(mode)
    # try:
    #     mode = vunit(mode)
    # except:
    #     pass
    if amode:
        cellt = p1.get_cell() + np.dot(p1.get_cell(), mode[-3:]/jacob)
        p1.set_cell(cellt, scale_atoms=True)
        p1.set_positions(p1.get_positions() + mode[:-3])

    pst = itin.press / 10.

    pstress = np.zeros((3, 3))
    pstress[0][0] = pst
    pstress[1][1] = pst
    pstress[2][2] = pst

    p1box = mushybox(p1, pstress)
    # print len(p1box)
    # print p1box.jacobian
    # print p1box.get_potential_energy()
    try:
        dyn = FIRE(p1box, dt=0.1, maxmove=0.2, dtmax=0.2)
        dyn.run(fmax=0.001, steps=2000)
        write("CONTCAR", p1, format='vasp', direct=True)
        pcell = set_cell_from_vasp("CONTCAR")
        e = p1box.get_potential_energy()
        h = itin.press * pcell.get_volume() / 1602.176487 + e
        pcell.set_e(h)
    except:
        pcell = set_cell_from_vasp("POSCAR")
        pcell.set_e(151206)

    f = open('pcell.bin', 'w')
    pick.dump(pcell, f)
    f.close()

    return pcell


if __name__ == '__main__':
    amode = True
    gopt_dftb(amode)


