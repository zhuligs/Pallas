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
from ase.calculators.lammpsrun import LAMMPS
from tsase.dimer import ssdimer
from tsase.dimer import lanczos
from tsase.neb.util import vunit, vrand

def gopt_lammps():
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
    # jacob = (vol/natom)**(1.0/3.0) * natom**0.5
    # mode = np.zeros((len(p)+3, 3))
    # mode = vrand(mode)
    # try:
    #     mode = vunit(mode)
    # except:
    #     pass
    # cellt = p1.get_cell() + np.dot(p1.get_cell(), mode[-3:]/jacob)
    # p1.set_cell(cellt, scale_atoms=True)
    # p1.set_positions(p1.get_positions() + mode[:-3])

    pstress = p1.get_cell()*0.0

    p1box = mushybox(p1, pstress)
    # print len(p1box)
    # print p1box.jacobian
    # print p1box.get_potential_energy()
    try:
        dyn = FIRE(p1box, dt=0.1, maxmove=0.2, dtmax=0.2)
        dyn.run(fmax=0.001, steps=2000)

        io.write("CONTCAR", p1, format='vasp')
        pcell = set_cell_from_vasp("CONTCAR")
        e = p1box.get_potential_energy()
        pcell.set_e(e)
    except:
        pcell = set_cell_from_vasp("POSCAR")
        pcell.set_e(151206)

    f = open('pcell.bin', 'w')
    pick.dump(pcell, f)
    f.close()

    return pcell

if __name__ == '__main__':
    gopt_lammps()


