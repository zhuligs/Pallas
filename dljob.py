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


def rundim_lammps():
    p = read('PRESAD.vasp', format='vasp')
    parameters = itin.parameters
    calc = LAMMPS(parameters=parameters)
    p.set_calculator(calc)
    # E0 = p.get_potential_energy()
    natom = len(p)
    vol = p.get_volume()
    jacob = (vol/natom)**(1.0/3.0) * natom**0.5
    # mode = np.zeros((len(p)+3, 3))
    # mode = vrand(mode)
    f = open('tmode.zf')
    mode = pick.load(f)
    f.close()
    cellt = p.get_cell() + np.dot(p.get_cell(), mode[-3:]/jacob)
    p.set_cell(cellt, scale_atoms=True)
    p.set_positions(p.get_positions() + mode[:-3])
    d = lanczos.lanczos_atoms(p, mode=mode, rotationMax=4,
                              ss=True, phi_tol=15)
    dyn = FIRE(d, dt=0.1, maxmove=0.2, dtmax=0.2)
    try:
        dyn.run(fmax=0.05, steps=2000)
        E1 = p.get_potential_energy()
        write("CDCAR.vasp", d.R0, format='vasp', direct=True)
        pcell = set_cell_from_vasp("CDCAR.vasp")
        pcell.set_e(E1)
    except:
        pcell = set_cell_from_vasp("PRESAD.vasp")
        pcell.set_e(151206)

    f = open('pcell.bin', 'w')
    pick.dump(pcell, f)
    f.close()
    return pcell


if __name__ == '__main__':
    rundim_lammps()
