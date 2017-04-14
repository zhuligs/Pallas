#!/usr/bin/env python

'''
d v job : Dimer VASP job
'''

from tsase.dimer import ssdimer
from tsase.neb.util import vunit, vrand
from ase.io import read, write
# import os
# import sys
import numpy as np
from tsase.calculators.vasp_ext import Vasp
import cPickle as pick

import itin
from zfunc import set_cell_from_vasp


# read geometry and set calculator
p = read('PRESAD.vasp', format='vasp')
# p = read('CdSe_hex',format='vasp')

ccwd = os.getcwd().split('/')[-1]
if 'x' in ccwd:
    calc = Vasp(prec='Normal',
                ediff=1e-5,
                kpts=(3, 3, 3),
                gamma=True,
                lcharg=False,
                npar=4,
                ismear=0,
                sigma=0.08,
                isym=0,
                )
else:
    calc = Vasp(prec='Normal',
                ediff=1e-5,
                kpts=(3, 3, 3),
                gamma=True,
                lcharg=False,
                npar=4,
                ismear=0,
                sigma=0.08,
                isym=0,
                )

p.set_calculator(calc)

pst = itin.press / 1602.176487

pstress = np.zeros((3, 3))
pstress[0][0] = pst
pstress[1][1] = pst
pstress[2][2] = pst

# calculate the jacobian for the tangent
natom = len(p)
vol = p.get_volume()
jacob = (vol / natom)**(1.0 / 3.0) * natom**0.5

# use the known final state to set the initial mode direction
# pfin = read('CdSe_sq',format='vasp')
# mode = np.zeros((len(p)+3,3))
# mode[-3:] = pfin.get_cell()-p.get_cell()
# icell= np.linalg.inv(p.get_cell())
# mode[-3:] = np.dot(icell, mode[-3:]) * jacob

#######################################
# set the initial mode randomly
# mode = np.zeros((len(p)+3,3))
# mode = vrand(mode)
# constrain 3 redundant freedoms
# mode[0]    *=0
# mode[-3,1:]*=0
# mode[-2,2] *=0
#######################################

# read mode file
f = open('mode.zf')
mode = pick.load(f)
f.close()

# displace along the initial mode direction
mode = vunit(mode)
cellt = p.get_cell() + np.dot(p.get_cell(), mode[-3:] / jacob)
p.set_cell(cellt, scale_atoms=True)
p.set_positions(p.get_positions() + mode[:-3])

E0 = p.get_potential_energy()

# set a ssdimer_atoms object
d = ssdimer.SSDimer_atoms(p, mode=mode, rotationMax=4,
                          phi_tol=15, ss=True, express=pstress)

#################################################
# use quickmin optimizer in ssdimer
# d.search(minForce=0.01, movie="dimer2.movie", interval=20)
d.search(minForce=0.01, maxForceCalls=3000)
#################################################
# use MDMin optimizer in ase
# dyn = MDMin(d)
# dyn.run(fmax=0.01)
##################################################
# use FIRE optimizer in ase
# dyn = FIRE(d,dt=0.05, maxmove=0.1, dtmax=0.05)
# dyn.run(fmax=0.01)
#################################################

write("dimer1.vasp", d.R0, format='vasp', direct=True)
E1 = p.get_potential_energy()
print "TTENERGY ", E1
# print "barrier:", E1 - E0

try:
    pcell = set_cell_from_vasp("dimer1.vasp")
    h = itin.press * pcell.get_volume() / 1602.176487 + E1
    pcell.set_e(h)
except:
    pcell = set_cell_from_vasp("PRESAD.vasp")
    h = 31118.
    pcell.set_e(h)

f = open('pcell.bin', 'w')
pick.dump(pcell, f)
f.close()
