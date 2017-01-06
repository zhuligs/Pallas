#!/usr/bin/env python

'''
d v job : Dimer VASP job
'''

from tsase.dimer import ssdimer
from tsase.neb.util import vunit, vrand
from ase.io import read, write
import os
import sys
import numpy as np
from tsase.calculators.vasp_ext import Vasp
import cPickle as pick


# read geometry and set calculator
p = read('TSCELL',format='vasp')
#p = read('CdSe_hex',format='vasp')

calc = Vasp(prec = 'Normal',
            ediff = 1e-5,
            kpts = (3,3,3),
            gamma= True,
            lcharg = False,
            npar = 4,
            encut= 520,
            ismear = 0,
            sigma  = 0.08,
            isym = 0,
              )
p.set_calculator(calc)

pstress = np.zeros((3,3))
pstress[0][0] = .093622
pstress[1][1] = .093622
pstress[2][2] = .093622

# calculate the jacobian for the tangent
natom = len(p)
vol   = p.get_volume()
jacob = (vol/natom)**(1.0/3.0) * natom**0.5

# use the known final state to set the initial mode direction
#pfin = read('CdSe_sq',format='vasp')
#mode = np.zeros((len(p)+3,3))
#mode[-3:] = pfin.get_cell()-p.get_cell()
#icell= np.linalg.inv(p.get_cell())
#mode[-3:] = np.dot(icell, mode[-3:]) * jacob

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
f = open('tmode')
mode = pick.load(f)
f.close()

# displace along the initial mode direction
mode = vunit(mode)
cellt = p.get_cell()+np.dot(p.get_cell(), mode[-3:]/jacob)
p.set_cell(cellt, scale_atoms=True)
p.set_positions(p.get_positions() + mode[:-3])

E0 = p.get_potential_energy()

# set a ssdimer_atoms object
d = ssdimer.SSDimer_atoms(p, mode = mode, rotationMax = 4, phi_tol=15, ss = True, express=pstress)

#################################################
# use quickmin optimizer in ssdimer
d.search(minForce = 0.01, movie = "dimer2.movie", interval = 20 )
#################################################
# use MDMin optimizer in ase
#dyn = MDMin(d)
#dyn.run(fmax=0.01)
##################################################
# use FIRE optimizer in ase
#dyn = FIRE(d,dt=0.05, maxmove=0.1, dtmax=0.05)
#dyn.run(fmax=0.01)
#################################################

write("dimer1.vasp", d.R0, format='vasp')
E1   = p.get_potential_energy()
print "TTENERGY ", E1
# print "barrier:", E1 - E0
