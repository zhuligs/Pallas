import numpy as np
import itdbase
from itdbase import Cell
import itin
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


def gopt(xcell, mode):
    if itin.interface == 'lammps':
        return gopt_lammps(xcell, mode)
    elif itin.interface == 'vasp':
        return gopt_vasp(xcell, mode)
    else:
        print 'ERROR: WRONG INTERFACE'
        sys.exit(1)


def gopt_vasp(xcell, mode):
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


def gopt_lammps(xcell, mode):
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
    if itin.interface == 'lammps':
        return rundim_lammps(xcell, mode)
    elif itin.interface == 'vasp':
        return rundim_vasp(xcell, mode)
    else:
        print 'ERROR: WRONG INTERFACE'
        sys.exit(1)


def rundim_vasp(xcell, mode):
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


def rundim_ts(xcell, mode):
    write_cell_to_vasp('TSCELL', 'w')
    f = open('tmode', 'w')
    pick.dump(mode, 'f')
    f.close()
    os.system('python -u dvjob.py > zout')
    os.system('rm WAVECAR')
    e = float(os.popen("awk '/TTENERGY/{print $2}' zout").read())
    pcell = set_cell_from_vasp('dimer1.con')
    h = itin.press * pcell.get_volume() / 1602.2 + e
    pcell.set_e(h)
    return pcell



def rundim_lammps(xcell, mode):
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
    try:
        mode = vunit(mode)
    except:
        mode = z_rmode()
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
    try:
        typt = np.array(buff[5], int)
    except:
        del(buff[5])
        typt = np.array(buff[5], int)
    pos = np.array(buff[7:7 + itin.nat], float)
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
    for iz in itin.znucl:
        f.write(itdbase.atom_data[iz][1])
        f.write('  ')
    f.write('\n')
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
    try:
        mode = vunit(mode)
    except:
        mode = np.zeros((itin.nat + 3, 3))
    return mode


def z_rmode():
    mode = np.zeros((itin.nat + 3, 3))
    mode = vrand(mode)
    mode = vunit(mode)
    return mode
