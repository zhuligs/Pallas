#!/usr/bin/python -u 
# 

import sys, os, time
from copy import deepcopy as cp
import numpy as np
import cPickle as pick

# local
import itin
import itdbase
import sdata
from util import vunit, vrand
from itin import npop
from itin import instep as total_step
from zfunc import write_cell_to_vasp, set_cell_from_vasp



def initrun():
    for istep in range(itin.instep):
        stepdata = []
        for i in range(itin.npop):
            stepdata.append(itdbase.ConnectObj())
        sdata.evodata.append(stepdata)

    for ip in range(itin.npop):
        cdir = 'Cal' + str(ip)
        os.system('mkdir -p ' + cdir)


def foo(xend, istep):
    energy_xend = xend.get_e()
    stepdata = sdata.evodata[istep]
    for ip in range(itin.npop):
        cdir = 'Cal' + str(ip)
        tmode = np.zeros((itin.nat + 3, 3))
        tmode = vrand(tmode)
        tmode = vunit(tmode)
        f = open(cdir + '/tmode.zf', 'w')
        pick.dump(tmode, f)
        f.close()
        stepdata[ip].leftmin = cp(xend)
        preminima = apply_mode(xend, tmode)
        stepdata[ip].premin = cp(preminima)
        # opt preminima
        # f = open(cdir + 'premin.zf', 'w')
        # pick.dump(preminima, f)
        # f.close()
    
    get_optfile_ready(istep)
    jobids = pushjob(istep) 
    if checkjob(jobids) == 0:
        pulljob('opt', istep)
    else:
        return 100


    for ip in range(itin.npop):
        print 'energy', stepdata[ip].rightmin.get_e()
    sys.exit(1)

    for ip in range(itin.npop):
        # set stepdata[ip].rightmin
        cdir = 'Cal' + str(ip)
        (sadcell, sadmode) = interpolatesad(stepdata[ip].leftmin, stepdata[ip].rightmin)
        stepdata[ip].sadmode = cp(sadmode)
        f = open(cdir + 'sadmode.zf', 'w')
        pick.dump(stepdata)
        f.close()

    get_dimfile_ready()
    jobids = pushjob()
    if checkjob(jobids) == 0:
        pulljob()
    else:
        return 100


    # RETURN stepdata

 ##############  END FOO  ######################   

    # for ip in range(itin.npop):
    #     stepdata[ip].rightmin = cp(saddle_points[ip])

    # # connect the left and right
    # for ip in range(itin.npop):
    #     # interpolate the saddle 



    # # local opt
    # get_opt_ready()
    # jobids = pushjob()

    # finished = False
    # while not finished:
    #     if check_job(jobids):
    #         locale_minima = pulljob()
    #         finished = True
    #         print 'MINIMA JOB COMPLETE'
    #     time.sleep(10)

    # stepdata[ip].rightmin = cp(locale_minima[ip])


def interpolatesad(leftmin, rightmin):
    write_cell_to_vasp(leftmin, 'leftpos')
    write_cell_to_vasp(rightmin, 'rightpos')
    os.system(itin.makenebcom + ' leftpos rightpos 1')
    nebcell = set_cell_from_vasp('01/POSCAR')
    neblat = nebcell.get_lattice()
    tmode = np.zeros((itin.nat+3, 3))
    tmode[:-3] = rightmin.get_positions() - leftmin.get_positions()
    tcell = cp(leftmin)
    tcell.set_lattice(neblat)
    return(tcell, tmode)


def apply_mode(xcell, mode):
    natom = itin.nat
    vol = xcell.get_volume()
    jacob = (vol / natom)**(1.0/3.0) * natom**0.5
    latt = xcell.get_lattice() + np.dot(xcell.get_lattice(), mode[-3:] / jacob)
    xcell.set_lattice(latt)
    xcell.set_positions(xcell.get_positions() + mode[:-3])
    return xcell


def get_optfile_ready(istep):
    stepdata = sdata.evodata[istep]
    for ip in range(itin.npop):
        cdir = 'Cal' + str(ip)
        pcar = cdir + '/POSCAR'
        write_cell_to_vasp(stepdata[ip].premin, pcar)
        os.system('cp pbs_opt.sh ' + cdir + '/pbs.sh')
        # get the incar, potcar, kpoints ready in the pbs.sh file




def pushjob(istep):
    jobids = []
    for ip in range(itin.npop):
        cdir = 'Cal' + str(ip)
        jbuff = os.popen('cd ' + cdir + '; qsub pbs.sh').read()
        jid = jbuff.strip()
        jobids.append(jid)
    return jobids


def checkids(jobids):
    while True:
        jstat = os.system("squeue > qbuff")
        if jstat == 0:
            break
        else:
            time.sleep(3)

    jbuff = []
    with open("qbuff") as f:
        for line in f:
            jbuff.append(line)

    if len(jbuff) > 0:
        reminds = []
        try:
            for x in jbuff[1:]:
                reminds.append(x.split()[0])
        except:
            print 'except 1'
            return True 
    else:
        print 'except 2'
        return True
    finished = True
    for id in jobids:
        if id in reminds:
            print 'id in reminds'
            finished = False
    return finished


def checkjob(jobids):
    finished = False
    while not finished:
        finished = checkids(jobids)
        print time.ctime(), 'JOB FINISHED:', finished
        time.sleep(3)
        if os.path.isfile('CSTOP'):
            os.system('rm -f CSTOP')
            print 'CSTOP'
            for ii in idpool:
                print 'qdel ' + ii
                os.system('qdel ' + ii)
            return 100
    return 0


def pulljob(otyp, istep):
    if otyp == 'opt':
        for ip in range(itin.npop):
            cdir = 'Cal' + str(ip)
            pcell = set_cell_from_vasp(cdir + '/CONTCAR')
            e = float(os.popen("awk '/free  energy/{print $5}' " + cdir +
                "/OUTCAR|tail -1").read())
            h = itin.press * pcell.get_volume() / 1602.2 + e 
            pcell.set_e(h)
            sdata.evodata[istep][ip].rightmin = cp(pcell)

    return 0



# def evolution():
#     for istep in range(total_step):
#         foo(reac, istep)

#         for i in range(npop):
#             pass()



def test():
    initrun()
    f = open('xend.bin')
    xend = pick.load(f)
    f.close()
    foo(xend, 0)




if __name__ == '__main__':
    test()
            
























