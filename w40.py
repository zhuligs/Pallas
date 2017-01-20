#!/usr/bin/python -u
#

import sys, os, time
from copy import deepcopy as cp
import numpy as np
import cPickle as pick

# local
import fppy
import itin
import itdbase
import sdata
from util import vunit, vrand
from itin import npop
from itin import instep as total_step
from zfunc import write_cell_to_vasp, set_cell_from_vasp, getx


def initrun():
    f = open('prod.bin')
    sdata.product = pick.load(f)
    f.close()
    f = open('reac.bin')
    sdata.reactant = pick.load(f)
    f.close()
    lp = sdata.product.get_lfp()
    ls = sdata.product.get_lfp()
    sdata.types = sdata.reactant.get_types()
    (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, lp, ls)
    for istep in range(itin.instep):
        stepdata = []
        for i in range(itin.npop):
            stepdata.append(itdbase.ConnectObj())
        sdata.evodata.append(stepdata)

    for ip in range(itin.npop):
        cdir = 'Cal' + str(ip)
        os.system('mkdir -p ' + cdir)
        sdata.pbests.append(cp(sdata.reactant))
        sdata.fitpbest.append(dist)


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
        stepdata[ip].psomode = cp(tmode)
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
    # sys.exit(1)

    for ip in range(itin.npop):
        # set stepdata[ip].rightmin
        cdir = 'Cal' + str(ip)
        (sadcell, sadmode) = interpolatesad(stepdata[ip].leftmin,
                                            stepdata[ip].rightmin)
        stepdata[ip].sadmode = cp(sadmode)
        stepdata[ip].presad = cp(sadcell)
        f = open(cdir + '/sadmode.zf', 'w')
        pick.dump(sadmode, f)
        f.close()

    get_dimfile_ready(istep)
    # sys.exit(1)
    jobids = pushjob(istep)
    if checkjob(jobids) == 0:
        pulljob('sad', istep)
    else:
        return 100

    # DEBUG
    for ip in range(itin.npop):
        print 'IP ', ip
        print stepdata[ip].leftmin.get_e()
        print stepdata[ip].saddle.get_e()
        print stepdata[ip].rightmin.get_e()
    print '###### END #####'

    # RETURN stepdata
    # #############  END FOO  ######################

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

    # use the lattice of leftmin, and atomic positions of neb


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


def get_dimfile_ready(istep):
    stepdata = sdata.evodata[istep]
    for ip in range(itin.npop):
        cdir = 'Cal' + str(ip)
        pcar = cdir + '/PRESAD.vasp'
        write_cell_to_vasp(stepdata[ip].presad, pcar)
        os.system('cp dvjob.py ' + cdir)
        os.system('cp pbs_dim.sh ' + cdir + '/pbs.sh')


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
            for ii in jobids:
                print 'qdel ' + ii
                os.system('qdel ' + ii)
            return 100
    return 0


def pulljob(otyp, istep):
    if otyp == 'opt':
        for ip in range(itin.npop):
            cdir = 'Cal' + str(ip)
            pcell = set_cell_from_vasp(cdir + '/CONTCAR')
            try:
                e = float(os.popen("awk '/free  energy/{print $5}' " + cdir +
                    "/OUTCAR|tail -1").read())
                h = itin.press * pcell.get_volume() / 1602.2 + e
            except:
                h = 31118.
            pcell.set_e(h)
            sdata.evodata[istep][ip].rightmin = cp(pcell)
            sdata.evodata[istep+1][ip].leftmin = cp(pcell)
    elif otyp == 'sad':
        for ip in range(itin.npop):
            cdir = 'Cal' + str(ip)
            pcell = set_cell_from_vasp(cdir + '/CONTCAR')
            try:
                e = float(os.popen("grep DIMERENERGY " + cdir + "/dvjob.out |tail -1").read().split()[1])
                h = itin.press * pcell.get_volume() / 1602.2 + e
            except:
                h = 31118.
            pcell.set_e(h)
            sdata.evodata[istep][ip].saddle = cp(pcell)
    else:
        print '# ZLOG: otyp error'

    return 0


def evolution():
    reac = sdata.reactant
    prod = sdata.product
    # fingerprint of product: fpp
    fpp = prod.get_lfp()
    # for istep in range(total_step):
    istep = 0
    foo(reac, istep)
    f = open('sdata.bin', 'w')
    pick.dump(sdata.evodata, f)
    f.close()
    # rank the fitness
    dists = []
    for ip in range(itin.npop):
        fpx = sdata.evodata[istep][ip].rightmin.get_lfp()
        (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpp, fpx)
        print "# ZLOG: STEP %d IP %d DIST %7.4E" % (istep, ip, dist)
        dists.append(dist)

        if dist < sdata.fitpbest[ip]:
            sdata.fitpbest[ip] = dist
            sdata.pbests[ip] = cp(sdata.evodata[istep][ip].rightmin)
    mindist = min(dists)
    minid = dists.index(mindist)
    # if sdata.gbest is None:
    sdata.gbest = cp(sdata.evodata[istep][minid].rightmin)
    sdata.bestdist = mindist
    f = open('sdata.bin', 'w')
    pick.dump(sdata.evodata, f)
    f.close()
    # else:
    #     if mindist < sdata.bestdist:
    #         sdata.gbest = cp(sdata.evodata[istep][minid].rightmin)
    #         sdata.bestdist = mindist

    for istep in range(1, total_step):
        stepdata = sdata.evodata[istep]
        for ip in range(itin.npop):
            psomode = gen_psomode(stepdata[ip].leftmin, istep, ip)
            cdir = 'Cal' + str(ip)
            f = open(cdir + '/pmode.zf', 'w')
            pick.dump(psomode, f)
            f.close()
            preminima = apply_mode(stepdata[ip].leftmin, psomode)
            stepdata[ip].premin = cp(preminima)
            stepdata[ip].psomode = cp(psomode)

        get_optfile_ready(istep)
        jobids = pushjob(istep)
        if checkjob(jobids) == 0:
            pulljob('opt', istep)
        else:
            return 100

        for ip in range(itin.npop):
            print 'energy', stepdata[ip].rightmin.get_e()

        for ip in range(itin.npop):
            cdir = 'Cal' + str(ip)
            (sadcell, sadmode) = interpolatesad(stepdata[ip].leftmin,
                                                stepdata[ip].rightmin)
            stepdata[ip].sadmode = cp(sadmode)
            stepdata[ip].presad = cp(sadcell)
            f = open(cdir + '/sadmode.zf', 'w')
            pick.dump(sadmode, f)
            f.close()

        get_dimfile_ready(istep)
        # sys.exit(1)
        jobids = pushjob(istep)
        if checkjob(jobids) == 0:
            pulljob('sad', istep)
        else:
            return 100

        # DEBUG
        for ip in range(itin.npop):
            print 'IP ', ip
            print stepdata[ip].leftmin.get_e()
            print stepdata[ip].saddle.get_e()
            print stepdata[ip].rightmin.get_e()

        dists = []
        for ip in range(itin.npop):
            fpx = sdata.evodata[istep][ip].rightmin.get_lfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpp, fpx)
            print "# ZLOG: STEP %d IP %d DIST %7.4E" % (istep, ip, dist)
            dists.append(dist)

            if dist < sdata.fitpbest[ip]:
                sdata.fitpbest[ip] = dist
                sdata.pbests[ip] = cp(sdata.evodata[istep][ip].rightmin)
        mindist = min(dists)
        minid = dists.index(mindist)
        if mindist < sdata.bestdist:
            sdata.gbest = cp(sdata.evodata[istep][minid].rightmin)
            sdata.bestdist = mindist

        f = open('sdata.bin', 'w')
        pick.dump(sdata.evodata, f)
        f.close()


def gen_psomode(x0, istep, ip):
    c1 = 2.0
    c2 = 2.0
    # c3 = 2.0
    w = 0.9 - 0.5 * (istep + 1) / itin.instep
    (r1, r2, r3) = np.random.rand(3)
    v0 = sdata.evodata[istep-1][ip].psomode
    pbest = sdata.pbests[ip]
    gbest = sdata.product
    difp = getx(pbest, x0)
    difg = getx(gbest, x0)
    print 'v0', v0
    print 'difp', difp
    print 'difg', difg
    v = v0 * w + c1 * r1 * difp + c2 * r2 * difg
    return v


def main():
    initrun()
    evolution()


def test1():
    initrun()
    reac = sdata.reactant
    prod = sdata.product
    print reac.get_e()
    print prod.get_e()
    f = open('sdata.bin')
    sdata.evodata = pick.load(f)
    f.close()
    for ip in range(itin.npop):
        print sdata.evodata[0][ip].rightmin.get_e()
        print sdata.evodata[0][ip].psomode



def test():
    initrun()
    f = open('reac.bin')
    xend = pick.load(f)
    f.close()
    foo(xend, 0)


if __name__ == '__main__':
    main()

























