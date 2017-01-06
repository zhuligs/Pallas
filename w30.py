#!/usr/bin/python -u 
# 

import sys, os, time
from copy import deepcopy as cp
import numpy as np
import cPickle as pick

# local
import itin
import itdbase
from util import vunit, vrand
from itin import npop
from itin import instep as total_step


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
        # (xsaddle, xvel) = gen_rsaddle(xend)
        # generate the input file for tsase job
        tmode = np.zeros((itin.nat + 3, 3))
        tmode = vrand(mode)
        tmode = vunit(mode)
        f = open(cdir + 'tmode', 'w')
        pick.dump(tmode, f)
        f.close()
        stepdata[ip].leftmin = cp(xend)
    
    get_sad_ready()
    jobids = pushjob()  

    finished = False
    while not finished:
        if check_job(jobids):
            saddle_points = pulljob()
            finished = True
            print 'SADDLE JOB COMPLETE'
        time.sleep(10)

    # saddle points
    # check_identical(saddle_points)

    for ip in range(itin.npop):
        stepdata[ip].saddle = cp(saddle_points[ip])

    # local opt
    get_opt_ready()
    jobids = pushjob()

    finished = False
    while not finished:
        if check_job(jobids):
            locale_minima = pulljob()
            finished = True
            print 'MINIMA JOB COMPLETE'
        time.sleep(10)

    stepdata[ip].rightmin = cp(locale_minima[ip])


def evolution():
    for istep in range(total_step):
        foo(reac, istep)

        for i in range(npop):
            
























