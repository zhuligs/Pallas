#!/usr/bin/python -u
# encoding: utf-8

import sys
import os
from copy import deepcopy as cp
import numpy as np
import cPickle as pick
from pele.storage import Database

from treelib import Node, Tree
import networkx as nx

import itin
import sdata
import fppy
import itdbase
from wrapdimer import get_mode, get_0mode, get_rmode
from zfunc import gopt, rundim, set_cell_from_vasp, write_cell_to_vasp, getx
from w2 import gen_rsaddle, initrun
# from w40 import get_dimfile_ready, get_optfile_ready
from w40 import checkjob


sys.setrecursionlimit(100000)


def w20init():
    for ip in range(itin.npop):
        xdir = 'Calx' + str(ip)
        ydir = 'Caly' + str(ip)
        os.system('mkdir -p ' + xdir)
        os.system('mkdir -p ' + ydir)
        sdata.xdirs.append(xdir)
        sdata.ydirs.append(ydir)
    # for istep in range(itin.instep):
    #     stepx = []
    #     stepy = []
    #     for i in range(itin.npop):
    #         stepx.append(itdbase.Cobj())
    #         stepy.append(itdbase.Cobj())
    #     sdata.evox.append(stepx)
    #     sdata.evoy.append(stepy)


def create_stepdata():
    stepdata = []
    for i in range(itin.npop):
        stepdata.append(itdbase.Cobj())
    return stepdata


def wo():
    # istep = 0
    reace = sdata.reactant.get_e()
    print 'ZLOG, reace', reace
    sdata.reactant.add_left(-1)
    sdata.product.add_left(-1)
    sdata.xllist.append(sdata.reactant)
    sdata.yllist.append(sdata.product)
    v = sdata.reactant.get_volume() / itin.nat
    e = sdata.reactant.get_e() - reace
    sdata.G.add_node('xl0', energy=e, volume=v)
    v = sdata.product.get_volume() / itin.nat
    e = sdata.product.get_e() - reace
    sdata.G.add_node('yl0', energy=e, volume=v)
    # sdata.xslist.append(sdata.reactant)
    # sdata.yslist.append(sdata.product)
    # DATABASE

    # prepdim(istep)
    # stepx = sdata.evox[istep]
    # stepy = sdata.evoy[istep]
    stepx = create_stepdata()
    stepy = create_stepdata()
    for ip in range(itin.npop):
        xdir = sdata.xdirs[ip]
        ydir = sdata.ydirs[ip]
        xpcar = xdir + '/PRESAD.vasp'
        ypcar = ydir + '/PRESAD.vasp'
        write_cell_to_vasp(sdata.reactant, xpcar)
        write_cell_to_vasp(sdata.product, ypcar)
        xmode = get_rmode()
        ymode = get_rmode()
        stepx[ip].v = cp(xmode)
        stepy[ip].v = cp(ymode)
        f = open(xdir + '/mode.zf', 'w')
        pick.dump(xmode, f)
        f.close()
        f = open(ydir + '/mode.zf', 'w')
        pick.dump(ymode, f)
        f.close()
        os.system('cp pbs_dim.sh ' + xdir + '/pbs.sh')
        os.system('cp pbs_dim.sh ' + ydir + '/pbs.sh')
    xkeep = [0] * itin.npop
    ykeep = [0] * itin.npop
    jobids = pushjob(xkeep, ykeep)
    if checkjob(jobids) == 0:
        (xsets, ysets) = pulljob(xkeep, ykeep)
    else:
        return 100

    itry = 0
    while itry < 3:
        itry += 1
        xkeep = get_keep(xsets)
        ykeep = get_keep(ysets)

        if not((0 in xkeep) or (0 in ykeep)):
            print 'ZLOG: no keep in xkeep and ykeep'
            break

        for ip in range(itin.npop):
            if xkeep[ip] == 0:
                xdir = sdata.xdirs[ip]
                xpcar = xdir + '/PRESAD.vasp'
                write_cell_to_vasp(sdata.reactant, xpcar)
                xmode = get_rmode()
                stepx[ip].v = cp(xmode)
                f = open(xdir + '/mode.zf', 'w')
                pick.dump(xmode, f)
                f.close()
                os.system('cp pbs_dim.sh ' + xdir + '/pbs.sh')

            if ykeep[ip] == 0:
                ydir = sdata.ydirs[ip]
                ypcar = ydir + '/PRESAD.vasp'
                write_cell_to_vasp(sdata.product, ypcar)
                ymode = get_rmode()
                stepy[ip].v = cp(ymode)
                f = open(ydir + '/mode.zf', 'w')
                pick.dump(ymode, f)
                f.close()
                os.system('cp pbs_dim.sh ' + xdir + '/pbs.sh')
        jobids = pushjob(xkeep, ykeep)
        if checkjob(jobids) == 0:
            (xsets_tmp, ysets_tmp) = pulljob(xkeep, ykeep)
        else:
            return 100
        for ip in range(itin.npop):
            if xkeep[ip] == 0:
                xsets[ip] = cp(xsets_tmp[ip])
            if ykeep[ip] == 0:
                ysets[ip] = cp(ysets_tmp[ip])


    # preploc(istep, xsets, ysets)
    # stepx = sdata.evox[istep]
    # stepy = sdata.evoy[istep]
    for ip in range(itin.npop):
        stepx[ip].sad = cp(xsets[ip])
        stepy[ip].sad = cp(ysets[ip])
        xid = update_iden(sdata.xslist, stepx[ip].sad)
        yid = update_iden(sdata.xslist, stepy[ip].sad)
        stepx[ip].sad.set_iden(xid)
        stepy[ip].sad.set_iden(yid)
        stepx[ip].sad.set_sm('S')
        stepy[ip].sad.set_sm('S')
        stepx[ip].sad.add_left(sdata.reactant.get_iden())
        stepy[ip].sad.add_left(sdata.product.get_iden())
        sdata.reactant.add_right(xid)
        sdata.product.add_right(yid)
        sdata.xslist.append(stepx[ip].sad)
        sdata.yslist.append(stepy[ip].sad)

        xdir = sdata.xdirs[ip]
        ydir = sdata.ydirs[ip]
        xpcar = xdir + '/POSCAR'
        ypcar = ydir + '/POSCAR'
        write_cell_to_vasp(stepx[ip].sad, xpcar)
        write_cell_to_vasp(stepy[ip].sad, ypcar)
        os.system('cp pbs_opt.sh ' + xdir + '/pbs.sh')
        os.system('cp pbs_opt.sh ' + ydir + '/pbs.sh')
        xnode_name = 'xs' + str(xid)
        xvol = xsets[ip].get_volume() / itin.nat
        xe = xsets[ip].get_e() - reace
        sdata.G.add_node(xnode_name, energy=xe, volume=xvol)
        ynode_name = 'ys' + str(yid)
        yvol = ysets[ip].get_volume() / itin.nat
        ye = xsets[ip].get_e() - reace
        sdata.G.add_node(ynode_name, energy=ye, volume=yvol)
        sdata.G.add_edge('xl0', xnode_name)
        sdata.G.add_edge('yl0', ynode_name)

    del(xsets)
    del(ysets)

    xkeep = [0] * itin.npop
    ykeep = [0] * itin.npop
    jobids = pushjob(xkeep, ykeep)
    if checkjob(jobids) == 0:
        (xsets, ysets) = pulljob(xkeep, ykeep)
    else:
        return 100

    # prepso(istep, xsets, ysets)
    # stepx = sdata.evox[istep]
    # stepy = sdata.evoy[istep]
    for ip in range(itin.npop):
        stepx[ip].loc = cp(xsets[ip])
        stepy[ip].loc = cp(ysets[ip])
        xid = update_iden(sdata.xllist, stepx[ip].loc)
        yid = update_iden(sdata.yllist, stepy[ip].loc)
        stepx[ip].loc.set_iden(xid)
        stepy[ip].loc.set_iden(yid)
        stepx[ip].loc.set_sm('M')
        stepy[ip].loc.set_sm('M')
        stepx[ip].loc.add_left(stepx[ip].sad.get_iden())
        stepy[ip].loc.add_left(stepy[ip].sad.get_iden())

        stepx[ip].sad.add_right(xid)
        stepy[ip].sad.add_right(yid)
        print "ZLOG: INIT STEP, IP %4d X SAD EN: %8.7E, X LOC EN: %8.7E" %\
              (ip, stepx[ip].sad.get_e(), stepx[ip].loc.get_e())
        print "ZLOG: INIT STEP, IP %4d Y SAD EN: %8.7E, Y LOC EN: %8.7E" %\
              (ip, stepy[ip].sad.get_e(), stepy[ip].loc.get_e())
        sdata.xllist.append(stepx[ip].loc)
        sdata.yllist.append(stepy[ip].loc)

        sdata.pbestx.append(stepx[ip].loc)
        sdata.pbesty.append(stepy[ip].loc)
        xnode_name = 'xl' + str(xid)
        xvol = xsets[ip].get_volume() / itin.nat
        xe = xsets[ip].get_e() - reace
        sdata.G.add_node(xnode_name, energy=xe, volume=xvol)
        ynode_name = 'yl' + str(yid)
        yvol = ysets[ip].get_volume() / itin.nat
        ye = xsets[ip].get_e() - reace
        sdata.G.add_node(ynode_name, energy=ye, volume=yvol)
        pxnode_name = 'xs' + str(stepx[ip].sad.get_iden())
        pynode_name = 'ys' + str(stepy[ip].sad.get_iden())
        sdata.G.add_edge(pxnode_name, xnode_name)
        sdata.G.add_edge(pynode_name, ynode_name)

    del(xsets)
    del(ysets)
    dumpdata()

    # stepx = sdata.evox[istep]
    # stepy = sdata.evoy[istep]
    xydist = []
    for ix in range(itin.npop):
        fpx = stepx[ix].loc.get_sfp()
        # ex = Xsad[ix].get_e() - reace
        # ex = get_barrier(sdata.xllist, sdata.xslist, sdata.reactant,
        #                  stepx[ix].loc)
        for iy in range(itin.npop):
            fpy = stepy[iy].loc.get_sfp()
            # ey = Ysad[iy].get_e() - reace
            # ey = get_barrier(sdata.yllist, sdata.yslist, sdata.product,
            #                  stepy[iy].loc)
            # ee = max(ex, ey)
            ee = 0.0
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            # mdist for multiobjective opt
            if dist < 1e-4:
                xdist = 1e-4
            else:
                xdist = dist
            mdist = np.log10(xdist) + ee
            # print 'ZLOG: mdist, dist, log(dist), ee',\
            #       mdist, dist, np.log(dist), ee
            xydist.append((mdist, (ix, iy), dist, ee))

            if dist < itin.dist:
                xnode_name = 'xl' + str(stepx[ix].loc.get_iden())
                ynode_name = 'yl' + str(stepy[iy].loc.get_iden())
                sdata.G.add_edge(xnode_name, ynode_name)

    xydistSort = sorted(xydist, key=lambda x: x[2])
    sdata.bestdist = xydistSort[0][2]
    sdata.bestmdist = xydistSort[0][0]
    # sdata.gbestx = cp(Xloc[xydistSort[0][1][0]])
    # sdata.gbesty = cp(Yloc[xydistSort[0][1][1]])
    sdata.gbestx = cp(stepx[xydistSort[0][1][0]].loc)
    sdata.gbesty = cp(stepy[xydistSort[0][1][1]].loc)
    print "ZLOG: INIT STEP, bestDist: %8.7E, bestmD: %8.7E, X-Y: %4d %4d" %\
          (xydistSort[0][2], xydistSort[0][0],
           xydistSort[0][1][0], xydistSort[0][1][1])

    # update pdist x
    for ix in range(itin.npop):
        xytdist = []
        # fpx = Xloc[ix].get_sfp()
        fpx = stepx[ix].loc.get_sfp()
        # ex = Xsad[ix].get_e() - reace
        # ex = get_barrier(sdata.xllist, sdata.xslist, reac, Xloc[ix])
        # ex = get_barrier(sdata.xllist, sdata.xslist, sdata.reactant,
                         # stepx[ix].loc)
        for iy in range(len(sdata.yllist)):
            fpy = sdata.yllist[iy].get_sfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            # ey = sdata.yslist[iy].get_e() - reace
            # ey = get_barrier(sdata.yllist, sdata.yslist, sdata.product,
                             # sdata.yllist[iy])
            # ee = max(ex, ey)
            ee = 0
            if dist < 1e-4:
                xdist = 1e-4
            else:
                xdist = dist
            # mdist = np.log10(xdist) + ee
            # print 'ZLOG: mdist, dist, log(dist), ee',\
            #       mdist, dist, np.log(dist), ee
            xytdist.append(dist)
        xytdistb = sorted(xytdist)[0]
        sdata.pdistx.append(xytdistb)

    # update pdist y
    for iy in range(itin.npop):
        yxtdist = []
        # fpy = Yloc[iy].get_sfp()
        fpy = stepy[iy].loc.get_sfp()
        # ey = Ysad[iy].get_e() - reace
        # ey = get_barrier(sdata.yllist, sdata.yslist, prod, Yloc[iy])
        # ey = get_barrier(sdata.yllist, sdata.yslist, sdata.product,
                         # stepy[iy].loc)
        for ix in range(len(sdata.xllist)):
            fpx = sdata.xllist[ix].get_sfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            # ex = sdata.xslist[ix].get_e() - reace
            # ex = get_barrier(sdata.xllist, sdata.xslist, sdata.reactant,
                             # sdata.xllist[ix])
            # ee = max(ex, ey)
            ee = 0
            if dist < 1e-4:
                xdist = 1e-4
            else:
                xdist = dist
            # mdist = np.log10(xdist) + ee
            # print 'ZLOG: mdist, dist, log(dist), ee',\
            #       mdist, dist, np.log(dist), ee
            yxtdist.append(dist)
        yxtdistb = sorted(yxtdist)[0]
        sdata.pdisty.append(yxtdistb)

    # iratio = int(itin.psoration * itin.npop)
    # if iratio >= itin.npop:
    #     iratio = itin.npop - 1
    for ip in range(itin.npop):
        sdata.ifpsox[ip] = True
        sdata.ifpsoy[ip] = True
    dumpdata()

    dataf = 'stepx_0.bin'
    f = open(dataf, 'w')
    pick.dump(stepx, f)
    f.close()
    dataf = 'stepy_0.bin'
    f = open(dataf, 'w')
    pick.dump(stepy, f)
    f.close()
    del(stepx)
    del(stepy)

    showpath()


def get_keep(xsets):
    dij = np.zeros((itin.npop, itin.npop))
    for i in range(itin.npop - 1):
        fpi = xsets[i].get_sfp()
        for j in range(i + 1, itin.npop):
            fpj = xsets[j].get_sfp()
            (d, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpi, fpj)
            dij[i][j] = d
            dij[j][i] = d

    xkeep = [1] * itin.npop

    for i in range(itin.npop - 1):
        for j in range(i + 1, itin.npop):
            if dij[i][j] < itin.dist:
                xkeep[j] = 0
    return xkeep


def dumpdata():
    # f = open('evox.bin', 'w')
    # pick.dump(sdata.evox, f)
    # f.close()
    # f = open('evoy.bin', 'w')
    # pick.dump(sdata.evoy, f)
    # f.close()
    f = open('xslist.bin', 'w')
    pick.dump(sdata.xslist, f)
    f.close()
    f = open('xllist.bin', 'w')
    pick.dump(sdata.xllist, f)
    f.close()
    f = open('yslist.bin', 'w')
    pick.dump(sdata.yslist, f)
    f.close()
    f = open('yllist.bin', 'w')
    pick.dump(sdata.yllist, f)
    f.close()


# @profile
def woo():
    # init step
    reace = sdata.reactant.get_e()
    wo()
    # pso step
    xkeep = [0] * itin.npop
    ykeep = [0] * itin.npop
    for istep in range(1, itin.instep):
        # stepx = sdata.evox[istep]
        # stepy = sdata.evoy[istep]
        # pstepx = sdata.evox[istep - 1]
        # pstepy = sdata.evoy[istep - 1]
        stepx = create_stepdata()
        stepy = create_stepdata()
        dataf = 'stepx_' + str(istep - 1) + '.bin'
        f = open(dataf)
        pstepx = pick.load(f)
        f.close()
        dataf = 'stepy_' + str(istep - 1) + '.bin'
        f = open(dataf)
        pstepy = pick.load(f)
        f.close()
        for ip in range(itin.npop):
            if os.path.isfile('CSTOP'):
                dumpdata()
                return
            # xloc = xlocs[ip]
            # yloc = ylocs[ip]
            if sdata.ifpsox[ip]:
                # (xsad, vx) = gen_psaddle('x', xloc, istep, ip)
                v0 = pstepx[ip].v
                xmode = gen_pmode('x', pstepx[ip].loc, v0, istep, ip)
            else:
                # (xsad, vx) = gen_rsaddle(xloc)
                xmode = get_rmode()
            if sdata.ifpsoy[ip]:
                # (ysad, vy) = gen_psaddle('y', yloc, istep, ip)
                v0 = pstepy[ip].v
                ymode = gen_pmode('y', pstepy[ip].loc, v0, istep, ip)
            else:
                # (ysad, vy) = gen_rsaddle(yloc)
                ymode = get_rmode()

            xdir = sdata.xdirs[ip]
            ydir = sdata.ydirs[ip]
            xpcar = xdir + '/PRESAD.vasp'
            ypcar = ydir + '/PRESAD.vasp'
            write_cell_to_vasp(pstepx[ip].loc, xpcar)
            write_cell_to_vasp(pstepy[ip].loc, ypcar)
            stepx[ip].v = cp(xmode)
            stepy[ip].v = cp(ymode)
            f = open(xdir + '/mode.zf', 'w')
            pick.dump(xmode, f)
            f.close()
            f = open(ydir + '/mode.zf', 'w')
            pick.dump(ymode, f)
            f.close()
            os.system('cp pbs_dim.sh ' + xdir + '/pbs.sh')
            os.system('cp pbs_dim.sh ' + ydir + '/pbs.sh')

        jobids = pushjob(xkeep, ykeep)
        if checkjob(jobids) == 0:
            (xsets, ysets) = pulljob(xkeep, ykeep)
        else:
            return 100

        for ip in range(itin.npop):
            stepx[ip].sad = cp(xsets[ip])
            stepy[ip].sad = cp(ysets[ip])
            xid = update_iden(sdata.xslist, stepx[ip].sad)
            yid = update_iden(sdata.xslist, stepy[ip].sad)
            stepx[ip].sad.set_iden(xid)
            stepy[ip].sad.set_iden(yid)
            pstepx[ip].loc.add_right(xid)
            pstepy[ip].loc.add_right(yid)
            stepx[ip].sad.set_sm('S')
            stepy[ip].sad.set_sm('S')
            stepx[ip].sad.add_left(pstepx[ip].loc.get_iden())
            stepy[ip].sad.add_left(pstepy[ip].loc.get_iden())
            sdata.xslist.append(stepx[ip].sad)
            sdata.yslist.append(stepy[ip].sad)

            xdir = sdata.xdirs[ip]
            ydir = sdata.ydirs[ip]
            xpcar = xdir + '/POSCAR'
            ypcar = ydir + '/POSCAR'
            write_cell_to_vasp(stepx[ip].sad, xpcar)
            write_cell_to_vasp(stepy[ip].sad, ypcar)
            os.system('cp pbs_opt.sh ' + xdir + '/pbs.sh')
            os.system('cp pbs_opt.sh ' + ydir + '/pbs.sh')
            xnode_name = 'xs' + str(xid)
            xvol = stepx[ip].sad.get_volume() / itin.nat
            xe = stepx[ip].sad.get_e() - reace
            sdata.G.add_node(xnode_name, energy=xe, volume=xvol)
            ynode_name = 'ys' + str(yid)
            yvol = stepy[ip].sad.get_volume() / itin.nat
            ye = stepy[ip].sad.get_e() - reace
            sdata.G.add_node(ynode_name, energy=ye, volume=yvol)
            pxnode_name = 'xl' + str(pstepx[ip].loc.get_iden())
            pynode_name = 'yl' + str(pstepy[ip].loc.get_iden())
            sdata.G.add_edge(pxnode_name, xnode_name)
            sdata.G.add_edge(pynode_name, ynode_name)



        del(xsets)
        del(ysets)

        jobids = pushjob(xkeep, ykeep)
        if checkjob(jobids) == 0:
            (xsets, ysets) = pulljob(xkeep, ykeep)
        else:
            return 100

        for ip in range(itin.npop):
            stepx[ip].loc = cp(xsets[ip])
            stepy[ip].loc = cp(ysets[ip])
            xid = update_iden(sdata.xllist, stepx[ip].loc)
            yid = update_iden(sdata.yllist, stepy[ip].loc)
            stepx[ip].loc.set_iden(xid)
            stepy[ip].loc.set_iden(yid)
            stepx[ip].loc.set_sm('M')
            stepy[ip].loc.set_sm('M')
            stepx[ip].loc.add_left(stepx[ip].sad.get_iden())
            stepy[ip].loc.add_left(stepy[ip].sad.get_iden())
            stepx[ip].sad.add_right(xid)
            stepy[ip].sad.add_right(yid)

            print "ZLOG: STEP %4d, IP %4d X SAD EN: %8.7E, X LOC EN: %8.7E" %\
                  (istep, ip, stepx[ip].sad.get_e(), stepx[ip].loc.get_e())
            print "ZLOG: STEP %4d, IP %4d Y SAD EN: %8.7E, Y LOC EN: %8.7E" %\
                  (istep, ip, stepy[ip].sad.get_e(), stepy[ip].loc.get_e())
            # Xen.append(xsad.get_e())
            # Yen.append(ysad.get_e())
            sdata.xllist.append(stepx[ip].loc)
            sdata.yllist.append(stepy[ip].loc)
            xnode_name = 'xl' + str(xid)
            xvol = stepx[ip].loc.get_volume() / itin.nat
            xe = stepx[ip].loc.get_e() - reace
            sdata.G.add_node(xnode_name, energy=xe, volume=xvol)
            ynode_name = 'yl' + str(yid)
            yvol = stepy[ip].loc.get_volume() / itin.nat
            ye = stepy[ip].loc.get_e() - reace
            sdata.G.add_node(ynode_name, energy=ye, volume=yvol)
            pxnode_name = 'xs' + str(stepx[ip].sad.get_iden())
            pynode_name = 'ys' + str(stepy[ip].sad.get_iden())
            sdata.G.add_edge(pxnode_name, xnode_name)
            sdata.G.add_edge(pynode_name, ynode_name)


        del(xsets)
        del(ysets)
        dumpdata()

        # update gbest
        xyldist = []
        for ix in range(len(sdata.xllist)):
            fpx = sdata.xllist[ix].get_sfp()
            # ex = sdata.xslist[ix].get_e() - reace
            # ex = get_barrier(sdata.xllist, sdata.xslist, sdata.reactant,
            #                  sdata.xllist[ix])
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_sfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                # ey = sdata.yslist[iy].get_e() - reace
                # ey = get_barrier(sdata.yllist, sdata.yslist, sdata.product,
                #                  sdata.yllist[iy])
                # ee = max(ex, ey)
                ee = 0.0
                if dist < 1e-4:
                    xdist = 1e-4
                else:
                    xdist = dist
                mdist = np.log10(xdist) + ee
                # print 'ZLOG: mdist, dist, log(dist), ee,',\
                #       mdist, dist, np.log(dist), ee, 'G', ix, iy
                xyldist.append([dist, [ix, iy], dist, ee])
                if dist < itin.dist:
                    xnode_name = 'xl' + str(sdata.xllist[ix].get_iden())
                    ynode_name = 'yl' + str(sdata.yllist[iy].get_iden())
                    sdata.G.add_edge(xnode_name, ynode_name)
        xyldistSort = sorted(xyldist, key=lambda x: x[2])
        ix = xyldistSort[0][1][0]
        iy = xyldistSort[0][1][1]
        sdata.gbestx = cp(sdata.xllist[ix])
        sdata.gbesty = cp(sdata.yllist[iy])
        bestdist = xyldistSort[0][2]
        # print 'ZLOG: STEP %4d, bestDist: %8.7E' % (istep, bestdist)
        # print "ZLOG: DEBUG: ix", ix, "iy", iy, "len xs", len(xslist), "len xl", len(xllist),\
        #       "len ys", len(yslist), "len yl", len(yllist)
        print "ZLOG: STEP %4d, bestDist: %8.7E, X-Y: %d %d" %\
              (istep, bestdist, ix, iy)
        # print "ZLOG: X %d SAD-E: %8.7E LOC-E: %8.7E" % \
        #       (ix, sdata.xslist[ix].get_e(), sdata.xllist[ix].get_e())
        # print "ZLOG: Y %d SAD-E: %8.7E LOC-E: %8.7E" % \
        #       (iy, sdata.yslist[iy].get_e(), sdata.yllist[iy].get_e())

        write_de(xyldist, sdata.reactant.get_e())

        # updata pbest
        xydist = []
        for ix in range(itin.npop):
            fpx = stepx[ix].loc.get_sfp()
            # ex = Xsad[ix].get_e() - reace
            # ex = get_barrier(sdata.xllist, sdata.xslist, sdata.reactant,
                             # stepx[ix].loc)
            xytdist = []
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_sfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                # ey = sdata.yslist[iy].get_e() - reace
                # ey = get_barrier(sdata.yllist, sdata.yslist, sdata.product,
                                 # sdata.yllist[iy])
                # ee = max(ex, ey)
                ee = 0
                if dist < 1e-4:
                    xdist = 1e-4
                else:
                    xdist = dist
                # mdist = np.log10(xdist) + ee
                # print 'ZLOG: mdist, dist, log(dist), ee',\
                      # mdist, dist, np.log(dist), ee, 'X', ix, iy, 'IP', istep
                xytdist.append([dist, iy])
            xytdistSort = sorted(xytdist, key=lambda x: x[0])
            xytbestdist = xytdistSort[0][0]
            # iy = xytbestdist[0][1]
            # get the best dist for each particle in the pop
            xydist.append(xytbestdist)
            if xytbestdist < sdata.pdistx[ix]:
                sdata.pbestx[ix] = cp(stepx[ix].loc)  # cp(Xloc[ix])
                sdata.pdistx[ix] = xytbestdist

        yxdist = []
        for iy in range(itin.npop):
            fpy = stepy[iy].loc.get_sfp()
            # ey = Ysad[iy].get_e() - reace
            # ey = get_barrier(sdata.yllist, sdata.yslist, sdata.product,
                             # stepy[iy].loc)
            yxtdist = []
            for ix in range(len(sdata.xllist)):
                fpx = sdata.xllist[ix].get_sfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                # ex = sdata.xslist[ix].get_e() - reace
                # ex = get_barrier(sdata.xllist, sdata.xslist, sdata.reactant,
                                 # sdata.xllist[ix])
                # ee = max(ex, ey)
                ee = 0
                if dist < 1e-4:
                    xdist = 1e-4
                else:
                    xdist = dist
                # mdist = np.log10(xdist) + ee
                # print 'ZLOG: mdist, dist, log(dist), ee',\
                      # mdist, dist, np.log(dist), ee, 'Y', ix, iy, 'IP', istep
                yxtdist.append([dist, ix])
            yxtdistSort = sorted(yxtdist, key=lambda x: x[0])
            yxtbestdist = yxtdistSort[0][0]
            yxdist.append(yxtbestdist)
            if yxtbestdist < sdata.pdisty[iy]:
                sdata.pbesty[iy] = cp(stepy[iy].loc)  # cp(Yloc[iy])
                sdata.pdisty[iy] = yxtbestdist

        # sdata.xlocs.append(Xloc)
        # sdata.ylocs.append(Yloc)

        # if abs(bestdist) < itin.dist:
        #     print "ZLOG: CONVERGED!"
        #     break
        dumpdata()

        dataf = 'stepx_' + str(istep) + '.bin'
        f = open(dataf, 'w')
        pick.dump(stepx, f)
        f.close()
        dataf = 'stepy_' + str(istep) + '.bin'
        f = open(dataf, 'w')
        pick.dump(stepy, f)
        f.close()
        del(stepx)
        del(stepy)
        del(pstepx)
        del(pstepy)

        showpath()


def showpath():
    print 'ZLOG: start showpath'
    paths = nx.all_simple_paths(sdata.G, source='xl0', target='yl0', cutoff=15)
    pathdata = []
    for path in paths:
        # print 'ZLOG: PATH:', path
        ee = []
        for node in path:
            e = sdata.G.node[node]['energy']
            ee.append(e)
        be = max(ee)
        pathdata.append([be, path])
        # print 'ZLOG: NODE-E:', ee
        # print 'ZLOG: BARRIER:', be
    if len(pathdata) > 0:
        sorpathdata = sorted(pathdata, key = lambda x: x[0])
        print 'ZLOG: BARRIER: ', sorpathdata[0][0]
        print 'ZLOG: PATH:', sorpathdata[0][1]
    print 'ZLOG: end showpath'

def pushjob(xkeep, ykeep):
    jobids = []
    if itin.client == 'pbs':
        cdirs = sdata.xdirs + sdata.ydirs
        cwdn = os.getcwd()
        for cdir in cdirs:
            ddir = cwdn + '/' + cdir
            f = open('sub.sh', 'w')
            f.write("ssh memex.local << !\n")
            f.write("cd " + ddir + "\n")
            f.write("qsub pbs.sh\n")
            f.write("!\n")
            f.close()
            jbuff = os.popen('sh sub.sh').read()
            # this is desinged for memex cluster
            jid = jbuff.strip()
            jobids.append(jid)
    elif itin.client == 'local':
        # cdirs = sdata.xdirs + sdata.ydirs
        cdirs = []
        for ip in range(itin.npop):
            if xkeep[ip] == 0:
                cdirs.append(sdata.xdirs[ip])
        for ip in range(itin.npop):
            if ykeep[ip] == 0:
                cdirs.append(sdata.ydirs[ip])
        print 'ZLOG: cal dir', cdirs
        print 'ZLOG: len dir', len(cdirs)

        for cdir in cdirs:
            print 'ZLOG: START JOB in dir:', cdir
            os.system('cd ' + cdir + '; sh pbs.sh')
    else:
        print 'ZLOG: ERROR client'
        sys.exe(0)
    print 'ZLOG: jobids', jobids
    return jobids


def pulljob(xkeep, ykeep):
    xsets = []
    ysets = []
    for ip in range(itin.npop):
        if xkeep[ip] == 0:
            xdir = sdata.xdirs[ip]
            try:
                f = open(xdir + '/pcell.bin')
                xx = pick.load(f)
                f.close()
            except:
                print 'ZLOG: fail to pull x pcell.bin'
                xx = set_cell_from_vasp(xdir + '/POSCAR.F')
                xx.set_e(31118.)
        else:
            xx = 'null'

        if ykeep[ip] == 0:
            ydir = sdata.ydirs[ip]
            try:
                f = open(ydir + '/pcell.bin')
                yy = pick.load(f)
                f.close()
            except:
                print 'ZLOG: fail to pull y pcell.bin'
                yy = set_cell_from_vasp(ydir + '/POSCAR.F')
                yy.set_e(31118.)
        else:
            yy = 'null'

        xsets.append(xx)
        ysets.append(yy)
    return (xsets, ysets)


def update_iden(xlist, cell):
    fpc = cell.get_sfp()
    ec = cell.get_e()
    oldids = []
    if len(xlist) == 0:
        return 0
    for x in xlist:
        fpx = x.get_sfp()
        idx = x.get_iden()
        (d, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpc)
        ex = x.get_e()
        edf = np.abs(ec - ex)
        if d < itin.dist and edf < itin.ediff:
            idc = idx
            return idc
        oldids.append(idx)
    idc = max(oldids) + 1
    return idc


def gen_pmode(xy, xcell, v0, istep, ip):
    if xy is 'x':
        pbest = cp(sdata.pbestx[ip])
        gbest = cp(sdata.gbestx)
        # v0 = cp(sdata.evox[istep - 1][ip].v)
    elif xy is 'y':
        pbest = cp(sdata.pbesty[ip])
        gbest = cp(sdata.gbesty)
        # v0 = cp(sdata.evoy[istep - 1][ip].v)
    else:
        print 'ZOUT ERROR xy'

    c1 = 2.0
    c2 = 2.0
    w = 0.9 - 0.5 * (istep) / itin.instep
    (r1, r2) = np.random.rand(2)
    # print pbest.get_lattice()
    # print xcell.get_lattice()
    v = v0 * w + c1 * r1 * getx(pbest, xcell) + c2 * r2 * getx(gbest, xcell)
    return v


def gen_psaddle(xy, xcell, istep, ip):
    if xy is 'x':
        pbest = cp(sdata.pbestx[ip])
        gbest = cp(sdata.gbestx)
        v0 = sdata.vx[ip]
    elif xy is 'y':
        pbest = cp(sdata.pbesty[ip])
        gbest = cp(sdata.gbesty)
        v0 = sdata.vy[ip]
    else:
        print 'ERROR: xy'
        sys.exit(0)
    c1 = 2.0
    c2 = 2.0
    w = 0.9 - 0.5 * (istep + 1) / itin.instep
    (r1, r2) = np.random.rand(2)
    v = v0 * w + c1 * r1 * getx(pbest, xcell) + c2 * r2 * getx(gbest, xcell)
    try:
        scell = rundim(xcell, v)
    except:
        scell = cp(xcell)
        v = cp(v0)
    return (scell, v)


def connect_path(ine, mlisted, slisted, xm, xend, fatherids, xpath):
    # input: saddlelist, minimalist, npop, istep
    # reactant/product minimalist[0]
    # xend : the end point, either product or reactant
    # find the first neighbor saddle for xend

    # MERGE mlist and slist
    # mlisted = mergelr(mlist)
    # slisted = mergelr(mlist)

    # snode0 = []
    # snode0id = []
    # for xs in slist:
    #    if 0 in xs.get_nbor():
    #        snode0.append(xs)
    #        snode0id.append(xs.get_iden())

    # print 'm left', xm.get_left()

    for sp_id in xm.get_left():
        if sp_id not in fatherids and sp_id > -1:
            fatherids.append(sp_id)
            sp = getx_fromid(sp_id, slisted)
            e = sp.get_e() - ine
            xpath.create_node('Saddle' + str(sp.get_iden()) + 'E' + str(e),
                              sp.get_nid(), parent=xm.get_nid(), data=sp)
            for m_id in sp.get_left():
                if m_id == 0:
                    # connect the xend
                    xend.set_nid(xend.get_nid() - 1)
                    xpath.create_node('END', xend.get_nid(),
                                      parent=sp.get_nid(), data=xend)
                    # print '# ZLOG: CONNECTED MID', m_id
                else:
                    # print '# ZLOG: SON ID', m_id
                    mp = getx_fromid(m_id, mlisted)
                    e = mp.get_e() - ine
                    xpath.create_node('Minima' + str(mp.get_iden()) + 'E' +
                                      str(e), mp.get_nid(),
                                      parent=sp.get_nid(), data=mp)
                    connect_path(ine, mlisted, slisted, mp,
                                 xend, fatherids, xpath)
    return 0


def getx_fromid(xid, listed):
    for xterm in listed:
        if xterm.get_iden() == xid:
            sdata.nidp += 1
            xterm.set_nid(sdata.nidp)
            return cp(xterm)
    print 'ERROR: getx_fromid', xid
    exit(1)


def mergelist(xlist):
    # xlist = []
    # for xx in xlist0:
    #     if xx.get_e() < 1000:
    #         xlist.append(xx)
    xlisted = []
    xid = []
    for xc in xlist:
        xid.append(xc.get_iden())

    # print 'set(xid)', set(xid)
    for idt in set(xid):
        simit = []
        lt = []
        rt = []
        for xc in xlist:
            # print xc.get_iden()
            if xc.get_iden() == idt:
                simit.append(xc)
                lt += xc.get_left()
                rt += xc.get_right()
                xt = cp(xc)
        ltt = list(set(lt))
        rtt = list(set(rt))
        xt.set_left(ltt)
        xt.set_right(rtt)
        xlisted.append(xt)
    return xlisted


def write_de(xyldist, e0):
    k = 0
    # for ix in range(len(sdata.xllist)):
    #     for iy in range(len(sdata.yllist)):
    #         e1 = sdata.xslist[ix].get
    #         de = min()
    data = []
    for xs in sdata.xslist:
        for ys in sdata.yslist:
            ex = xs.get_e() - e0
            ey = ys.get_e() - e0
            ee = max(ex, ey)
            md = xyldist[k][0]
            d = xyldist[k][2]
            k += 1
            data.append([d, md, ee])

    f = open('de.dat', 'w')
    for x in data:
        f.write("%8.7E  %8.7E  %8.7E\n" % tuple(x))
    f.close()


def get_barrier(mlist, slist, startp, endp):
    # get barrier energy for startp (reactant/product)
    # to endp (one local minima)
    mlisted = mergelist(mlist)
    slisted = mergelist(slist)
    endp.set_nid(-1)
    ine = endp.get_e()
    if endp.get_iden() > 0:
        xpath = Tree()
        fatherids = []
        endp.set_nid(0)
        xpath.create_node('M' + str(endp.get_iden()), 0, data=endp)
        sdata.nidp = 0
        connect_path(ine, mlisted, slisted, endp, startp, fatherids, xpath)

        dd = []
        for xx in xpath.all_nodes():
            if xx.identifier < 0:
                d = []
                d.append(xx.data)
                nid = xx.bpointer
                nxx = xpath.get_node(nid)
                while True:
                    d.append(nxx.data)
                    if nxx.is_root():
                        break
                    nid = nxx.bpointer
                    nxx = xpath.get_node(nid)
                dd.append(d)
        xe = []
        for x in dd:
            ee = []
            for dx in x:
                exx = dx.get_e() - sdata.reace
                ee.append(exx)
            xe.append(max(ee))
        mxe = min(xe)
    else:
        mxe = 0.0

    return mxe


def main():
    (reac, prod) = initrun()
    sdata.reactant = reac
    sdata.product = prod
    w20init()
    woo()
    f = open('xm.dat', 'w')
    pick.dump(sdata.xllist, f)
    f.close()
    f = open('xs.dat', 'w')
    pick.dump(sdata.xslist, f)
    f.close()
    f = open('ym.dat', 'w')
    pick.dump(sdata.yllist, f)
    f.close()
    f = open('ys.dat', 'w')
    pick.dump(sdata.yslist, f)
    f.close()
    # outputw()


def utest1():
    f = open('xm.dat')
    xmlist = pick.load(f)
    f.close()
    f = open('xs.dat')
    xslist = pick.load(f)
    f.close()
    xmlisted = mergelist(xmlist)
    xslisted = mergelist(xslist)
    print 'nxmlist, nxmlisted', len(xmlist), len(xmlisted)
    print 'nxslist, nsmlisted', len(xslist), len(xslisted)
    xend = cp(xmlist[0])
    xm = cp(xmlist[-3])
    print xm.get_iden()
    fatherids = []
    xpath = Tree()
    xm.set_nid(0)
    xpath.create_node('Minima' + str(xm.get_iden()), 0, data=xm)
    xpath.show()
    sdata.nidp = 0
    xend.set_nid(-1)
    connect_path(xmlisted, xslisted, xm, xend, fatherids, xpath)
    xpath.show()

    dd = []
    for xxx in xpath.all_nodes():
        if xxx.identifier < 0:
            print 'WA'
            d = []
            d.append(xxx.data)
            nid = xxx.bpointer
            xx = xpath.get_node(nid)
            while True:
                d.append(xx.data)
                if xx.is_root():
                    break
                nid = xx.bpointer
                xx = xpath.get_node(nid)
            dd.append(d)

    print len(dd)

    for x in dd:
        print 'EE',
        ee = []
        for xx in x:
            print xx.get_e() - xend.get_e(),
            ee.append(xx.get_e() - xend.get_e())
        print
        print max(ee)


def utest2():
    f = open('xllist.bin')
    xmlist = pick.load(f)
    f.close()
    f = open('xslist.bin')
    xslist = pick.load(f)
    f.close()
    f = open('yllist.bin')
    ymlist = pick.load(f)
    f.close()
    f = open('yslist.bin')
    yslist = pick.load(f)
    f.close()

    xmlisted = mergelist(xmlist)
    xslisted = mergelist(xslist)
    ymlisted = mergelist(ymlist)
    yslisted = mergelist(yslist)

    xend = cp(xmlist[0])
    # print 'xend', xend.get_e()
    yend = cp(ymlist[0])

    sdata.types = xend.get_types()

    goodlist = []
    for xx in xmlist:
        fpx = xx.get_sfp()
        for yy in ymlist:
            fpy = yy.get_sfp()
            (d, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
            if d < itin.dist:
                print 'dd', d
                goodlist.append([xx, yy])
                print xx.get_iden(), yy.get_iden()

    print 'len', len(goodlist)

    ine = xend.get_e()

    xend.set_nid(-1)
    yend.set_nid(-1)
    kk = 0
    mxy = []
    for xyxy in goodlist:
        print 'kk', kk
        kk += 1
        xx = xyxy[0]
        yy = xyxy[1]
        print 'xx, yy', xx.get_iden(), yy.get_iden()
        if xx.get_iden() > 0:
            xpath = Tree()
            fatherids = []
            xx.set_nid(0)
            xpath.create_node('Minima' + str(xx.get_iden()), 0, data=xx)
            sdata.nidp = 0
            connect_path(ine, xmlisted, xslisted, xx, xend, fatherids, xpath)

        if yy.get_iden() > 0:
            ypath = Tree()
            fatherids = []
            yy.set_nid(0)
            ypath.create_node('Minima' + str(yy.get_iden()), 0, data=yy)
            sdata.nidp = 0
            connect_path(ine, ymlisted, yslisted, yy, yend, fatherids, ypath)

        if xx.get_iden() > 0:
            dd = []
            for xxx in xpath.all_nodes():
                if xxx.identifier < 0:
                    d = []
                    d.append(xxx.data)
                    nid = xxx.bpointer
                    nxx = xpath.get_node(nid)
                    while True:
                        d.append(nxx.data)
                        if nxx.is_root():
                            break
                        nid = nxx.bpointer
                        nxx = xpath.get_node(nid)
                    dd.append(d)
            xe = []
            for x in dd:
                ee = []
                for dx in x:
                    exx = dx.get_e() - xend.get_e()
                    ee.append(exx)
                xe.append(max(ee))
            mxe = min(xe)
        else:
            mxe = 0.0

        if yy.get_iden() > 0:
            dd = []
            for yyy in ypath.all_nodes():
                if yyy.identifier < 0:
                    d = []
                    d.append(yyy.data)
                    nid = yyy.bpointer
                    nyy = ypath.get_node(nid)
                    while True:
                        d.append(nyy.data)
                        if nyy.is_root():
                            break
                        nid = nyy.bpointer
                        nyy = ypath.get_node(nid)
                    dd.append(d)
            ye = []
            for y in dd:
                ee = []
                for dy in y:
                    exx = dy.get_e() - xend.get_e()
                    ee.append(exx)
                ye.append(max(ee))
            mye = min(ye)
        else:
            mye = 0.0

        print 'BAR', max(mxe, mye), mxe, mye
        mxy.append([max(mxe, mye), xx, yy])

    mxysort = sorted(mxy, key=lambda x: x[0])
    print mxysort[0][0]
    print mxysort[0][1].get_iden()
    print mxysort[0][2].get_iden()

    goodx = cp(mxysort[0][1])
    idenx = goodx.get_iden()
    fatherids = []
    xpath = Tree()
    goodx.set_nid(0)
    xpath.create_node('Minima' + str(goodx.get_iden()), 0, data=goodx)
    xpath.show()
    sdata.nidp = 0
    xend.set_nid(-1)
    connect_path(ine, xmlisted, xslisted, goodx, xend, fatherids, xpath)
    xpath.show()
    xpath.save2file('xp.dat')


if __name__ == "__main__":
    main()
    # utest2()
