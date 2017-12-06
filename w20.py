#!/usr/bin/python -u
# encoding: utf-8

import sys
import os
import socket
from copy import deepcopy as cp
import numpy as np
import cPickle as pick
import time
# from pele.storage import Database

# from treelib import Node, Tree
import networkx as nx

import itin
import sdata
import fppy
import itdbase
from wrapdimer import get_rmode
from zfunc import rundim, set_cell_from_vasp, write_cell_to_vasp, getx
from w2 import initrun
# from w40 import get_dimfile_ready, get_optfile_ready
# from w40 import checkjob


# sys.setrecursionlimit(100000)


def w20init():
    sdata.servername = itin.servername
    with open('PORT.txt') as f:
        portstring = f.readline()
    sdata.serverport = int(portstring)
    for ip in range(itin.npop):
        xdir = 'Calx' + str(ip)
        ydir = 'Caly' + str(ip)
        os.system('mkdir -p ' + xdir)
        os.system('mkdir -p ' + ydir)
        sdata.xdirs.append(xdir)
        sdata.ydirs.append(ydir)


def create_stepdata():
    stepdata = []
    for i in range(itin.npop):
        stepdata.append(itdbase.Cobj())
    return stepdata


def wo():
    # istep = 0
    reace = sdata.reactant.get_e()
    print 'ZOUT: reace', reace
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
            print 'ZOUT: no keep in xkeep and ykeep'
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
        yid = update_iden(sdata.yslist, stepy[ip].sad)
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
        ye = ysets[ip].get_e() - reace
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
        print "ZOUT: INIT STEP, IP %4d X SAD EN: %8.7E, X LOC EN: %8.7E" %\
              (ip, stepx[ip].sad.get_e(), stepx[ip].loc.get_e())
        print "ZOUT: INIT STEP, IP %4d Y SAD EN: %8.7E, Y LOC EN: %8.7E" %\
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
        ye = ysets[ip].get_e() - reace
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
        ex = stepx[ix].loc.get_e()
        # ex = Xsad[ix].get_e() - reace
        # ex = get_barrier(sdata.xllist, sdata.xslist, sdata.reactant,
        #                  stepx[ix].loc)
        for iy in range(itin.npop):
            fpy = stepy[iy].loc.get_sfp()
            ey = stepy[iy].loc.get_e()
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

            if dist < itin.dist and abs(ex - ey) < itin.ediff:
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
    print "ZOUT: INIT STEP, bestDist: %8.7E, bestmD: %8.7E, X-Y: %4d %4d" %\
          (xydistSort[0][2], xydistSort[0][0],
           xydistSort[0][1][0], xydistSort[0][1][1])

    # update pdist x
    for ix in range(itin.npop):
        xytdist = []
        fpx = stepx[ix].loc.get_sfp()
        for iy in range(len(sdata.yllist)):
            fpy = sdata.yllist[iy].get_sfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
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
        fpy = stepy[iy].loc.get_sfp()
        for ix in range(len(sdata.xllist)):
            fpx = sdata.xllist[ix].get_sfp()
            (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
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
    f = open('G.bin', 'w')
    pick.dump(sdata.G, f)
    f.close()


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
            yid = update_iden(sdata.yslist, stepy[ip].sad)
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

            print "ZOUT: STEP %4d, IP %4d X SAD EN: %8.7E, X LOC EN: %8.7E" %\
                  (istep, ip, stepx[ip].sad.get_e(), stepx[ip].loc.get_e())
            print "ZOUT: STEP %4d, IP %4d Y SAD EN: %8.7E, Y LOC EN: %8.7E" %\
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
        bconnected = False
        xyldist = []
        for ix in range(len(sdata.xllist)):
            fpx = sdata.xllist[ix].get_sfp()
            ex = sdata.xllist[ix].get_e()
            # ex = sdata.xslist[ix].get_e() - reace
            # ex = get_barrier(sdata.xllist, sdata.xslist, sdata.reactant,
            #                  sdata.xllist[ix])
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_sfp()
                ey = sdata.yllist[iy].get_e()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                # ey = sdata.yslist[iy].get_e() - reace
                # ey = get_barrier(sdata.yllist, sdata.yslist, sdata.product,
                #                  sdata.yllist[iy])
                # ee = max(ex, ey)
                ee = 0.0
                # if dist < 1e-4:
                #     xdist = 1e-4
                # else:
                #     xdist = dist
                # mdist = np.log10(xdist) + ee
                # print 'ZLOG: mdist, dist, log(dist), ee,',\
                #       mdist, dist, np.log(dist), ee, 'G', ix, iy
                mdist = 10000.
                if dist < itin.dist and abs(ex - ey) < itin.ediff:
                    bconnected = True
                    xnode_name = 'xl' + str(sdata.xllist[ix].get_iden())
                    ynode_name = 'yl' + str(sdata.yllist[iy].get_iden())
                    sdata.G.add_edge(xnode_name, ynode_name)
                    # if the x - y is connected, set the x y with lowest
                    # barrier energy as the global best
                    xbarrier = get_barrier('x', ix, istep, sdata.xllist[ix], 10)
                    ybarrier = get_barrier('y', iy, istep, sdata.yllist[iy], 10)
                    ebar = max(xbarrier, ybarrier)
                    mdist = ebar
                xyldist.append([dist, [ix, iy], mdist, ee])

        if bconnected:
            xyldistSort = sorted(xyldist, key=lambda x: x[2])
        else:
            xyldistSort = sorted(xyldist, key=lambda x: x[0])
        ix = xyldistSort[0][1][0]
        iy = xyldistSort[0][1][1]
        sdata.gbestx = cp(sdata.xllist[ix])
        sdata.gbesty = cp(sdata.yllist[iy])
        bestmdist = xyldistSort[0][2]
        bfpdist = xyldistSort[0][0]
        print "ZOUT: STEP %4d, fpDist: %8.7E, mDist: %8.7E, X-Y: %d %d" %\
              (istep, bfpdist, bestmdist, ix, iy)
        # print "ZLOG: X %d SAD-E: %8.7E LOC-E: %8.7E" % \
        #       (ix, sdata.xslist[ix].get_e(), sdata.xllist[ix].get_e())
        # print "ZLOG: Y %d SAD-E: %8.7E LOC-E: %8.7E" % \
        #       (iy, sdata.yslist[iy].get_e(), sdata.yllist[iy].get_e())

        write_de(xyldist, sdata.reactant.get_e())

        # updata pbest
        xydist = []
        for ix in range(itin.npop):
            fpx = stepx[ix].loc.get_sfp()
            xytdist = []
            for iy in range(len(sdata.yllist)):
                fpy = sdata.yllist[iy].get_sfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                ee = 0
                # if dist < 1e-4:
                #     xdist = 1e-4
                # else:
                #     xdist = dist
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
            yxtdist = []
            for ix in range(len(sdata.xllist)):
                fpx = sdata.xllist[ix].get_sfp()
                (dist, m) = fppy.fp_dist(itin.ntyp, sdata.types, fpx, fpy)
                ee = 0
                # if dist < 1e-4:
                #     xdist = 1e-4
                # else:
                #     xdist = dist
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
    print 'ZOUT: start showpath'
    paths = nx.all_simple_paths(sdata.G, source='xl0', target='yl0', cutoff=15)
    pathdata = []
    bes = []
    for path in paths:
        # print 'ZLOG: PATH:', path
        ee = []
        for node in path:
            e = sdata.G.node[node]['energy']
            ee.append(e)
        be = max(ee)
        if be < sdata.bestbe:
            pathdata.append([be, path])
            bes.append(be)
        # print 'ZLOG: NODE-E:', ee
        # print 'ZLOG: BARRIER:', be
    if len(pathdata) > 0:
        sorpath = []
        sdata.bestbe = min(bes)
        for pt in pathdata:
            if abs(pt[0] - sdata.bestbe) < 0.0001:
                sorpath.append([len(pt[1]), pt[1]])

        sorpathdata = sorted(sorpath, key=lambda x: x[0])
        sdata.bestpath = sorpathdata[0][1]
        # print 'ZLOG: BARRIER: ', sorpathdata[0][0]
        # print 'ZLOG: PATH:', sorpathdata[0][1]
    print 'ZOUT: BARRIER: ', sdata.bestbe
    print 'ZOUT: PATH:', sdata.bestpath
    print 'ZOUT: end showpath'


def get_barrier(xy, ic, istep, xcell, cutoff):
    reace = sdata.reactant.get_e()
    if xy == 'x':
        sour = 'xl0'
        targ = 'xl' + str(xcell.get_iden())
        paths = nx.all_simple_paths(sdata.G, source=sour, target=targ,
                                    cutoff=cutoff)
        ees = []
        for path in paths:
            for node in path:
                ee = []
                e = sdata.G.node[node]['energy']
                ee.append(e)
            ees.append(max(ee))
        if len(ees) > 0:
            be = min(ees)
        else:
            ee = []
            for i in range(istep):
                ix = ic - i * itin.npop
                ee.append(sdata.xllist[ix].get_e())
                ee.append(sdata.xslist[ix - 1].get_e())
            be = max(ee) - reace
    if xy == 'y':
        sour = 'yl0'
        targ = 'yl' + str(xcell.get_iden())
        paths = nx.all_simple_paths(sdata.G, source=sour, target=targ,
                                    cutoff=cutoff)
        ees = []
        for path in paths:
            for node in path:
                ee = []
                e = sdata.G.node[node]['energy']
                ee.append(e)
            ees.append(max(ee))
        if len(ees) > 0:
            be = min(ees)
        else:
            ee = []
            for i in range(istep):
                ix = ic - i * itin.npop
                ee.append(sdata.yllist[ix].get_e())
                ee.append(sdata.yslist[ix - 1].get_e())
            be = max(ee) - reace
    return be


def pushjob(xkeep, ykeep):
    jobids = []
    if itin.client == 'pbs':
        # cdirs = sdata.xdirs + sdata.ydirs
        cdirs = []
        for ip in range(itin.npop):
            if xkeep[ip] == 0:
                cdirs.append(sdata.xdirs[ip])
        for ip in range(itin.npop):
            if ykeep[ip] == 0:
                cdirs.append(sdata.ydirs[ip])
        print 'ZOUT: cal dir', cdirs
        print 'ZOUT: len dir', len(cdirs)
        cwdn = os.getcwd()
        for cdir in cdirs:
            ddir = cwdn + '/' + cdir
            consock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            consock.connect((sdata.servername, sdata.serverport))
            msg = 'subjob%' + ddir
            consock.send(msg)
            jbuff = consock.recv(2048)
            consock.close()
            # f = open('sub.sh', 'w')
            # f.write("ssh memex.local << !\n")
            # f.write("cd " + ddir + "\n")
            # f.write("qsub pbs.sh\n")
            # f.write("!\n")
            # f.close()
            # jbuff = os.popen('sh sub.sh').read()
            # this is desinged for memex cluster
            jid = jbuff.strip()
            print ('* ZLOG: received job id:', jid)
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


def checkjob(jobids):
    if itin.client == 'local':
        return 0
    finished = False
    while not finished:
        finished = checkids(jobids)
        print time.ctime(), 'ZLOG: JOB FINISHED:', finished
        time.sleep(10)
        if os.path.isfile('CSTOP'):
            os.system('rm -f CSTOP')
            print 'ZLOG: CSTOP'
            allid = ' '.join(jobids)
            consock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            consock.connect((sdata.servername, sdata.serverport))
            msg = 'qdeljob%' + allid
            consock.send(msg)
            jbuff = consock.recv(2048)
            consock.close()
            return 100
    return 0


def checkids(jobids):
    while True:
        consock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        consock.connect((sdata.servername, sdata.serverport))
        msg = 'checkjob%0'
        consock.send(msg)
        print ("* ZLOG: send the checkjob message")
        jbuff = consock.recv(2048)
        consock.close()
        rbuff = jbuff.split('.')
        print ("* ZLOG: received the checkjob message")
        if rbuff[0] == 'DONE':
            return True
        else:
            finished = True
            reid = rbuff[1:]
            for id in jobids:
                if id in reid:
                    print ('id in reid')
                    finished = False
        return finished


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
                xx.set_e(151206.)
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
                yy.set_e(151206.)
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


if __name__ == "__main__":
    main()
    # utest2()

