#!/usr/bin/python
import numpy as np
import math
import fppy


class Cell:

    def __init__(self,
                 name=None,
                 symbols=None,
                 typt=None,
                 types=None,
                 lattice=None,
                 positions=None,
                 cart_positions=None,
                 znucl=None,
                 latv=None,
                 atomv=None,
                 stress=None,
                 e=None,
                 sfp=None,
                 lfp=None):

        if name is None:
            self.name = None
        else:
            self.name = str(name)

        if typt is None:
            self.typt = None
        else:
            self.typt = np.array(typt)

        if lattice is None:
            self.lattice = None
        else:
            self.lattice = np.array(lattice)

        if positions is None:
            self.positions = None
        else:
            self.positions = np.array(positions)

    def set_name(self, name):
        self.name = str(name)

    def get_name(self):
        return self.name

    def set_lattice(self, lat):
        self.lattice = np.array(lat)

    def get_lattice(self):
        return self.lattice

    def set_znucl(self, znucl):
        self.znucl = np.array(znucl)

    def set_types(self):
        types = []
        for i in range(len(self.typt)):
            types += [i+1] * self.typt[i]
        self.types = np.array(types, int)

    def get_types(self):
        return self.types

    def set_e(self, e):
        self.e = float(e)

    def get_e(self):
        return self.e

    def set_positions(self, pos):
        # for ia in range(len(pos)):
        #     for ib in range(3):
        #         if pos[ia][ib] >= 1: pos[ia][ib] -= int(pos[ia][ib])
        #         if pos[ia][ib] < 0:  pos[ia][ib] -= (int(pos[ia][ib]) - 1)
        self.positions = np.array(pos)
        self.cart_positions = np.dot(pos, self.lattice)

    def set_cart_positions(self, rxyz):
        self.cart_positions = np.array(rxyz)
        self.positions = np.dot(rxyz, np.linalg.inv(self.lattice))

    def get_positions(self):
        return self.positions

    def get_cart_positions(self):
        return self.cart_positions

    def set_typt(self, typ):
        self.typt = np.array(typ)

    def get_typt(self):
        return self.typt

    def set_symbols(self, symb):
        self.symbols = symb

    def get_symbols(self):
        return self.symbols

    def set_latv(self, v):
        self.latv = np.array(v)

    def get_latv(self):
        return self.latv

    def set_atomv(self, v):
        self.atomv = np.array(v)

    def get_atomv(self):
        return self.atomv

    def get_volume(self):
        return np.linalg.det(self.lattice)

    def set_stress(self, stres):
        self.stress = np.array(stres)

    def get_stress(self):
        return self.stress

    def cal_fp(self, cutoff, lmax, natx=300):
        lat = self.lattice
        rxyz = self.get_cart_positions()
        types = self.types
        znucl = self.znucl
        (sfp, lfp) = fppy.fp_periodic(lat, rxyz, types, znucl, lmax, natx,
                                      cutoff)
        self.sfp = sfp
        self.lfp = lfp

    def get_sfp(self):
        return self.sfp

    def get_lfp(self):
        return self.lfp


def lat2vec(lat):
    return np.array([lat[0][0], lat[1][1], lat[2][2],
                     lat[1][0], lat[2][0], lat[2][1]], float)


def vec2lat(vec):
    return np.array([[vec[0], 0., 0.],
                     [vec[3], vec[1], 0.],
                     [vec[4], vec[5], vec[2]]], float)


def lat2lcons(lat):
    ra = math.sqrt(lat[0][0]**2 + lat[0][1]**2 + lat[0][2]**2)
    rb = math.sqrt(lat[1][0]**2 + lat[1][1]**2 + lat[1][2]**2)
    rc = math.sqrt(lat[2][0]**2 + lat[2][1]**2 + lat[2][2]**2)

    cosa = (lat[1][0]*lat[2][0] + lat[1][1]*lat[2][1] +
            lat[1][2]*lat[2][2])/rb/rc
    cosb = (lat[0][0]*lat[2][0] + lat[0][1]*lat[2][1] +
            lat[0][2]*lat[2][2])/ra/rc
    cosc = (lat[0][0]*lat[1][0] + lat[0][1]*lat[1][1] +
            lat[0][2]*lat[1][2])/rb/ra

    alpha = math.acos(cosa)
    beta = math.acos(cosb)
    gamma = math.acos(cosc)

    return np.array([ra, rb, rc, alpha, beta, gamma], float)


def lcons2lat(cons):
    (a, b, c, alpha, beta, gamma) = cons

    bc2 = b**2 + c**2 - 2*b*c*math.cos(alpha)

    h1 = a
    h2 = b * math.cos(gamma)
    h3 = b * math.sin(gamma)
    h4 = c * math.cos(beta)
    h5 = ((h2 - h4)**2 + h3**2 + c**2 - h4**2 - bc2)/(2 * h3)
    h6 = math.sqrt(c**2 - h4**2 - h5**2)

    lattice = [[h1, 0., 0.], [h2, h3, 0.], [h4, h5, h6]]
    return lattice


def get_cutoff(lat):
    volume = np.linalg.det(lat)
    (a, b, c, alpha, beta, gamma) = lat2lcons(lat)
    area_ab = a * b * np.sin(gamma)
    area_ac = a * c * np.sin(beta)
    area_bc = b * c * np.sin(alpha)
    h_ab = volume / area_ab
    h_ac = volume / area_ac
    h_bc = volume / area_bc
    h = np.array([h_ab, h_ac, h_bc], float)
    return h.min() * 0.75 / 2.

