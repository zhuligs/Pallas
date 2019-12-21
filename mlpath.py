#!/usr/bin/env python 

import numpy as np
import pickle as pick
from sklearn import svm
import fppy

import itin

class Mlpath(object):
    def __init__(self):
        self.num_minima = 0
        self.num_saddle = 0
        self.xl = None
        self.yl = None
        self.xs = None
        self.ys = None
        self.d_e = None
        self.clf = None
    
    def loadbin(self):
        f = open('xllist.bin')
        self.xl = pick.load(f)
        f.close()
        f = open('yllist.bin')
        self.yl = pick.load(f)
        f.close()
        f = open('xslist.bin')
        self.xs = pick.load(f)
        f.close()
        f = open('yslist.bin')
        self.ys = pick.load(f)
        f.close()
        self.types = self.xl[0].get_types()

    def get_de(self):
        xx = []
        yy = []
        numx = len(self.xl)
        for i in range(numx-1):
            ei = self.xl[i].get_e()
            fpi = self.xl[i].get_lfp()
            for j in range(i+1, numx):
                ej = self.xl[j].get_e()
                fpj = self.xl[j].get_lfp()
                (dij, m) = fppy.fp_dist(itin.ntyp, self.types, fpi, fpj)
                de = abs(ei - ej)
                xx.append(dij)
                yy.append(de)
        
        numy = len(self.yl)
        for i in range(numy-1):
            ei = self.yl[i].get_e()
            fpi = self.yl[i].get_lfp()
            for j in range(i+1, numy):
                ej = self.yl[j].get_e()
                fpj = self.yl[j].get_lfp()
                (dij, m) = fppy.fp_dist(itin.ntyp, self.types, fpi, fpj)
                de = abs(ei - ej)
                xx.append(dij)
                yy.append(de)
        return (xx, yy)

    def training(self):
        xx = []
        yy = []
        numx = len(self.xl)
        for i in range(numx-1):
            ei = self.xl[i].get_e()
            fpi = self.xl[i].get_lfp()
            for j in range(i+1, numx):
                ej = self.xl[j].get_e()
                fpj = self.xl[j].get_lfp()
                (dij, m) = fppy.fp_dist(itin.ntyp, self.types, fpi, fpj)
                de = abs(ei - ej)
                xx.append(dij)
                yy.append(de)
        
        numy = len(self.yl)
        for i in range(numy-1):
            ei = self.yl[i].get_e()
            fpi = self.yl[i].get_lfp()
            for j in range(i+1, numy):
                ej = self.yl[j].get_e()
                fpj = self.yl[j].get_lfp()
                (dij, m) = fppy.fp_dist(itin.ntyp, self.types, fpi, fpj)
                de = abs(ei - ej)
                xx.append(dij)
                yy.append(de)           

        self.clf = svm.SVR(gamma='scale',tol=0.0001, verbose=True)
        self.clf.fit(xx, yy)

    def predict(self):
        pass




    def run(self):
        self.loadbin(self)
        self.training(self)
        self.predict(self)

def ptest():
    ml = Mlpath()
    ml.loadbin()
    (xx, yy) = ml.get_de()
    f = open('d_e.txt')
    for i in range(len(xx)):
        f.write("%g  %g\n" % (xx[i], yy[i]))
    f.close()

def main():
    mlpath = Mlpath()
    mlpath.run()

if __name__ == '__main__':
    ptest()