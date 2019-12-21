#!/usr/bin/env python 

import numpy as np
import glob
from zfunc import set_cell_from_vasp, write_cell_to_vasp
import fppy
import itin
import time
import joblib
from sklearn import manifold
from sklearn.metrics import euclidean_distances

class Pallasmap(object):
    def __init__(self, *args, **kwargs):
        self.nstruct = 0
        self.structlist = []
        self.natlist = []
        self.maxnat = 0
        self.delist = []
        # return super().__init__(*args, **kwargs)    

    def readfiles(self, dirs):
        for d in dirs:
            fenergy = d + '/energy.dat'
            buff = []
            with open(fenergy) as f:
                for line in f:
                    buff.append(line.split())
            # ens = np.array(buff, float)
            ens = np.array([x[0] for x in buff], float)
            structs = glob.glob(d + '/*.vasp')
            structs = sorted(structs, key=lambda x:int(x.split('_')[1]))
            print len(ens)
            print len(structs)
            for i in range(len(ens)):
                xcell = set_cell_from_vasp(structs[i])
                xcell.set_e(ens[i])
                xcell.set_name(d + '.' + str(i+1) + '_' + buff[i][1])
                self.structlist.append(xcell)
                self.natlist.append(xcell.get_nat())
        self.nstruct = len(self.structlist)
        self.maxnat = max(self.natlist)
        
    def refp(self, fp, maxnat):
        if maxnat == len(fp):
            return fp
        else:
            nd = maxnat/len(fp)
            refp = []
            # retypes = []
            # retypt = np.array(typt, int) * nd
            # for i in range(len(retypt)):
            #     retypes += [i + 1] * retypt[i]
            # retypes = np.array(retypes, int)
            j = 0
            for i in range(len(fp)):
                for j in range(nd):
                    refp.append(fp[i])
            refp = np.array(refp, float)
            return refp
    
    def remove_duplicates(self):
        tol = 0.001
        newlist = []
        declist = [1] * self.nstruct
        for i in range(self.nstruct - 1):
            for j in range(i+1, self.nstruct):
                if self.dijmatrix[i][j] < tol:
                    declist[j] = 0
        self.nstructp = sum(declist)
        self.dijmatrixp = np.array(np.zeros((self.nstructp, self.nstructp)))
        oldindex = []
        self.penergy=[]
        for i in range(self.nstruct):
            if declist[i] == 1:
                newlist.append(self.structlist[i])
                self.penergy.append(self.structlist[i].e)
                oldindex.append(i)
        for i in range(self.nstructp - 1):
            for j in range(i+1, self.nstructp):
                self.dijmatrixp[i][j] = self.dijmatrix[oldindex[i]][oldindex[j]]
                self.dijmatrixp[j][i] = self.dijmatrixp[i][j] 
        self.pstructures = newlist
    
    def get_dijmatrix(self):
        self.dijmatrix = np.array(np.zeros((self.nstruct, self.nstruct)))
        for i in range(self.nstruct - 1):
            fpi = self.structlist[i].get_lfp()
            nati = self.structlist[i].get_nat()
            namei = self.structlist[i].name
            # typesi = self.structlist[i].get_types()
            # typti = self.structlist[i].get_typt()
            # refpi = self.refp(fpi)   
            for j in range(i+1, self.nstruct):
                fpj = self.structlist[j].get_lfp()
                natj = self.structlist[j].get_nat()
                namej = self.structlist[j].name
                # typesj = self.structlist[j].get_types()
                # typtj = self.structlist[j].get_typt()
                maxnat = lcm(nati, natj)
                refpi = self.refp(fpi, maxnat)
                refpj = self.refp(fpj, maxnat)
                # for k in range(len(retypesi)):
                #     retypesi[k] = 1
                types = np.array([1]*maxnat, int)
                # print i, j
                (d, m) = fppy.fp_dist(1, types, refpi, refpj)
                ediff = abs(self.structlist[i].e - self.structlist[j].e)
                self.dijmatrix[i][j] = d
                self.dijmatrix[j][i] = d
                self.delist.append((d, ediff, i, j, namei, namej))
            print ('%s finish %d of %d' % (time.ctime(), i, self.nstruct-1))
    
    def pmds(self):
        mds = manifold.MDS(n_components=2, max_iter=5000, eps=1e-3,
                   dissimilarity="precomputed", n_jobs=1, verbose=3)
        print 'break0'
        pos = self.dijmatrixp.copy()
        print 'break1'
        npos = mds.fit_transform(pos)
        print 'break3'
        self.mdspos = npos
        print 'break4'


def gcd(a,b):
    """Compute the greatest common divisor of a and b"""
    while b > 0:
        a, b = b, a % b
    return a
    
def lcm(a, b):
    """Compute the lowest common multiple of a and b"""
    return a * b / gcd(a, b)


if __name__ == '__main__':
    pmap = Pallasmap()
    pmap.readfiles(['SR'])
    print ('nstruct, ', pmap.nstruct)
    pmap.get_dijmatrix()
    joblib.dump(pmap, 'pmap.bin')
    with open('pmapdij.txt', 'w') as f:
        for x in pmap.delist:
            f.write("%15.9f %15.9f %4d %4d %s %s\n" % tuple(x))
    
    # pmap = joblib.load('pmap.bin')
    pmap.remove_duplicates()
    print (pmap.nstructp)
    pmap.pmds()
    f = open('mdsmap.txt', 'w')
    for i in range(len(pmap.mdspos)):
        x = pmap.mdspos[i]
        name = pmap.pstructures[i].name
        f.write("%15.9f  %15.9f  %15.9f  %s\n" % (x[0], x[1], pmap.penergy[i], name))
    f.close()








            
            
            
