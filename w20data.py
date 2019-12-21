import fppy
class Wstorage(object):
    def __init__(self, ediff=0.001, fpdiff=0.001, 
                 ntyp=None):
        self.minima = []
        self.saddle = []
        self.ediff = float(ediff)
        self.fpdiff = float(fpdiff)
        self.ntyp = ntyp
        # self.types = types

    def add_minima(self, new):
        for m in self.minima:
            if self.compare(new, m):
                return m
        newid = len(self.minima) + 1
        new.set_iden(newid)
        new.set_name('M' + str(newid))
        self.minima.append(new)
        return new
    
    def add_saddle(self, new, m1, m2):
        id1 = m1.get_iden()
        id2 = m2.get_iden()
        for s in self.saddle:
            if self.compare(new, s):
                s.add_conn(id1)
                s.add_conn(id2)
                retrun s
        newid = len(self.saddle) + 1
        new.set_iden(newid)
        new.add_conn(id1)
        new.add_conn(id2)
        new.set_name('S' + str(newid))
        self.saddle.append(new)
        return new

    def compare(self, cell1, cell2):
        fp1 = cell1.get_lfp()
        fp2 = cell2.get_lfp()
        e1 = cell1.get_e()
        e2 = cell2.get_e()
        types = cell1.get_types()
        (dist, m) = fppy.fp_dist(self.ntyp, types, fp1, fp2)
        if dist < self.fpdiff and abs(e1-e2) < self.ediff:
            return True
        else:
            return False
