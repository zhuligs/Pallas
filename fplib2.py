import numpy as np
from munkres import Munkres


def fpdist(ntyp, types, fp1, fp2):
    nat = len(types)
    fpd = 0.0
    for ityp in range(ntyp):
        itype = ityp + 1
        M = []
        for iat in range(nat):
            m = []
            if types[iat] == itype:
                for jat in range(nat):
                    if types[jat] == itype:
                        tfpd = fp1[iat] - fp2[jat]
                        tt = np.sqrt(np.dot(tfpd, tfpd))
                        m.append(tt)
            M.append(m)

        Mu = Munkres()
        indexes = Mu.compute(M)
        total = 0.0
        print len(M), len(M[0])
        print indexes
        for i, j in indexes:
            v = M[i][j]
            total += v

        fpd += total

    fpd = fpd / np.sqrt(nat)
    return fpd
