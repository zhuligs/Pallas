import numpy as np
from zfunc import set_cell_from_vasp, write_cell_to_vasp
import fppy

x1 = set_cell_from_vasp('P1.vasp')
x2 = set_cell_from_vasp('P2.vasp')
fpx = x1.get_sfp()
fpy = x2.get_sfp()
types = x1.get_types()
print len(fpx)
print types

(dist, m) = fppy.fp_dist(2, types, fpx, fpy)

print dist
print m