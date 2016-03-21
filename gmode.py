import numpy as np
import fppy
import itin

from tsase.neb.util import vunit, vrand

# TODO: re-implenment vunit vrand


def get_mode(cell1, cell2):
    types = cell1.get_types()
    fp1 = np.array(cell1.get_lfp())
    fp2 = np.array(cell2.get_lfp())
    (d, m) = fppy.fp_dist(itin.ntyp, types, fp1, fp2)
    mode = np.zeros((itin.nat + 3, 3))
    mode[-3:] = cell2.get_lattice() - cell1.get_lattice()
    ilat = np.linalg.inv(cell1.get_lattice())
    vol = cell1.get_volume()
    jacob = (vol/itin.nat)**(1.0/3.0) * itin.nat**0.5
    mode[-3:] = np.dot(ilat, mode[-3:]) * jacob
    pos1 = cell1.get_cart_positions()
    pos2 = cell2.get_cart_positions()
    for i in range(itin.nat):
        mode[i] = pos2[m[i]] - pos1[i]
    # mode[:-3] = cell2.get_cart_positions() - cell1.get_cart_positions()
    mode = vunit(mode)
    return mode


def get_rmode():
    mode = np.zeros((itin.nat + 3, 3))
    mode = vrand(mode)
    mode = vunit(mode)
    return mode


def get_0mode():
    return np.zeros((itin.nat + 3, 3))
