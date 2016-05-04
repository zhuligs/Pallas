#!/usr/bin/env python

import numpy as np
import itin
import sdata

from wrapdimer import get_rmode, get_0mode, get_mode

def con(reac, prod):
	mode = get_mode(reac, prod)
	sdd = runvdim(reac, mode)
	mode = -1.0 * get_mode(sdd, reac)
	newmin = goptv(sdd, mode)
	types = sdata.types()
	fp0 = reac.get_lfp()
