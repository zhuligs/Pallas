#!/usr/bin/env python
# encoding: utf-8

import numpy as np

saddlehistory = []
localhistory = []

types = np.array([0])

wdis=[]
wmodes = []
wsads = []
wlocs = []
wvs = []
pmodes = []
gmodes = []
bestdist = 1.0
bestmdist = 100.
gbests = []
leader = None
pbests = []
fitpbest = []
x0s = []

ifpso = []


xslist = []
yslist = []
xllist = []
yllist = []

gbestx = None
gbesty = None
pbestx = []
pbesty = []
pdistx = []
pdisty = []
ifpsox = []
ifpsoy = []
xlocs = []
ylocs = []

vx = []
vy = []

gdir = ''
ddir = ''

xmdb = []  # minima database from x side
xsdb = []  # saddle points database from x side
ymdb = []  # minima database from y side
ysdb = []  # saddle points database from y side

nidp = 0
