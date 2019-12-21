import numpy as np
import networkx as nx
import w20data
import itin

servername = 'memex.locale'
serverport = 18000

G = nx.Graph()
wstore = w20data.Wstorage(ediff=itin.ediff, fpdiff=itin.dist, ntyp=itin.ntyp)

saddlehistory = []
localhistory = []

types = np.array([0])

wdis = []
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
gbest = None
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
reace = 0.0

evodata = []
product = None
reactant = None
xdirs = []
ydirs = []

evox = []
evoy = []

stepx = []
stepy = []
bestbe = 10000.
bestpath = []

xen = []
yen = []
