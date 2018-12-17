#!/usr/bin/env python

## ./post-proc.py post-proc/data_
## ./post-proc.py post-proc/data_ > data_x

import numpy as np
from scipy.interpolate import interp1d
import glob, sys

#fnametpl = "post-proc/data_"

if len(sys.argv)==1:
    print "missing file name pattern"
    sys.exit(0)

fnametpl = "%s*" % sys.argv[1]
datas =  glob.glob(fnametpl)
curve = []
for f in datas:
    with open(f) as ff:
        heat = float(ff.readline())
        head = ff.readline()
    inst,prob,energie,expansion = np.loadtxt(f,usecols=[0,1,4,5], skiprows=2,unpack=True)
    t0 = np.interp(.5, prob, inst)
    e0 = np.interp(t0, inst, energie)
    x0 = np.interp(t0, inst, expansion)
    curve.append((heat, t0, e0, x0))
curve.sort() # sort by 0-column

#from tabulate import tabulate
#print tabulate(curve, tablefmt='plain',floatfmt='.14g', headers=['HEAT','TIME','ENERGIE','EXPANSION'])
#sys.exit(0)



print "HEAT   TIME              ENERGIE      EXPANSION"
for cc in curve:
    for rr in cc:
        print rr,
    print
