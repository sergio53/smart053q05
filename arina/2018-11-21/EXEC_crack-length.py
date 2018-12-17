#!/usr/bin/env python2
import os, shutil

sought_dir = "crack_XX"
if os.path.isdir(sought_dir): shutil.rmtree(sought_dir)
os.mkdir(sought_dir)

for _ in os.listdir('.'):
    if not os.path.isdir(_): continue
    if _[:-2]<>'calc_': continue
    medFileName = "./%s/CT_results.rmed" % _ 
    fnout = "./%s/crack-length%s.dat" % (sought_dir, _[-3:])
    if os.path.isfile(medFileName):
        execfile('crack-length.py') 