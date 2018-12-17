#!/usr/bin/env python2

sought = 'crack-length*.dat'

import os, shutil, fnmatch

sought_dir = "crack_XX"
if os.path.isdir(sought_dir): shutil.rmtree(sought_dir)
os.mkdir(sought_dir)

for _ in os.listdir('.'):
    if _== sought_dir :
        continue
    calc = "./%s" % _
    if os.path.isdir(calc):
        for data in os.listdir(calc):
            if fnmatch.fnmatch(data, sought):
                src = ".%s/%s" % (calc,data)
                lnk = "./%s/%s" % (sought_dir, data)
                if os.path.islink(lnk): os.unlink(lnk)            
                os.symlink(src,lnk)
                print src, lnk