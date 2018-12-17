#!/usr/bin/env python2

sought = 'REPE_OUT/data'

import os, shutil

sought_dir = "./data_XX"
if os.path.isdir(sought_dir): shutil.rmtree(sought_dir)
os.mkdir(sought_dir)

for _ in os.listdir('.'):
    if os.path.isfile("./%s/%s" % (_, sought)):
        data = "../%s/%s" % (_, sought)
        lnk = "%s/%s%s" % (sought_dir, sought.split('/')[-1], _[-3:])
        if os.path.islink(lnk): os.unlink(lnk)            
        os.symlink(data,lnk)
        print data, lnk