#!/usr/bin/env python

## ./filter01.py < force_depl > force_depl_2
## ./filter01.py force_depl ==> force_depl.filtered

import sys

if len(sys.argv)==1:
    inp = sys.stdin
    out = sys.stdout
else:    
    inp = open(sys.argv[1],'r')
    out = open(sys.argv[1] + '.filtered','w') 

inp.readline()
inp.readline()
arr = eval(inp.readline())[0]
for r in arr:
    out.write(' '.join(str(e) for e in r) + '\n')
    