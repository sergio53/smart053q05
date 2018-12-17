#! /usr/bin/python

import  os
import  string
import  sys


def main () :

    file = sys.argv[1]
    fileout = file + '.derive'

    try :
       txt = open(file,'r')
    except :
       sys.exit(0)

    Xold = 0.
    Yold = 0.
    dYdXold = 0.

    out = '#t eps sig depsdt\n'
    txt.readline()
    lines = txt.readlines()

    for i in range(len(lines)):

        try :
                line=lines[i].split()
        except :
                line=''

#        time = float(line[0])
        X = float(line[0])
        Y = float(line[2])
        Z = float(line[1])
        dYdX = dYdXold

        if (X != Xold) and (Y > Yold) : 
           dYdX = (Y-Yold)/(X-Xold)

        out = out + str(X) + ' ' +  str(Y) + ' ' + str(Z) + ' ' +str(dYdX) + '\n'

        Xold = X
        Yold = Y
        dYdXold = dYdX

    output=open(fileout,"w")
    output.write(out)
    output.close() 

main()
