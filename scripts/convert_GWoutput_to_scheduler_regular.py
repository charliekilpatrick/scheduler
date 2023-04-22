#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
infile=sys.argv[1]

outfile=sys.argv[2]

def main(infile='{0}'.format(infile),
         addstring='1,2019-04-12,12,GW'):

    #name,ra,dec,p1,p2 = np.loadtxt('tmp.txt',unpack=True,delimiter=',')

    fout = open('{0}.csv'.format(outfile),'w')
    with open(infile) as fin:
        print("Name,RA,DEC,Priority,DiscDates,DiscMags,Type",file=fout)
        for line in fin:
            # print(line)
            line = line.replace('\n','')
            line = line.replace(' ','')
            if line == '' or line.startswith('F') or line.startswith('#'):
                continue
            
            lineparts = line.split(',')
            lineparts[3] = addstring
           
            print(",".join(lineparts[:4]),file=fout)

            
            
    fout.close()
            
if __name__ == """__main__""":
    main()