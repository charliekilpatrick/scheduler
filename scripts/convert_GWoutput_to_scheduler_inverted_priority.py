#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np

from optparse import OptionParser
parser = OptionParser()


parser.add_option("-t", "--tiles_file",
                  action="store", type="string", dest="tiles_file")

parser.add_option("-i", "--input_scheduler",
                  action="store", type="string", dest="input_scheduler")

(options, args) = parser.parse_args()


infile=options.tiles_file

outfile=options.input_scheduler

def main(infile='{0}'.format(infile),
         addstring='2019-04-12,12,GW'):

    #name,ra,dec,p1,p2 = np.loadtxt('tmp.txt',unpack=True,delimiter=',')

    fout = open('{0}'.format(outfile),'w')
    with open(infile) as fin:
        print("Name,RA,DEC,Priority,DiscDates,DiscMags,Type",file=fout)
        for line in fin:
            # print(line)
            line = line.replace('\n','')
            line = line.replace(' ','')
            if line == '' or line.startswith('F') or line.startswith('#'):
                continue
            
            lineparts = line.split(',')
            lineparts[4] = addstring
            # print (1/np.float(lineparts[6])/1000)
            lineparts[3] = np.str(1/np.float(lineparts[6])/1000)
           
            print(",".join(lineparts[:5]),file=fout)

            
    # print ("Tiles file converted into Scheduler input format")
    fout.close()
            
if __name__ == """__main__""":
    main()