#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
from astropy.io import ascii

from optparse import OptionParser
parser = OptionParser()

parser.add_option("-s", "--scheduler_file",
                  action="store", type="string", dest="scheduler_file")

parser.add_option("-t", "--tiles_file",
                  action="store", type="string", dest="tiles_file")

parser.add_option("-c", "--fc_file",
                  action="store", type="string", dest="fc_file")

(options, args) = parser.parse_args()

infile=options.scheduler_file

outfile=options.tiles_file

outputfile=options.fc_file



tile_name=np.loadtxt('{0}'.format(infile),delimiter=',',unpack=True,usecols=0,skiprows=1,
                     dtype=str)

# print (len(tile_name))

scheduled_tiles=[]
for i in tile_name:
    if i != "":
        # print (i
        scheduled_tiles.append(i)

print ("Number of tiles scheduled: ",len(scheduled_tiles))

print ("ID of scheduled tiles :", scheduled_tiles)


lines_to_keep = []
with open('{0}'.format(outfile), "r") as f:

     for line in f.readlines():
     	for i in scheduled_tiles:
         if i in line:
             lines_to_keep.append(line)

# Write all the links in our list to the file
with open("{0}".format(outputfile), "w") as f:

    for line in lines_to_keep:
        f.write(line)
f.close()