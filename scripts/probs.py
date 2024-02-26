import numpy as np
from astropy.io import ascii

from optparse import OptionParser
parser = OptionParser()

parser.add_option("-f", "--tiles_file",
                  action="store", type="string", dest="tiles_file")


parser.add_option("-d", "--date",
                  action="store", type="string", dest="date")

parser.add_option("-t", "--telescope",
                  action="store", type="string", dest="telescope",help="Options: Swope, Thacher, Nickel")

(options, args) = parser.parse_args()

date=options.date
telescope=options.telescope

if telescope=='Swope':
    schedule=np.loadtxt('LCO_Swope_{0}_GoodSchedule.csv'.format(date),delimiter=',',unpack=True,usecols=0,skiprows=1,
                     dtype=str)
elif telescope=='Thacher':
    schedule=np.loadtxt('Thacher_Thacher_{0}_GoodSchedule.csv'.format(date),delimiter=',',unpack=True,usecols=0,skiprows=1,
                     dtype=str)

elif telescope=='Nickel':
    schedule=np.loadtxt('Lick_Nickel_{0}_GoodSchedule.csv'.format(date),delimiter=',',unpack=True,usecols=0,skiprows=1,
                     dtype=str)

elif telescope=='Keck':
  schedule=np.loadtxt('Keck_Keck_{0}_GoodSchedule.csv'.format(date),delimiter=',',unpack=True,usecols=0,skiprows=1,
                     dtype=str)

scheduled_tiles=[]
for i in schedule:
    if i != "":
        
        scheduled_tiles.append(i)

print ("Number of tiles scheduled: ",len(scheduled_tiles))

tiles_file=ascii.read('{0}'.format(options.tiles_file),data_start=0, delimiter=',')


coincidence=[]
ps=[]
for name in scheduled_tiles:
    for j,val in enumerate(tiles_file['FieldName']):
        if name==val:
            coincidence.append(name)
            ps.append(tiles_file['Priority'][j])

print ("Coincidence", len(coincidence), "(must be equal to number above): it's # overlapping tiles between original tiles file and scheduled tiles)")
print("Sum of Priorities: ", np.sum(ps))
print("Probability: ", 100*np.sum(ps))


