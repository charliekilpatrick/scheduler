import os,sys
import numpy as np

import argparse #optparse deprecated


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--tiles_file",
                  action="store", type=str, dest="tiles_file")

# parser.add_option("-i", "--input_scheduler",
#                   action="store", type="string", dest="input_scheduler")

parser.add_argument("-d", "--date",
                  action="store", type=str, dest="date")

# parser.add_option("-o", "--output_scheduler",
#                   action="store", type="string", dest="output_scheduler")

# parser.add_option("-c", "--fc_file",
#                   action="store", type="string", dest="fc_file")

parser.add_argument("-t", "--telescope",
                  action="store", type=str, dest="telescope",help="Options: Swope, Thacher, Nickel")
parser.add_argument("-a", "--now", help="Start Now -- True or False")
parser.add_argument("-b", "--start", help="Desired Start Time in the format of HHMM")
parser.add_argument("-c", "--end", help="Desired End Time in the format of HHMM")
parser.add_argument("-A", "--asap", action='store_true', default=False)



args = parser.parse_args()

tiles_file=args.tiles_file

# input_scheduler=args.input_scheduler
date=args.date
telescope=args.telescope


input_scheduler='{0}_{1}_targets.csv'.format(date,telescope)

# output_scheduler=args.output_scheduler

# fc_file=args.fc_file

date=args.date

# print (telescope)

print("telescope is {0}".format(telescope))

print ("Tiles file converted into Scheduler input format")

os.system('python convert_GWoutput_to_scheduler_inverted_priority.py -t {0} -i {1}'.format(tiles_file,input_scheduler))

print ("***************************")


print ("Running Scheduler")

obs = {'Swope': 'LCO', 'Thacher': 'Thacher', 'Nickel': 'Lick', 'Keck': 'Keck'}


create_schedule_cmd = 'python CreateSchedule.py -d {0} -f {1} --obstele {2}:{3} --now {4} --start {5} --end {6}'.format(date,input_scheduler, obs[telescope], telescope, args.now, args.start, args.end)
if args.asap:
	create_schedule_cmd += " -A"

os.system(create_schedule_cmd)

print ("***************************")

print ("Creating fields-center formatted file")

convert_scheduler_cmd = 'python convert_scheduler_output_to_GWoutput.py -s {0}_{1}_{2}_GoodSchedule.csv -t {3} -c GWSearch_output_{2}_{4}.txt'.format(obs[telescope], telescope, date,tiles_file,telescope)
os.system(convert_scheduler_cmd)
os.system('python probs.py -f {0} --telescope {1} --date {2}'.format(tiles_file, telescope, date))

"""

if telescope=='Swope':

	oscommand = 'python CreateSchedule.py -d {0} -f {1} --obstele LCO:Swope --now {2} --start {3} --end {4}'.format(date,input_scheduler, args.now, args.start, args.end)
	if args.asap:
		oscommand += " -A"
	os.system(oscommand)

	print ("***************************")

	print ("Creating fields-center formatted file")


	os.system('python convert_scheduler_output_to_GWoutput.py -s LCO_Swope_{0}_GoodSchedule.csv -t {1} -c FC_{0}_{2}.txt'.format(date,tiles_file,telescope))

elif telescope=='Thacher':

	os.system('python CreateSchedule.py -d {0} -f {1} --obstele Thacher:Thacher --now {2} --start {3} --end {4}'.format(date,input_scheduler, args.now, args.start, args.end))

	print ("***************************")


	print ("Creating fields-center formatted file")


	os.system('python convert_scheduler_output_to_GWoutput.py -s Thacher_Thacher_{0}_GoodSchedule.csv -t {1} -c FC_{0}_{2}.txt'.format(date,tiles_file,telescope))

elif telescope=='Nickel':

	os.system('python CreateSchedule.py -d {0} -f {1} --obstele Lick:Nickel --now {2} --start {3} --end {4}'.format(date,input_scheduler, args.now, args.start, args.end))

	print ("***************************")


	print ("Creating fields-center formatted file")


	os.system('python convert_scheduler_output_to_GWoutput.py -s Lick_Nickel_{0}_GoodSchedule.csv -t {1} -c FC_{0}_{2}.txt'.format(date,tiles_file,telescope))

elif telescope=='Keck':

	os.system('python CreateSchedule.py -d {0} -f {1} --obstele Keck:Keck --now {2} --start {3} --end {4}'.format(date,input_scheduler, args.now, args.start, args.end))

	print ("***************************")


	print ("Creating fields-center formatted file")


	os.system('python convert_scheduler_output_to_GWoutput.py -s Keck_Keck_{0}_GoodSchedule.csv -t {1} -c FC_{0}_{2}.txt'.format(date,tiles_file,telescope))

os.system('python probs.py -f {0} --telescope {1} --date {2}'.format(tiles_file,telescope,date))"""
