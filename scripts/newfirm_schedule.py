#!/usr/bin/env python
import os
import sys
from dateutil.parser import parse

def add_options(parser=None):
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('date', type=str,
        help='Date to make NEWFIRM schedule.')
    parser.add_argument('--first', default=False, action='store_true',
        help='Only make for first half.')
    parser.add_argument('--second', default=False, action='store_true',
        help='Only make for second half.')

    options = parser.parse_args()

    return(options)

def main():

    options = add_options()

    date = parse(options.date)
    datestr = date.strftime('%Y%m%d')

    ndir = os.environ['NEWFIRM_DIR']

    cmd = 'python CreateSchedule.py '
    cmd += f'-f {ndir}/targlists/newfirm_good_templates_{datestr}.ecsv '
    cmd += f'--date {options.date} --newfirm --minimize-slew '
    cmd += f'--outdir {ndir}/obsplans --obstele CTIO:Blanco '

    if options.first:
        cmd += '--first'
    elif options.second:
        cmd += '--second'

    print(cmd)
    os.system(cmd)

if __name__=="__main__": main()
