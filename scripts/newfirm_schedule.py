#!/usr/bin/env python
import os
import sys
from dateutil.parser import parse
import make_finders

def add_options(parser=None):
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('date', type=str,
        help='Date to make NEWFIRM schedule.')
    parser.add_argument('--first', default=False, action='store_true',
        help='Only make for first half.')
    parser.add_argument('--second', default=False, action='store_true',
        help='Only make for second half.')
    parser.add_argument('--compress', default=False, action='store_true',
        help='Compress output of NEWFIRM directory into a tarball.')
    parser.add_argument('--finders', default=False, action='store_true',
        help='Make finder charts from target file in NEWFIRM output dir.')
    parser.add_argument('--supernova-targets', default=None, type=str,
        help='Add file with supernova targets to list of input targets.')

    options = parser.parse_args()

    return(options)

def main():

    options = add_options()

    date = parse(options.date)
    datestr = date.strftime('%Y%m%d')

    ndir = os.environ['NEWFIRM_DIR']

    cmd = 'python CreateSchedule.py '
    cmd += f'-f {ndir}/targlists/newfirm_good_templates_{datestr}.ecsv '
    if options.supernova_targets:
        if os.path.exists(options.supernova_targets):
            print(f'Adding {options.supernova_targets}')
            cmd += f' {options.supernova_targets} '
    cmd += f'--date {options.date} --newfirm --minimize-slew '
    cmd += f'--outdir {ndir}/obsplans --obstele CTIO:Blanco '
    cmd += '--one-off '

    if options.first:
        cmd += '--first'
    elif options.second:
        cmd += '--second'

    print(cmd)
    os.system(cmd)

    if options.compress:
        outdir = f'newfirm_{datestr}'
        os.chdir(f'{ndir}/obsplans')
        cmd = f'tar --exclude .DS_Store -czvf {outdir}.tar.gz {outdir}'
        print(cmd)
        os.system(cmd)

    if options.finders:
        fulloutdir = os.path.join(ndir, 'obsplans', f'newfirm_{datestr}')
        fulloutfile = os.path.join(fulloutdir, f'newfirm_schedule_{datestr}.txt')

        outdir = os.path.join(ndir, 'finders', f'newfirm_{datestr}')

        print('Make finders')
        make_finders.make_finders(fulloutfile, directory=outdir,
            clobber=False)


if __name__=="__main__": main()
