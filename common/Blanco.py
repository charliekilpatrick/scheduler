import common.Constants as C
from common.Target import TargetType, Target
import common.Logs

from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np, csv, sys
from astropy.time import Time

import common.Telescope
import os

# Used with Las Campanas Observatory
class Blanco(common.Telescope.Telescope):

    def __init__(self):
        self.targets = None
        self.name = "Blanco"

        #change as of Jun 4 after mirror cleaning
        self.filters = {
            C.J_band:15.0,
            C.H_band:15.0,
            C.K_band:15.0,
        }

        self.exp_funcs = {
            TargetType.Supernova: self.compute_sn_exposure,
            TargetType.Template: self.compute_template_exposure,
            TargetType.Standard: self.compute_standard_exposure,
            TargetType.GW: self.compute_gw_exposure,
            TargetType.NEWFIRM: self.compute_survey_exposure,
        }

    def set_targets(self, targets):
        self.targets = targets

    def get_targets(self):
        return self.targets

    def compute_sn_exposure(self, sn):
        exposures = {}
        exposures.update({C.J_band: 80})
        exposures.update({C.H_band: 80})
        exposures.update({C.K_band: 80})

        sn.exposures = exposures

    def compute_standard_exposure(self, std):
        exposures = {}

        std.exposures = exposures

    def compute_template_exposure(self, tmp):
        exposures = {}

        tmp.exposures = exposures

    def compute_survey_exposure(self, tmp):
        exposures = {}
        exposures.update({C.J_band: 80})

        tmp.exposures = exposures

    def compute_gw_exposure(self, gw, s_n=5, filts=[C.J_band]):
        exposures = {}

        # Compute gw exposure
        for filt in filts:
            exp = self.time_to_S_N(s_n, gw.apparent_mag, self.filters[filt])
            mean_exp = self.round_to_num(C.round_to, exp)

            exposures.update({filt: mean_exp})

        gw.exposures = exposures

    def compute_exposures(self):
        for tgt in self.targets:

            total_possible_time = np.sum(np.where(tgt.raw_airmass_array <= C.airmass_threshold)[0])

            if total_possible_time > 0:
                tgt.total_observable_min = total_possible_time * C.round_to/60.

                self.exp_funcs[tgt.type](tgt) # Sets exposures for each target by target type

                # per observatory - LCO Swope
                fudge_factor = 100. if len(tgt.exposures) > 0 else 0 # Build in a fudge factor based on # of exps GW

                total_minutes = (sum(tgt.exposures.values()) + fudge_factor)/60. # Sum total minutes
                tgt.total_minutes = round(total_minutes * 60./C.round_to) * C.round_to/60.
                only_exposures = sum(tgt.exposures.values())/60.
                tgt.total_minutes_only_exposures = round(only_exposures * 60./C.round_to) * C.round_to/60.

            good_airmass = tgt.raw_airmass_array[np.asarray(np.where(tgt.raw_airmass_array <= C.airmass_threshold))]
            integrated_good_am = np.sum(good_airmass)

            if integrated_good_am > 0:
                tgt.total_good_air_mass = integrated_good_am

            if tgt.total_observable_min > 0:
                tgt.fraction_time_obs = float(tgt.total_minutes)/float(tgt.total_observable_min)

    def swope_filter_row(self, exp_name, exp_time):
        filter_row = []
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(exp_name)
        filter_row.append(exp_time)

        return filter_row

    def initialize_obsfile(self, num, outdir='newfirm', subdir='scripts',
        observer='Charlie Kilpatrick', propid='2024A-937812',
        comment='Survey Field for NEWFIRM Infrared Survey for Transients'):

        num_str = str(num).zfill(3)
        filename = os.path.join(outdir, f'newfirm_sequence_{num_str}')
        file = open(filename,'w')

        file.write(f'PAN set obs.observer \"{observer}\" \n')
        file.write(f'PAN set obs.propid \"{propid}\" \n')
        file.write(f'PAN set image.comment \"{comment}\" \n')
        file.write(f'PAN set image.basename obj \n')

        return(file)

    def make_obs_sequence_line(self, repeats=1, fowler=1, coadd=4,
        exptime=5.0, filt='J', offsRA=0.0, offsDec=0.0, focus=None,
        objname='default'):

        if focus is None: focus = 'NONE'

        line  = f'{repeats}       {fowler}       '
        line += f'{exptime}       {filt}          '
        line += f'{offsRA}        {offsDec}       '
        line += f'{focus}      {objname} \n'

        return(name)

    def get_relative_focus(self, currfilt, to_filt):

        if currfilt=='J' and to_filt=='J':
            return('NONE')
        elif currfilt=='J' and to_filt=='H':
            return('100')
        elif currfilt=='J' and to_filt=='Ks':
            return('220')
        elif currfilt=='H' and to_filt=='J':
            return('-100')
        elif currfilt=='H' and to_filt=='H':
            return('NONE')
        elif currfilt=='H' and to_filt=='Ks':
            return('120')
        elif currfilt=='Ks' and to_filt=='J':
            return('-220')
        elif currfilt=='Ks' and to_filt=='H':
            return('-100')
        elif currfilt=='Ks' and to_filt=='Ks':
            return('NONE')

    def make_survey_sequence(self, file, name, currfilt='J'):

        focus = self.get_relative_focus(currfilt, 'J')

        file.write(f'repeats\tfowler\tcoadds\texp.time\tfiltername\toffsRA\toffsDEC\tfocus\tobject\n')
        file.write(f'1\t1\t4\t   5.000\tJX       \t-75.00\t0.00\t{focus}\t{name}\n')
        file.write(f'1\t1\t4\t   5.000\tJX       \t0.00\t75.00\tNONE\t{name}\n')
        file.write(f'1\t1\t4\t   5.000\tJX       \t75.00\t0.00\tNONE\t{name}\n')
        file.write(f'1\t1\t4\t   5.000\tJX       \t0.00\t-75.00\tNONE\t{name}\n')

        return(file, 'J')

    def make_supernova_sequence(self, file, name, exptime_per_filt=80.0,
        filts=['J','H','Ks'], currfilt='J'):

        exptime = exptime_per_filt/16.0
        first_filter = True

        file.write(f'repeats\tfowler\tcoadds\texp.time\tfiltername\toffsRA\toffsDEC\tfocus\tobject\n')

        if 'J' in filts:
            focus = self.get_relative_focus(currfilt, 'J')
            currfilt = 'J'

            if first_filter:
                offRA='450.00'
                offDec='450.00'
                first_filter = False
            else:
                offRA='-75.00'
                offDec='0.00'

            file.write(f'1\t1\t4\t   {exptime}\tJX       \t{offRA}\t{offDec}\t{focus}\t{name}\n')
            file.write(f'1\t1\t4\t   {exptime}\tJX       \t0.00\t75.00\tNONE\t{name}\n')
            file.write(f'1\t1\t4\t   {exptime}\tJX       \t75.00\t0.00\tNONE\t{name}\n')
            file.write(f'1\t1\t4\t   {exptime}\tJX       \t0.00\t-75.00\tNONE\t{name}\n')

        if 'H' in  filts:

            focus = self.get_relative_focus(currfilt, 'H')
            currfilt = 'H'

            if first_filter:
                offRA='450.00'
                offDec='450.00'
                first_filter = False
            else:
                offRA='-75.00'
                offDec='0.00'

            file.write(f'1\t1\t4\t   {exptime}\tHX       \t{offRA}\t{offDec}\t{focus}\t{name}\n')
            file.write(f'1\t1\t4\t   {exptime}\tHX       \t0.00\t75.00\tNONE\t{name}\n')
            file.write(f'1\t1\t4\t   {exptime}\tHX       \t75.00\t0.00\tNONE\t{name}\n')
            file.write(f'1\t1\t4\t   {exptime}\tHX       \t0.00\t-75.00\tNONE\t{name}\n')

        if 'Ks' in  filts:

            focus = self.get_relative_focus(currfilt, 'Ks')
            currfilt = 'Ks'

            if first_filter:
                offRA='450.00'
                offDec='450.00'
                first_filter = False
            else:
                offRA='-75.00'
                offDec='0.00'

            file.write(f'1\t1\t4\t   {exptime}\tKXs       \t{offRA}\t{offDec}\t{focus}\t{name}\n')
            file.write(f'1\t1\t4\t   {exptime}\tKXs       \t0.00\t75.00\tNONE\t{name}\n')
            file.write(f'1\t1\t4\t   {exptime}\tKXs       \t75.00\t0.00\tNONE\t{name}\n')
            file.write(f'1\t1\t4\t   {exptime}\tKXs       \t0.00\t-75.00\tNONE\t{name}\n')

        return(file, currfilt)

    def write_schedule(self, observatory_name, obs_date, targets, outdir=None,
        output_files=None, fieldcenters=None, pointing=None, 
        newfirm_dir='newfirm', obs_subdir='scripts', batch_size=10):

        if outdir is None:
            outdir = '.'

        # Update name of NEWFIRM subdir
        obs_date_str = obs_date.strftime('%Y%m%d')
        newfirm_dir = f'{newfirm_dir}_{obs_date_str}'

        if newfirm_dir and not os.path.exists(os.path.join(outdir, newfirm_dir)):
            os.makedirs(os.path.join(outdir, newfirm_dir))

        if not os.path.exists(os.path.join(outdir, newfirm_dir, obs_subdir)):
            os.makedirs(os.path.join(outdir, newfirm_dir, obs_subdir))

        if output_files:
            photfilename = os.path.join(outdir, newfirm_dir, f'{output_files}.txt')
            phot = open(photfilename, 'w')
        else:
            photfilename = f"newfirm_schedule_{obs_date_str}.txt"
            photfilename = os.path.join(outdir, newfirm_dir, photfilename)
            phot = open(photfilename, 'w')

        num_scripts = int(np.ceil(float(len(targets))/batch_size))
        curr_script_num = 1
        curr_script = self.initialize_obsfile(curr_script_num, 
            outdir=os.path.join(outdir, newfirm_dir), subdir=obs_subdir)

        # Assume that we start in J-band
        currfilt = 'J'

        for i,t in enumerate(targets):
            ra = t.coord.ra.hms
            dec = t.coord.dec.dms

            ra_deg = t.coord.ra.deg
            dec_deg = t.coord.dec.deg

            hmsdms = t.coord.to_string(sep=':', style='hmsdms', precision=3)
            ra_hms = hmsdms.split()[0]
            dec_dms = hmsdms.split()[1]

            curr_script.write(f'TCS slew {ra_hms} {dec_dms} \n')

            obs_sequence_file = f'{t.name}_sequence.obs'
            fullobsfile = os.path.join(outdir, newfirm_dir, obs_subdir, 
                obs_sequence_file)
            obs_sequence = open(fullobsfile, 'w')

            if t.type==TargetType.NEWFIRM:
                obs_sequence, currfilt = self.make_survey_sequence(
                    obs_sequence, t.name, currfilt=currfilt)
            elif t.type==TargetType.Supernova:
                obs_sequence, currfilt = self.make_supernova_sequence(
                    obs_sequence, t.name, currfilt=currfilt, 
                    exptime_per_filt=80.0, filts=['J','H','Ks'])
            elif t.type==TargetType.GW:
                obs_sequence, currfilt = self.make_supernova_sequence(
                    obs_sequence, t.name, currfilt=currfilt, 
                    exptime_per_filt=80.0, filts=['J'])

            obs_sequence.close()

            curr_script.write(f'seqrun.py -run {obs_sequence_file} \n')

            if (i+1)%batch_size==0 and i+1 < len(targets):
                curr_script.close()
                curr_script_num += 1 
                curr_script = self.initialize_obsfile(curr_script_num, 
                    outdir=os.path.join(outdir, newfirm_dir), subdir=obs_subdir)

            phot_line = '{idx}_{name} {ra} {dec} 2000.0 \n'
            phot_line = phot_line.format(idx=str(i+1).zfill(3), name=t.name, 
                ra=ra_hms,dec=dec_dms)
            phot.write(phot_line)

        curr_script.close()
        phot.close()
