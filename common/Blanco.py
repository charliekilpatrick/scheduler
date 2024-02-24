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
            C.u_band:15.05285437,
            C.g_band:16.4313935,
            C.r_band:16.43413935,
            C.i_band:15.86755865,
            C.z_band:15.4,
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

        # Compute the current guess at apparent magnitude
        days_from_disc = (sn.obs_date - sn.disc_date).jd
        mag_reduction = days_from_disc*0.03
        adj_app_mag = sn.apparent_mag + mag_reduction
        fmt = "days: %0.3f, mag red: %0.3f, adj mag: %0.3f"
        print(fmt % (days_from_disc, mag_reduction, adj_app_mag))

        # Change S/N depending on phase...
        s_to_n = 30. # base signal to noise

        g_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.g_band])
        r_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.r_band])
        i_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.i_band])
        V_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.V_band])

        # Specific to Swope -- make Vgri the same length exposure...
        mean_exp = np.mean([V_exp,g_exp,r_exp,i_exp])
        mean_exp = self.round_to_num(C.round_to, mean_exp)

        exposures.update({C.g_band: mean_exp})
        exposures.update({C.r_band: mean_exp})
        exposures.update({C.i_band: mean_exp})

        u_exp = self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.u_band]))
        B_exp = self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.B_band]))

        # print (B_exp)

        if (mean_exp <= 540):
            print("Target Name: %s; u_exp: %s, mean_exp: %s" % (sn.name, u_exp, mean_exp))
            exposures.update({C.B_band: B_exp})
            exposures.update({C.V_band: mean_exp})

        if (mean_exp <= 300):
            exposures.update({C.u_band: u_exp})

        # Finally, adjust exposures
        for key, value in exposures.items():
            if exposures[key] < 45:
                exposures[key] = 45
            elif exposures[key] > 600:
                exposures[key] = 600

        sn.exposures = exposures

    def compute_standard_exposure(self, std):
        exposures = {}
        s_to_n = 100

        # Don't know what the apparent mag should be?
        exposures.update({C.u_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.u_band]))})
        exposures.update({C.g_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.g_band]))})
        exposures.update({C.r_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.r_band]))})
        exposures.update({C.i_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.i_band]))})

        # Finally, for standards round exps and don't go less than 10s, don't go more than 600s on Swope
        # Round to nearest "exp_round_to" secs
        for key, value in exposures.items():
            if exposures[key] < 10:
                exposures[key] = 10
            elif exposures[key] > 600:
                exposures[key] = 600

        std.exposures = exposures

    def compute_template_exposure(self, tmp):
        exposures = {}
        exposures.update({C.u_band: 1800})
        exposures.update({C.g_band: 1200})
        exposures.update({C.r_band: 1200})
        exposures.update({C.i_band: 1200})

        tmp.exposures = exposures

    def compute_survey_exposure(self, tmp):
        exposures = {}
        exposures.update({C.J_band: 80})

        tmp.exposures = exposures

    def compute_gw_exposure(self, gw, s_n=5, filts=[C.r_band]):
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
            photfilename = "%s_%s_%s_GoodSchedule.csv" % (observatory_name, 
                self.name, obs_date.strftime('%Y%m%d'))
            photfilename = os.path.join(outdir, newfirm_dir, photfilename)
            phot = open(photfilename, 'w')

        num_scripts = int(np.ceil(float(len(targets))/batch_size))
        curr_script_num = 1
        curr_script = self.initialize_obsfile(curr_script_num, 
            outdir=os.path.join(outdir, newfirm_dir), subdir=obs_subdir)

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

            if t.exposures[C.J_band]==80.0:
                obs_sequence.write(f'repeats\tfowler\tcoadds\texp.time\tfiltername\toffsRA\toffsDEC\tfocus\tobject\n')
                obs_sequence.write(f'1\t1\t4\t   5.000\tJX       \t-75.00\t0.00\tNONE\t{t.name}\n')
                obs_sequence.write(f'1\t1\t4\t   5.000\tJX       \t0.00\t75.00\tNONE\t{t.name}\n')
                obs_sequence.write(f'1\t1\t4\t   5.000\tJX       \t75.00\t0.00\tNONE\t{t.name}\n')
                obs_sequence.write(f'1\t1\t4\t   5.000\tJX       \t0.00\t-75.00\tNONE\t{t.name}\n')

            obs_sequence.close()

            curr_script.write(f'seqrun.py -run {obs_sequence_file} \n')

            if (i+1)%batch_size==0 and i+1 < len(targets):
                curr_script.close()
                curr_script_num += 1 
                curr_script = self.initialize_obsfile(curr_script_num, 
                    outdir=os.path.join(outdir, newfirm_dir), subdir=obs_subdir)

            if output_files:
                phot_line = '{idx}_{name} {ra} {dec} 2000.0 \n'
                phot_line = phot_line.format(idx=str(i+1).zfill(3), 
                    name=t.name, ra=ra_hms,dec=dec_dms)
                phot.write(phot_line)

        curr_script.close()

        if output_files:
            phot.close()
