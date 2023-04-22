import Constants as C
from Target import TargetType, Target
import Logs

from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np, csv, sys
from astropy.time import Time

from Telescope import Telescope

# Used with Lick Observatory
class Nickel(Telescope):

    def __init__(self):
        self.targets = None
        self.name = "Nickel"
        # Filter name: Zero-point
        self.filters = {
            C.B_band:14.40086907,
            C.V_band:14.94364672,
            C.r_prime:15.22269534,
            C.i_prime:15.01731372
        }
        self.exp_funcs = {
            TargetType.Supernova: self.compute_sn_exposure,
            TargetType.Template: self.compute_template_exposure,
            TargetType.Standard: self.compute_standard_exposure,
            TargetType.GW: self.compute_gw_exposure,
            TargetType.Force: self.compute_force_exposure
        }

    def set_targets(self, targets):
        self.targets = targets

    def get_targets(self):
        return self.targets

    def compute_force_exposure(self, sn):
        sn.exposures = sn.fixed_exp

    def compute_sn_exposure(self, sn):
        exposures = {}

        # Compute the current guess at apparent magnitude
        days_from_disc = (sn.obs_date - sn.disc_date).jd
        mag_reduction = days_from_disc*0.03
        adj_app_mag = sn.apparent_mag + mag_reduction

        # Change S/N depending on phase...
        s_to_n = 15. # base signal to noise

        r_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.r_prime])
        i_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.i_prime])
        V_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.V_band])
        B_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.B_band])

        if r_exp <= 900:
            exposures.update({C.r_prime: self.round_to_num(C.round_to, r_exp)})
        if i_exp <= 900:
            exposures.update({C.i_prime: self.round_to_num(C.round_to, i_exp)})
        if V_exp <= 900:
            exposures.update({C.V_band: self.round_to_num(C.round_to, V_exp)})
        if B_exp <= 900:
            exposures.update({C.B_band: self.round_to_num(C.round_to, B_exp)})

        # Finally, don't go less than 45s (~ readout time), don't go more than 600s on Swope
        for key, value in exposures.items():

            if exposures[key] < 180:
                exposures[key] = 180
            elif exposures[key] > 900:
                exposures[key] = 900

        sn.exposures = exposures

    def compute_standard_exposure(self, std):
        exposures = {}
        s_to_n = 100.

        # Don't know what the apparent mag should be?
        exposures.update({C.r_prime: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.ApparentMag, self.filters[C.r_prime]))})
        exposures.update({C.i_prime: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.ApparentMag, self.filters[C.i_prime]))})
        exposures.update({C.B_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.ApparentMag, self.filters[C.B_band]))})
        exposures.update({C.V_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.ApparentMag, self.filters[C.V_band]))})

        # Finally, don't go less than 10s for Nickel std, don't go more than 600s on Nickel
        for key, value in exposures.items():

            if exposures[key] < 10:
                exposures[key] = 10
            elif exposures[key] > 600:
                exposures[key] = 600

        std.exposures = exposures

    def compute_template_exposure(self, tmp):
        exposures = {}
        exposures.update({C.B_band: 1800})
        exposures.update({C.V_band: 1200})
        exposures.update({C.r_prime: 1200})
        exposures.update({C.i_prime: 1200})

        tmp.exposures = exposures


    def compute_gw_exposure(self, gw):
        exposures = {}
        exposures.update({C.r_prime: 300})

        gw.exposures = exposures

    def compute_exposures(self):
        for tgt in self.targets:

            if tgt.coord.dec.degree > 65:
                continue

            total_possible_time = np.sum(np.where(tgt.raw_airmass_array <= C.airmass_threshold)[0])

            if total_possible_time > 0:
                tgt.total_observable_min = total_possible_time * C.round_to/60.

                self.exp_funcs[tgt.type](tgt) # Sets exposures for each target by target type

                # per observatory - Nickel
                fudge_factor = 120. if len(tgt.exposures) > 0 else 0 # Build in a fudge factor based on # of exps GW

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

    def filter_row(self, exp_name, exp_time):
        filter_row = []
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(exp_name)
        filter_row.append(exp_time)

        return filter_row

    def write_schedule(self, observatory_name, obs_date, targets,
        output_files=None, fieldcenters=None, pointing=None):

        if output_files:
            file_to_write = output_files + '.csv'
            phot_file_to_write = output_files + '.phot'
            phot = open(phot_file_to_write, 'w')
        else:
            file_to_write = "%s_%s_%s_GoodSchedule.csv" % (observatory_name, self.name, obs_date.strftime('%Y%m%d'))
        if fieldcenters:
            fc = open(fieldcenters, 'w')
            fc.write('# field ampl ra dec epoch raD decD RAoffset DecOffset \n')

        # Initialize Google sheets methods for Nickel
        params = Logs.gsheets_params('nickel')
        sheet = Logs.initiate_gsheet(params['gsheets_token'])
        observers_sheet = params['OBSERVERS_SHEET']
        observer, night = Logs.check_if_tel_on_date(sheet, observers_sheet,
            obs_date)

        # Make an empty Nickel log
        success = Logs.copy_log(sheet, params['TEMPLATE'], 'Nickel Log',
            params['CURRENT_SHEET'], obs_date.strftime('%Y%m%d'))

        if not success:
            return(1)

        with open(file_to_write,'w') as csvoutput:
            writer = csv.writer(csvoutput, lineterminator="\n")

            output_rows = []

            # Append lines for pointing and focusing
            if pointing:
                for k,point in enumerate(pointing):
                    pointing_row = []
                    pointing_row.append('pointing'+str(k+1))
                    pointing_row.append('\''+point['ra'])
                    pointing_row.append('\''+point['dec'])
                    pointing_row.append('')
                    pointing_row.append('r\'')
                    pointing_row.append('1')
                    output_rows.append(pointing_row)

            focus_row = []
            focus_row.append('focus')
            focus_row.append('')
            focus_row.append('')
            focus_row.append('')
            focus_row.append('r\'')
            focus_row.append('')
            output_rows.append(focus_row)

            header_row = []
            header_row.append("Object Name")
            header_row.append("Right Ascension")
            header_row.append("Declination")
            header_row.append("Estimated Magnitude")
            header_row.append("Filter")
            header_row.append("Exposure Time")
            output_rows.append(header_row)

            for t in targets:

                # Ignore if we did not calculate any exposures for this target
                if len(t.exposures)==0:
                    continue

                hmsdms = t.coord.to_string(sep=':', style='hmsdms', precision=3)
                ra = hmsdms.split()[0]
                dec = hmsdms.split()[1]

                if t.type is not TargetType.GW:
                    tgt_row = []
                    tgt_row.append('\''+t.name.lower())
                    tgt_row.append('\''+ra)
                    tgt_row.append('\''+dec)
                    tgt_row.append('')
                    tgt_row.append('\''+C.r_prime)
                    tgt_row.append(10) # Acquisition in r'

                    if output_files:
                        for filt in t.exposures.keys():

                            phot_line = '{name} {name} {ra} {dec} {filt} '+\
                                '{exptime} {m3sigma} \n'

                            zeropoint = self.filters[filt]
                            exptime = t.exposures[filt]
                            m3sigma = self.limiting_magnitude(zeropoint, exptime, 3)

                            phot_line = phot_line.format(name=t.name, ra=ra,
                                dec=dec, filt=filt, exptime=exptime,
                                m3sigma=m3sigma)

                            phot.write(phot_line)


                    if fieldcenters:
                        fc_line = '{name:<40} {ampl} {ra_hms} {dec_dms} J2000  '+\
                            '{ra:>11}  {dec:>11}    0.0000000    0.0000000 \n'
                        fc.write(fc_line.format(name=t.name.lower(), ampl=1,
                            ra_hms=ra, dec_dms=dec, ra=t.coord.ra.degree,
                            dec=t.coord.dec.degree))

                    output_rows.append(tgt_row)

                    if C.B_band in t.exposures.keys():
                        output_rows.append(self.filter_row(C.B_band, t.exposures[C.B_band]))
                    if C.V_band in t.exposures.keys():
                        output_rows.append(self.filter_row(C.V_band, t.exposures[C.V_band]))
                    if C.r_prime in t.exposures.keys():
                        output_rows.append(self.filter_row(C.r_prime, t.exposures[C.r_prime]))
                    if C.i_prime in t.exposures.keys():
                        output_rows.append(self.filter_row(C.i_prime, t.exposures[C.i_prime]))

                else:
                    tgt_row = []
                    tgt_row.append('\''+t.name.lower())
                    tgt_row.append('\''+ra)
                    tgt_row.append('\''+dec)
                    tgt_row.append('')
                    tgt_row.append('\''+C.r_prime)
                    tgt_row.append(t.exposures[C.r_prime]) # Acquisition in r'

                    if output_files:
                        for filt in t.exposures.keys():

                            phot_line = '{name} {name} {ra} {dec} {filt} '+\
                                '{exptime} {m3sigma} \n'

                            zeropoint = self.filters[filt]
                            exptime = t.exposures[filt]
                            m3sigma = self.limiting_magnitude(zeropoint, exptime, 3)

                            phot_line = phot_line.format(name=t.name, ra=ra,
                                dec=dec, filt=filt, exptime=exptime,
                                m3sigma=m3sigma)

                            phot.write(phot_line)


                    if fieldcenters:
                        fc_line = '{name:<40} {ampl} {ra_hms} {dec_dms} J2000  '+\
                            '{ra:>11}  {dec:>11}    0.0000000    0.0000000 \n'
                        fc.write(fc_line.format(name=t.name.lower(), ampl=1,
                            ra_hms=ra, dec_dms=dec, ra=t.coord.ra.degree,
                            dec=t.coord.dec.degree))

                    output_rows.append(tgt_row)

            writer.writerows(output_rows)

            # Output log in Nickel log format to output dir
            target_list = 'nickel_targets_{0}'.format(obs_date.strftime('%Y%m%d'))
            target_list = params['target_dir'] + target_list
            with open(target_list, 'w') as f:
                for row in output_rows:
                    if row[0] and 'object' in row[0].lower():
                        continue
                    if row[0] and row[1] and row[2]:
                        target_list_name=str(row[0]).replace('\'','')
                        l='{0} {1} {2} {3} \n'.format(target_list_name,
                            row[1].replace('\'',''),
                            row[2].replace('\'',''),'2000')
                        f.write(l)



            # Need to append a dummy to the start of each row to account for PIC
            for i,row in enumerate(output_rows):
                output_rows[i] = [''] + row

            print('Adding schedule to Google sheet under:',obs_date.strftime('%Y%m%d'))
            print('Observer for tonight is:',observer)

            Logs.populate_nickel_log(sheet, params['CURRENT_SHEET'],
                obs_date.strftime('%Y%m%d'), output_rows,
                date=obs_date.strftime('%Y/%m/%d'), observer=observer,
                start_number=str(night * 1000 + 1))

        if fieldcenters:
            fc.close()

        if output_files:
            phot.close()


    def post_schedule(self, observatory_name, obs_date, targets):
        writer = csv.writer(csvoutput, lineterminator="\n")

        output_rows = []
        header_row = []
        header_row.append("Object Name")
        header_row.append("Right Ascension")
        header_row.append("Declination")
        header_row.append("Estimated Magnitude")
        header_row.append("Filter")
        header_row.append("Exposure Time")
        output_rows.append(header_row)

        for t in targets:

            # Ignore if we did not calculate any exposures for this target
            if len(t.exposures)==0:
                continue

            ra = t.coord.to_string(style='hmsdms', sep=':').split()[0]
            dec = t.coord.to_string(style='hmsdms', sep=':').split()[1]

            tgt_row = []
            tgt_row.append(t.name.lower())
            tgt_row.append(ra)
            tgt_row.append(dec)
            tgt_row.append(None)
            tgt_row.append(C.r_prime)
            tgt_row.append(10) # Acquisition in r'

            output_rows.append(tgt_row)

            if C.B_band in t.exposures.keys():
                output_rows.append(self.filter_row(C.B_band, t.exposures[C.B_band]))
            if C.V_band in t.exposures.keys():
                output_rows.append(self.filter_row(C.V_band, t.exposures[C.V_band]))
            if C.r_prime in t.exposures.keys():
                output_rows.append(self.filter_row(C.r_prime, t.exposures[C.r_prime]))
            if C.i_prime in t.exposures.keys():
                output_rows.append(self.filter_row(C.i_prime, t.exposures[C.i_prime]))

        writer.writerows(output_rows)
