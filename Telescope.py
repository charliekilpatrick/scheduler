import Constants as C
from Target import TargetType, Target
import Logs

from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np, csv, sys
from astropy.time import Time

class Telescope(metaclass=ABCMeta):

    @abstractmethod
    def set_targets(self, targets):
        pass

    @abstractmethod
    def get_targets(self):
        pass

    @abstractmethod
    def compute_exposures(self):
        pass

    @abstractmethod
    def write_schedule(self, observatory_name, obs_date, good_targets,
        output_files=None, fieldcenters=None, pointing=None):
        pass

    def round_to_num(self, round_to_num, input_to_round):
        return int(round_to_num*round(float(input_to_round)/round_to_num))

    def time_to_S_N(self, desired_s_n, apparent_mag, zeropoint, min_exp=45.):
        term1 = (desired_s_n / 3.)**2
        term2 = 0.4*(apparent_mag - zeropoint)
        exp_time = term1*10**term2

        if exp_time < min_exp:
            exp_time = min_exp

        return exp_time

    def limiting_magnitude(self, zeropoint, exptime, s_n):
        return(2.5 * np.log10(exptime / (s_n/3.)**2) + zeropoint)

    def compute_net_priorities(self):

        targets = self.get_targets()

        total_p = np.sum([t.priority for t in targets])
        print("Total Priority: %s" % total_p)

        total_good_time = np.sum([t.total_observable_min for t in targets])
        print("Total Good Time: %s" % total_good_time)

        total_exp_time = np.sum([t.total_minutes for t in targets])
        print("Total Exposure Time: %s" % total_exp_time)

        total_prob = 0
        if (total_p > 0 and total_good_time > 0 and total_exp_time > 0):
            for t in targets:
                frac_p = float(t.priority) / float(total_p)
                frac_time = float(t.total_observable_min)/float(total_good_time)
                frac_exp_time = (1.0-float(t.total_minutes)/float(total_exp_time))

                if (frac_exp_time == 0.0):
                    frac_exp_time = 1.0

                total_prob += frac_p*frac_time*frac_exp_time

            for t in targets:

                frac_p = float(t.priority) / float(total_p)
                frac_time = float(t.total_observable_min)/float(total_good_time)
                frac_exp_time = (1.0-float(t.total_minutes)/float(total_exp_time))

                t.net_priority = t.priority+((frac_p*frac_time*frac_exp_time)/total_prob)
                print("Nat: %s; Net: %0.5f" % (t.priority, t.net_priority))
        else:
            print("No valid targets...")

# Used with Las Campanas Observatory
class Swope(Telescope):

    def __init__(self):
        self.targets = None
        self.name = "Swope"

        #change as of Jun 4 after mirror cleaning
        self.filters = {
            C.u_band:15.05285437,
            C.B_band:16.24963648,
            C.V_band:16.2719428,
            C.g_band:16.4313935,
            C.r_band:16.43413935,
            C.i_band:15.86755865
        }

        self.exp_funcs = {
            TargetType.Supernova: self.compute_sn_exposure,
            TargetType.Template: self.compute_template_exposure,
            TargetType.Standard: self.compute_standard_exposure,
            TargetType.GW: self.compute_gw_exposure
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
        exposures.update({C.B_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.B_band]))})
        exposures.update({C.V_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.V_band]))})
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
        exposures.update({C.B_band: 1800})
        exposures.update({C.V_band: 1200})
        exposures.update({C.g_band: 1200})
        exposures.update({C.r_band: 1200})
        exposures.update({C.i_band: 1200})

        tmp.exposures = exposures

    def compute_gw_exposure(self, gw):
        exposures = {}

        # Compute gw exposure
        s_to_n = 5.
        r_exp = self.time_to_S_N(s_to_n, gw.apparent_mag, self.filters[C.r_band])
        mean_exp = self.round_to_num(C.round_to, r_exp)

        exposures.update({C.r_band: mean_exp})

        gw.exposures = exposures

    def compute_exposures(self):
        for tgt in self.targets:

            total_possible_time = np.sum(np.where(tgt.raw_airmass_array <= C.airmass_threshold)[0])

            if total_possible_time > 0:
                tgt.total_observable_min = total_possible_time * C.round_to/60.

                self.exp_funcs[tgt.type](tgt) # Sets exposures for each target by target type

                # per observatory - LCO Swope
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

    def swope_filter_row(self, exp_name, exp_time):
        filter_row = []
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(exp_name)
        filter_row.append(exp_time)

        return filter_row

    def write_schedule(self, observatory_name, obs_date, targets, output_files=None,
        fieldcenters=None, pointing=None):

        if output_files:
            file_to_write = output_files + '.csv'
            phot_file_to_write = output_files + '.phot'
            phot = open(phot_file_to_write, 'w')
        else:
            file_to_write = "%s_%s_%s_GoodSchedule.csv" % (observatory_name, self.name, obs_date.strftime('%Y%m%d'))
        if fieldcenters:
            fc = open(fieldcenters, 'w')
            fc.write('# field ampl ra dec epoch raD decD RAoffset DecOffset \n')
        with open(file_to_write,"w") as csvoutput:
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

            last_filter = C.r_band
            for t in targets:
                ra = t.coord.ra.hms
                dec = t.coord.dec.dms

                hmsdms = t.coord.to_string(sep=':', style='hmsdms', precision=3)
                ra_hms = hmsdms.split()[0]
                dec_dms = hmsdms.split()[1]

                if output_files:
                    for filt in t.exposures.keys():

                        phot_line = '{name} {name} {ra} {dec} {filt} '+\
                            '{exptime} {m3sigma} \n'

                        zeropoint = self.filters[filt]
                        exptime = t.exposures[filt]
                        m3sigma = self.limiting_magnitude(zeropoint, exptime, 3)

                        phot_line = phot_line.format(name=t.name, ra=ra_hms,
                            dec=dec_dms, filt=filt, exptime=exptime,
                            m3sigma=m3sigma)

                        phot.write(phot_line)

                if fieldcenters:
                    fc_line = '{name:<40} {ampl} {ra_hms} {dec_dms} J2000  '+\
                        '{ra:>11}  {dec:>11}    0.0000000    0.0000000 \n'
                    fc.write(fc_line.format(name=t.name.lower(), ampl='1',
                        ra_hms=ra_hms, dec_dms=dec_dms, ra=t.coord.ra.degree,
                        dec=t.coord.dec.degree))

                tgt_row = []
                tgt_row.append(t.name.lower())
                tgt_row.append("=\"%02d:%02d:%0.1f\"" % (ra[0],ra[1],ra[2]))

                dec_field = ("=\"%02d:%02d:%0.1f\"" % (dec[0],np.abs(dec[1]),np.abs(dec[2])))
                # Python has a -0.0 object. If the deg is this (because object lies < 60 min south), the string formatter will drop the negative sign
                if t.coord.dec < 0.0 and dec[0] == 0.0:
                    dec_field = ("=\"-0:%02d:%0.1f\"" % (np.abs(dec[1]),np.abs(dec[2])))

                tgt_row.append(dec_field)
                tgt_row.append(None)

                # Last criterion: if previous obj had full 6 filters, but this target only has 3
                if (last_filter == C.r_band) or \
                   (last_filter == C.i_band) or \
                    (last_filter == C.g_band) or \
                    (last_filter == C.B_band and len(t.exposures) < 6):

                    # tgt_row.append(C.r_band)   #comment out for GW
                    # tgt_row.append(10) # Acquisition in r #comment out for GW
                    output_rows.append(tgt_row)  #do not uncomment

                    # Start in riguVB order
                    output_rows.append(self.swope_filter_row(C.r_band, t.exposures[C.r_band]))  #comment out for GW
                    # output_rows.append(self.swope_filter_row(C.i_band, t.exposures[C.i_band]))
                    # output_rows.append(self.swope_filter_row(C.g_band, t.exposures[C.g_band]))  #comment out for GW
                    last_filter = C.g_band

                    if len(t.exposures) > 3:
                        if C.u_band in t.exposures:
                            output_rows.append(self.swope_filter_row(C.u_band, t.exposures[C.u_band]))

                        # output_rows.append(self.swope_filter_row(C.u_band, t.exposures[C.u_band]))
                        output_rows.append(self.swope_filter_row(C.V_band, t.exposures[C.V_band]))
                        output_rows.append(self.swope_filter_row(C.B_band, t.exposures[C.B_band]))
                        last_filter = C.B_band

                # Flip order: BVugir
                else:
                    tgt_row.append(C.B_band)
                    tgt_row.append(20) # Acquisition in B
                    output_rows.append(tgt_row)

                    output_rows.append(self.swope_filter_row(C.B_band, t.exposures[C.B_band]))
                    output_rows.append(self.swope_filter_row(C.V_band, t.exposures[C.V_band]))

                    if C.u_band in t.exposures:
                        output_rows.append(self.swope_filter_row(C.u_band, t.exposures[C.u_band]))

                    output_rows.append(self.swope_filter_row(C.g_band, t.exposures[C.g_band]))
                    output_rows.append(self.swope_filter_row(C.i_band, t.exposures[C.i_band]))
                    output_rows.append(self.swope_filter_row(C.r_band, t.exposures[C.r_band]))
                    last_filter = C.r_band

            writer.writerows(output_rows)

            if output_files:
                phot.close()

            if fieldcenters:
                fc.close()


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

class Thacher(Telescope):

    def __init__(self):
        self.targets = None
        self.name = "Thacher"
        # Filter name: Zero-point
        self.filters = {
            C.V_band:14.93966783,
            C.r_prime:15.1197765,
            C.i_prime:14.52342636,
            C.z_prime:13.68546838,
            C.g_prime:15.29722848
        }
        self.exp_funcs = {
            TargetType.Supernova: self.compute_sn_exposure,
            TargetType.Template: self.compute_template_exposure,
            TargetType.Standard: self.compute_standard_exposure,
            TargetType.GW: self.compute_gw_exposure
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

        # Change S/N depending on phase...
        s_to_n = 10 # base signal to noise
        if days_from_disc <= 10:
            s_to_n = 30
        elif days_from_disc > 10 and days_from_disc <= 60:
            s_to_n = 20

        exposures.update({C.r_prime: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.r_prime]))})
        exposures.update({C.i_prime: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.i_prime]))})

        # Only include these exposures if a relatively new SN
        if days_from_disc < 60:
            exposures.update({C.B_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.B_band]))})
            exposures.update({C.V_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.V_band]))})

        # Finally, don't go less than 45s (~ readout time), don't go more than 600s on Swope
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

        s_to_n = 5.
        r_exp = self.time_to_S_N(s_to_n, gw.apparent_mag, self.filters[C.r_prime])
        mean_exp = self.round_to_num(C.round_to, r_exp)

        exposures.update({C.r_prime: mean_exp})

        gw.exposures = exposures

    def compute_exposures(self):
        for tgt in self.targets:

            total_possible_time = np.sum(np.where(tgt.raw_airmass_array <= C.airmass_threshold)[0])

            if total_possible_time > 0:
                tgt.total_observable_min = total_possible_time * C.round_to/60.

                self.exp_funcs[tgt.type](tgt) # Sets exposures for each target by target type

                # per observatory
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

    def write_schedule(self, observatory_name, obs_date, targets, output_files=None,
        fieldcenters=None, pointing=None):

        if output_files:
            file_to_write = output_files + '.csv'
            phot_file_to_write = output_files + '.phot'
            phot = open(phot_file_to_write, 'w')
        else:
            file_to_write = "%s_%s_%s_GoodSchedule.csv" % (observatory_name, self.name, obs_date.strftime('%Y%m%d'))
        if fieldcenters:
            fc = open(fieldcenters, 'w')
            fc.write('# field ampl ra dec epoch raD decD RAoffset DecOffset \n')

        with open(file_to_write,"w") as csvoutput:
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

            last_filter = C.r_prime
            for t in targets:

                hmsdms = t.coord.to_string(sep=':', style='hmsdms', precision=3)
                ra = hmsdms.split()[0]
                dec = hmsdms.split()[1]

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

                tgt_row = []
                tgt_row.append(t.name.lower())
                tgt_row.append(ra)
                tgt_row.append(dec)
                tgt_row.append(None)

                # Last criterion: if previous obj had full 4 filters, but this target only has 2
                if (last_filter == C.r_prime) or \
                   (last_filter == C.i_prime) or \
                    (last_filter == C.B_band and len(t.exposures) < 4):

                    # tgt_row.append(C.r_prime)
                    # tgt_row.append(10) # Acquisition in r'
                    output_rows.append(tgt_row)

                    # Start in r'i'VB order
                    output_rows.append(self.filter_row(C.r_prime, t.exposures[C.r_prime]))
                    last_filter = C.i_prime

                    if len(t.exposures) > 2:
                        output_rows.append(self.filter_row(C.V_band, t.exposures[C.V_band]))
                        output_rows.append(self.filter_row(C.B_band, t.exposures[C.B_band]))
                        last_filter = C.B_band

                # Flip order: BVi'r'
                else:
                    tgt_row.append(C.B_band)
                    tgt_row.append(20) # Acquisition in B
                    output_rows.append(tgt_row)

                    output_rows.append(self.filter_row(C.B_band, t.exposures[C.B_band]))
                    output_rows.append(self.filter_row(C.V_band, t.exposures[C.V_band]))
                    output_rows.append(self.filter_row(C.i_prime, t.exposures[C.i_prime]))
                    output_rows.append(self.filter_row(C.r_prime, t.exposures[C.r_prime]))
                    last_filter = C.r_prime

            writer.writerows(output_rows)

        if output_files:
            phot.close()

        if fieldcenters:
            fc.close()
