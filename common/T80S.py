import common.Constants as C
from common.Target import TargetType, Target
import common.Logs

from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np, csv, sys
from astropy.time import Time

import common.Telescope

# Used with CTIO
class T80S(common.Telescope.Telescope):

    def __init__(self):
        self.targets = None
        self.name = "T80S"

        # Base overhead on observations (in seconds)
        self.base_overhead = 60.

        self.filters = {
            C.u_band:15.05285437,
            C.g_band:16.4313935,
            C.r_band:16.43413935,
            C.i_band:15.86755865,
            C.z_band:15.4,
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
        if sn.obs_date is None or sn.ref_date is None:
            days_from_disc = 0.
        else:
            days_from_disc = (sn.obs_date - sn.ref_date).jd
        mag_reduction = days_from_disc*0.03
        adj_app_mag = sn.apparent_mag + mag_reduction
        fmt = "days: %0.3f, mag red: %0.3f, adj mag: %0.3f"
        print(fmt % (days_from_disc, mag_reduction, adj_app_mag))

        # Change S/N depending on phase...
        s_to_n = 30. # base signal to noise

        g_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.g_band])
        r_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.r_band])
        i_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.i_band])
        z_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.z_band])

        exposures.update({C.g_band: g_exp})
        exposures.update({C.r_band: r_exp})
        exposures.update({C.i_band: i_exp})
        exposures.update({C.z_band: z_exp})

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
        exposures.update({C.g_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.g_band]))})
        exposures.update({C.r_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.r_band]))})
        exposures.update({C.i_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.i_band]))})
        exposures.update({C.z_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.z_band]))})


        # Round to nearest "exp_round_to" secs
        for key, value in exposures.items():
            if exposures[key] < 10:
                exposures[key] = 10
            elif exposures[key] > 600:
                exposures[key] = 600

        std.exposures = exposures

    def compute_template_exposure(self, tmp):
        exposures = {}
        exposures.update({C.g_band: 1200})
        exposures.update({C.r_band: 1200})
        exposures.update({C.i_band: 1200})
        exposures.update({C.z_band: 1500})

        tmp.exposures = exposures

    def compute_gw_exposure(self, gw, s_n=5, filts=[C.i_band]):
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

                # per observatory
                fudge_factor = self.base_overhead if len(tgt.exposures) > 0 else 0 # Build in a fudge factor based on # of exps GW

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

    def t80s_filter_row(self, exp_name, exp_time):
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
            file_to_write = "%s_%s_Schedule.csv" % (self.name, obs_date.strftime('%Y%m%d'))
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

                output_rows.append(tgt_row)  #do not uncomment

                # griz order
                if C.g_band in t.exposures.keys(): output_rows.append(self.t80s_filter_row(C.g_band, t.exposures[C.g_band]))  #comment out for GW
                if C.r_band in t.exposures.keys(): output_rows.append(self.t80s_filter_row(C.r_band, t.exposures[C.r_band]))  #comment out for GW
                if C.i_band in t.exposures.keys(): output_rows.append(self.t80s_filter_row(C.i_band, t.exposures[C.i_band]))
                if C.z_band in t.exposures.keys(): output_rows.append(self.t80s_filter_row(C.z_band, t.exposures[C.z_band]))


            writer.writerows(output_rows)

            if output_files:
                phot.close()

            if fieldcenters:
                fc.close()
