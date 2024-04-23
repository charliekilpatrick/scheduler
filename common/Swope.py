import common.Constants as C
from common.Target import TargetType, Target
import common.Logs

from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np, csv, sys
from astropy.time import Time

import common.Telescope

# Used with Las Campanas Observatory
class Swope(common.Telescope.Telescope):

    def __init__(self):
        self.targets = None
        self.name = "Swope"

        self.fixed_exptime = None

        # Base overhead on observations (in seconds)
        self.base_overhead = 60.

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
        s_to_n = 3.
        i_exp = self.time_to_S_N(s_to_n, gw.apparent_mag, self.filters[C.i_band])
        mean_exp = self.round_to_num(C.round_to, i_exp)

        if self.fixed_exptime:
            mean_exp = self.fixed_exptime

        exposures.update({C.i_band: mean_exp})

        gw.exposures = exposures

    def compute_exposures(self):
        for tgt in self.targets:

            total_possible_time = np.sum(np.where(tgt.raw_airmass_array <= C.airmass_threshold)[0])

            if total_possible_time > 0:
                tgt.total_observable_min = total_possible_time * C.round_to/60.

                self.exp_funcs[tgt.type](tgt) # Sets exposures for each target by target type

                # per observatory - LCO Swope
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

    def swope_filter_row(self, exp_name, exp_time):
        filter_row = []
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(None)
        filter_row.append(exp_name)
        filter_row.append(exp_time)

        return filter_row

    def write_schedule(self, observatory_name, obs_date, targets, outdir=None,
        output_files=None, fieldcenters=None, pointing=None, **kwargs):

        if not output_files:
            output_files = "%s_%s_Schedule" % (self.name, obs_date.strftime('%Y%m%d'))

        output_rows = self.serialize_output_rows(targets, pointing=pointing)

        if fieldcenters:
            self.write_fieldcenters_file(targets, fieldcenters)

        self.write_csv_output(output_rows, output_files+'.csv')
        self.write_phot_file(targets, output_files+'.phot')
