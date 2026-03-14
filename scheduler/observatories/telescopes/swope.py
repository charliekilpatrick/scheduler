"""Swope telescope at Las Campanas Observatory."""

from scheduler.common import constants as C
from scheduler.common.target import TargetType
from scheduler.common.logging import get_logger
from scheduler.observatories.base import Telescope

import numpy as np

logger = get_logger(__name__)


class Swope(Telescope):

    def __init__(self):
        self.targets = None
        self.name = "Swope"
        self.fixed_exptime = None
        self.base_overhead = 60.
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
        if sn.obs_date is None or sn.ref_date is None:
            days_from_disc = 0.
        else:
            days_from_disc = (sn.obs_date - sn.ref_date).jd
        mag_reduction = days_from_disc*0.03
        adj_app_mag = sn.apparent_mag + mag_reduction
        logger.debug("days: %0.3f, mag red: %0.3f, adj mag: %0.3f", days_from_disc, mag_reduction, adj_app_mag)
        s_to_n = 30.
        g_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.g_band])
        r_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.r_band])
        i_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.i_band])
        V_exp = self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.V_band])
        mean_exp = np.mean([V_exp,g_exp,r_exp,i_exp])
        mean_exp = self.round_to_num(C.round_to, mean_exp)
        exposures.update({C.g_band: mean_exp})
        exposures.update({C.r_band: mean_exp})
        exposures.update({C.i_band: mean_exp})
        u_exp = self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.u_band]))
        B_exp = self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, adj_app_mag, self.filters[C.B_band]))
        if (mean_exp <= 540):
            logger.debug("Target Name: %s; u_exp: %s, mean_exp: %s", sn.name, u_exp, mean_exp)
            exposures.update({C.B_band: B_exp})
            exposures.update({C.V_band: mean_exp})
        if (mean_exp <= 300):
            exposures.update({C.u_band: u_exp})
        for key, value in exposures.items():
            if exposures[key] < 45:
                exposures[key] = 45
            elif exposures[key] > 600:
                exposures[key] = 600
        sn.exposures = exposures

    def compute_standard_exposure(self, std):
        exposures = {}
        s_to_n = 100
        exposures.update({C.u_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.u_band]))})
        exposures.update({C.B_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.B_band]))})
        exposures.update({C.V_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.V_band]))})
        exposures.update({C.g_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.g_band]))})
        exposures.update({C.r_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.r_band]))})
        exposures.update({C.i_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.i_band]))})
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
        s_to_n = 3.
        i_exp = self.time_to_S_N(s_to_n, gw.apparent_mag, self.filters[C.i_band])
        mean_exp = self.round_to_num(C.round_to, i_exp)
        if self.fixed_exptime:
            mean_exp = self.fixed_exptime
        exposures.update({C.i_band: mean_exp})
        gw.exposures = exposures

    def compute_exposures(self):
        fudge = self.base_overhead if len(self.targets) and getattr(self.targets[0], "exposures", None) else 60.0
        for tgt in self.targets:
            total_possible_time = np.sum(np.where(tgt.raw_airmass_array <= C.airmass_threshold)[0])
            if total_possible_time > 0:
                tgt.total_observable_min = total_possible_time * C.round_to / 60.0
                self.exp_funcs[tgt.type](tgt)
                fudge_factor = self.base_overhead if len(tgt.exposures) > 0 else 0
                total_minutes = (sum(tgt.exposures.values()) + fudge_factor) / 60.0
                tgt.total_minutes = round(total_minutes * 60.0 / C.round_to) * C.round_to / 60.0
                only_exposures = sum(tgt.exposures.values()) / 60.0
                tgt.total_minutes_only_exposures = round(only_exposures * 60.0 / C.round_to) * C.round_to / 60.0
            good_airmass = tgt.raw_airmass_array[np.asarray(np.where(tgt.raw_airmass_array <= C.airmass_threshold))]
            integrated_good_am = np.sum(good_airmass)
            if integrated_good_am > 0:
                tgt.total_good_air_mass = integrated_good_am
            if tgt.total_observable_min > 0:
                tgt.fraction_time_obs = float(tgt.total_minutes) / float(tgt.total_observable_min)

    def write_schedule(self, observatory_name, obs_date, targets, outdir=None,
        output_files=None, fieldcenters=None, pointing=None, one_off=False):
        if not output_files:
            output_files = "%s_%s_Schedule" % (self.name, obs_date.strftime('%Y%m%d'))
        output_rows = self.serialize_output_rows(targets, pointing=pointing)
        if fieldcenters:
            self.write_fieldcenters(targets, fieldcenters)
        self.write_csv_output(output_rows, output_files+'.csv')
        self.write_phot_file(targets, output_files+'.phot')
