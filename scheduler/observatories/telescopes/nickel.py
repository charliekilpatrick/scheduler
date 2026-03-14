"""Nickel telescope at Lick Observatory."""

from scheduler.common import constants as C
from scheduler.common.target import TargetType
from scheduler.common import logs
from scheduler.common import utilities
from scheduler.common.logging import get_logger
from scheduler.observatories.base import Telescope

import numpy as np
import os

from scheduler.common._datetime_parse import parse_date_flexible

logger = get_logger(__name__)


class Nickel(Telescope):

    def __init__(self):
        self.targets = None
        self.name = "Nickel"
        self.fixed_exptime = None
        self.base_overhead = 60.
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
        if sn.obs_date is None or sn.ref_date is None:
            days_from_disc = 0.
        else:
            days_from_disc = (sn.obs_date - sn.ref_date).jd
        mag_reduction = days_from_disc*0.03
        adj_app_mag = sn.apparent_mag + mag_reduction
        s_to_n = 15.
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
        for key, value in exposures.items():
            if exposures[key] < 90:
                exposures[key] = 90
        sn.exposures = exposures

    def compute_standard_exposure(self, std):
        exposures = {}
        s_to_n = 100.
        exposures.update({C.r_prime: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.r_prime]))})
        exposures.update({C.i_prime: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.i_prime]))})
        exposures.update({C.B_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.B_band]))})
        exposures.update({C.V_band: self.round_to_num(C.round_to, self.time_to_S_N(s_to_n, std.apparent_mag, self.filters[C.V_band]))})
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
                self.exp_funcs[tgt.type](tgt)
                fudge_factor = 120. if len(tgt.exposures) > 0 else 0
                total_minutes = (sum(tgt.exposures.values()) + fudge_factor)/60.
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

    def initiate_nickel_sheet(self, obs_date, is_gw=False):
        if ('NICKEL_SCHEDULE_SHEET' in os.environ.keys() and 
            'GSHEETS_TOKEN' in os.environ.keys() and
            'NICKEL_TEMPLATE_SHEET' in os.environ.keys() and
            'NICKEL_OBSERVERS_SHEET' in os.environ.keys()):
            params = Logs.gsheets_params(os.environ['NICKEL_SCHEDULE_SHEET'], 
                os.environ['NICKEL_TEMPLATE_SHEET'], 
                os.environ['GSHEETS_TOKEN'],
                observers_sheet=os.environ['NICKEL_OBSERVERS_SHEET'])
            sheet = logs.initiate_gsheet(params['gsheets_token'])
            if not is_gw:
                (observer, night) = logs.check_if_tel_on_date(
                    sheet,
                    params["OBSERVERS_SHEET"],
                    parse_date_flexible(obs_date.strftime("%Y%m%d")),
                )
                tab_name = obs_date.strftime('%Y%m%d')
                if night==0: night = 1
            else:
                observer = ''
                night = 1
                tab_name = obs_date.strftime('%Y%m%d')+'gw'
            success = logs.copy_log(sheet, params['TEMPLATE'], 'Nickel Log',
                params['CURRENT_SHEET'], tab_name)
        else:
            success = False
            observer = None
            params = None
            tab_name = None
            night = None
            sheet = None
        return(success, observer, params, tab_name, night, sheet)

    def write_schedule(self, observatory_name, obs_date, targets, outdir=None,
        output_files=None, fieldcenters=None, pointing=None, **kwargs):
        is_gw = all([t.type is TargetType.GW for t in targets])
        if not output_files:
            output_files = "%s_%s_Schedule" % (self.name, obs_date.strftime('%Y%m%d'))
        log_success, observer, params, tab_name, night, sheet = self.initiate_nickel_sheet(obs_date, is_gw=is_gw)
        output_rows = self.serialize_output_rows(targets, pointing=pointing,
            focus=True, do_acquisition=True)
        for i,row in enumerate(output_rows):
            output_rows[i] = [''] + row
        if params is not None:
            logger.info("Adding schedule to Google sheet under: %s", obs_date.strftime('%Y%m%d'))
            logger.info("Observer for tonight is: %s", observer)
            if log_success:
                logs.populate_nickel_log(sheet, params['CURRENT_SHEET'], tab_name, 
                    output_rows, date=obs_date.strftime('%Y/%m/%d'), 
                    observer=observer, start_number=str(night * 1000 + 1))
        if fieldcenters:
            self.write_fieldcenters(targets, fieldcenters)
        self.write_csv_output(output_rows, output_files+'.csv')
        self.write_phot_file(targets, output_files+'.phot')
