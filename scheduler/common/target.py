"""Target types and Target class for scheduling (coordinates, airmass, exposures)."""

from enum import Enum
import numpy as np
from .logging import get_logger

logger = get_logger(__name__)


class TargetType(Enum):
    Supernova = 1
    Template = 2
    Standard = 3
    GW = 4
    Force = 5
    NEWFIRM = 6


def handle_target_type(target):
    """Map row type string (STD, NEWFIRM, TMP, SN, GW) to TargetType and disc_date (None)."""
    if target['type'] == 'STD':
        return TargetType.Standard, None
    if target['type'] == 'NEWFIRM':
        return TargetType.NEWFIRM, None
    if target['type'] == 'TMP':
        return TargetType.Template, None
    if target['type'] == 'SN':
        return TargetType.Supernova, None
    if target['type'] == 'GW':
        return TargetType.GW, None
    logger.error("Unrecognized target type: %s", target.get('type'))
    raise ValueError('Unrecognized target type!')


class Target:
    """Single observable target: coordinates, airmass array, exposures, priority."""

    def __init__(self, name, coord, priority, target_type, observatory_lat,
                 sidereal_radian_array, ref_date=None, apparent_mag=None,
                 obs_date=None, fixed_exp=None, orig_priority=None):
        if fixed_exp is None:
            fixed_exp = {}
        self.name = name
        self.coord = coord
        self.orig_priority = orig_priority
        self.priority = priority
        self.type = target_type
        self.ref_date = ref_date
        self.apparent_mag = apparent_mag
        self.obs_date = obs_date
        self.utc_start_time = None
        self.raw_airmass_array = None
        self.peak_airmass_idx = None
        self.peak_airmass_val = None
        self.initialize_airmass(observatory_lat, sidereal_radian_array)
        self.net_priority = self.priority
        self.starting_index = 0
        self.exposures = {}
        self.fixed_exp = fixed_exp
        self.total_observable_min = 0
        self.total_minutes = 0
        self.total_minutes_only_exposures = 0
        self.fraction_time_obs = 9999
        self.total_good_air_mass = 9999
        self.scheduled_time_array = None
        self.scheduled_airmass_array = None

    def initialize_airmass(self, latitude, sidereal_radian_array):
        """Compute airmass over night; set peak index/value."""
        self.raw_airmass_array = self.compute_airmass(latitude, sidereal_radian_array)
        self.peak_airmass_idx = np.argmin(self.raw_airmass_array)
        self.peak_airmass_val = self.raw_airmass_array[self.peak_airmass_idx]

    def compute_airmass(self, latitude, sidereal_radian_array, halimit=6.0):
        """Airmass at each time step; 9999 where below horizon or beyond halimit."""
        n = len(sidereal_radian_array)
        RA = np.empty(n)
        RA.fill(self.coord.ra.radian)
        DEC = np.empty(n)
        DEC.fill(self.coord.dec.radian)
        HA = sidereal_radian_array - RA
        HAhour = np.array(HA) * 180.0 / np.pi * 24.0 / 360.0
        HAhour[np.where(HAhour > 12.0)] = HAhour[np.where(HAhour > 12.0)] - 24.0
        LAT = np.empty(n)
        LAT.fill(latitude)
        term1 = np.sin(DEC) * np.sin(LAT)
        term2 = np.cos(DEC) * np.cos(LAT) * np.cos(HA)
        am = 1.0 / (np.sin(np.arcsin(term1 + term2)))
        mask = (am > 10.0) | (am < 1.0)
        am[mask] = 9999
        am[np.abs(HAhour) > float(halimit)] = 9999
        return am
