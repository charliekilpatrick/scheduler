from enum import Enum

import numpy as np

class TargetType(Enum):
    Supernova = 1
    Template = 2
    Standard = 3
    GW = 4
    Force = 5
    NEWFIRM = 6

def handle_target_type(target):

    target_type = None
    disc_date = None
    if target['type'] == 'STD':
        target_type = TargetType.Standard
        disc_date = None
    elif target['type'] == 'NEWFIRM':
        target_type = TargetType.NEWFIRM
        disc_date = None
    elif target['type'] == 'TMP':
        target_type = TargetType.Template
    elif target['type'] == 'SN':
        target_type = TargetType.Supernova
    elif target['type'] == 'GW':
        target_type = TargetType.GW
    else:
        raise ValueError('Unrecognized target type!')

    return(target_type, disc_date)

class Target:
    def __init__(self, name, coord, priority, target_type, observatory_lat,
                 sidereal_radian_array, ref_date=None, apparent_mag=None,
                 obs_date=None,fixed_exp={}, orig_priority=None):
        # Provided by Constructor
        self.name = name
        self.coord = coord
        # If the target has a re-calculated priority, store the original value
        # as a separate quantity
        self.orig_priority = orig_priority
        self.priority = priority
        self.type = target_type
        # This is a reference date when the apparent magnitude applies
        self.ref_date = ref_date
        self.apparent_mag = apparent_mag
        self.obs_date = obs_date

        self.utc_start_time = None

        self.raw_airmass_array = None
        self.peak_airmass_idx = None
        self.peak_airmass_val = None

        self.initialize_airmass(observatory_lat, sidereal_radian_array)

        # Computed by Telescope
        self.net_priority = self.priority
        self.starting_index = 0 # Used to order net priority
        self.exposures = {} # Dictionary: filter:minutes
        self.fixed_exp = fixed_exp # Fixed exposures with filter:minutes
        self.total_observable_min = 0
        self.total_minutes = 0 # Total length of observation (inc. overhead)
        self.total_minutes_only_exposures = 0 # Total minutes minus overhead
        self.fraction_time_obs = 9999 # TotalMinutes / TotalObservableMin
        self.total_good_air_mass = 9999 # Proxy for elevation
        self.scheduled_time_array = None # Airmass plot abscissa
        self.scheduled_airmass_array = None # Airmass plot ordinate

    def initialize_airmass(self, latitude, sidereal_radian_array):

        # Computed by Constructor
        self.raw_airmass_array = self.compute_airmass(latitude,
            sidereal_radian_array)
        self.peak_airmass_idx = np.argmin(self.raw_airmass_array)
        self.peak_airmass_val = self.raw_airmass_array[self.peak_airmass_idx]

    # Compute airmass for an observatory at a specific latitude given an input
    # array of sidereal times (in radians).
    def compute_airmass(self, latitude, sidereal_radian_array, halimit=6.0):
        n = len(sidereal_radian_array)

        RA = np.empty(n)
        RA.fill(self.coord.ra.radian)

        DEC = np.empty(n)
        DEC.fill(self.coord.dec.radian)

        HA = sidereal_radian_array - RA
        HAhour = np.array(HA)  * 180.0/np.pi * 24.0/360.0
        HAhour[np.where(HAhour > 12.0)] = HAhour[np.where(HAhour > 12.0)]-24.0
        LAT = np.empty(n)
        LAT.fill(latitude)

        term1 = np.sin(DEC)*np.sin(LAT)
        term2 = np.cos(DEC)*np.cos(LAT)*np.cos(HA)
        am = 1.0/(np.sin(np.arcsin(term1+term2)))

        # Set airmass limit to some large value
        mask = (am > 10.0) | (am < 1.0)
        am[mask] = 9999

        # Hour angle limit
        am[np.abs(HAhour) > float(halimit)]=9999

        return am
