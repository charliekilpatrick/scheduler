"""
Observatory and Telescope base classes for the scheduler.

Observatory: site location, night boundaries, and scheduling.
Telescope: abstract base for exposure computation and schedule output.
"""

from scheduler.common import constants as C
from scheduler.common.utilities import UTC_Offset, download_ps1_catalog
from scheduler.common.logging import get_logger

import ephem
import requests
from datetime import timedelta, datetime, timezone
import numpy as np

from scheduler.common._datetime_parse import parse_datetime_ymd_hm
import operator
import copy
import os
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Column, Table
from astropy.time import Time, TimeDelta

logger = get_logger(__name__)


class Observatory:
    """
    An observatory site with telescopes and night/time arrays.

    Attributes
    ----------
    name : str
        Observatory key (e.g. LCO, Lick).
    lon, lat, elevation, horizon : str or float
        Site parameters for ephem.
    telescopes : dict
        Mapping telescope name -> Telescope instance.
    utc_offset, utc_offset_name : float, str
        Local time offset from UTC (hours) and name (e.g. PST).
    """

    def __init__(self, name, lon, lat, elevation, horizon, telescopes,
                 utc_offset, utc_offset_name, obs_date=None, start_time=None,
                 end_time=None, first=False, second=False):
        """
        Initialize an observatory and compute night boundaries for obs_date.

        Parameters
        ----------
        name : str
            Observatory key (e.g. 'LCO', 'Lick').
        lon, lat : str
            Longitude and latitude as strings for pyephem (e.g. '-70.69').
        elevation : float
            Elevation in meters.
        horizon : str
            Sun horizon in degrees for twilight (e.g. '-12').
        telescopes : dict
            Mapping telescope name -> Telescope instance.
        utc_offset : float
            Local time offset from UTC in hours (e.g. -8 for PST).
        utc_offset_name : str
            Timezone abbreviation (e.g. 'PST').
        obs_date : datetime or Time, optional
            Observation date; default is today UTC.
        start_time, end_time : str, optional
            Override night boundaries as 'HHMM'.
        first, second : bool, optional
            If True, use only first or second half of night.
        """
        self.name = name
        if obs_date is None:
            self.obs_date = datetime.now(timezone.utc)
        else:
            self.obs_date = Time(obs_date).datetime

        self.ephemeris = None
        self.telescopes = telescopes
        self.horizon = horizon
        self.lon = lon
        self.lat = lat
        self.elevation = elevation
        self.airmass_limit = C.airmass_threshold
        self.halimit = None
        self.utc_offset = utc_offset
        self.utc_offset_name = utc_offset_name
        self.pointing = []
        self.sidereal_string_array = []
        self.sidereal_radian_array = []
        self.utc_time_array = None
        self.local_time_array = None
        self.utc_begin_night = None
        self.utc_end_night = None
        self.utc_middle_night = None
        self.local_begin_night = None
        self.local_end_night = None
        self.local_middle_night = None
        self.length_of_night = 0.

        self.update_obs_date(self.obs_date, start_time=start_time,
            end_time=end_time, first=first, second=second)

    def update_obs_date(self, date, start_time=None, end_time=None,
                        first=False, second=False):
        """
        Recompute night boundaries and sidereal arrays for a given date.

        Parameters
        ----------
        date : datetime
            Observation date.
        start_time, end_time : str, optional
            Override night start/end as 'HHMM'.
        first, second : bool, optional
            If True, restrict to first or second half of night.
        """
        self.ephemeris = ephem.Observer()
        self.ephemeris.lon = self.lon
        self.ephemeris.lat = self.lat
        self.ephemeris.elevation = self.elevation
        self.ephemeris.horizon = self.horizon

        date_str = date.strftime("%Y-%m-%d")
        self.obs_date = parse_datetime_ymd_hm(f"{date_str} 12:00")

        self.ephemeris.date = (self.obs_date - timedelta(hours=self.utc_offset))

        self.utc_begin_night = self.ephemeris.next_setting(ephem.Sun(),
            use_center=True).datetime()
        self.utc_end_night = self.ephemeris.next_rising(ephem.Sun(),
            use_center=True).datetime()

        night_duration = self.utc_end_night - self.utc_begin_night
        self.utc_middle_night = self.utc_begin_night + 0.5*night_duration

        if first:
            self.utc_end_night = self.utc_middle_night
        elif second:
            self.utc_begin_night = self.utc_middle_night
        elif start_time is not None:
            hour = start_time[0:2]
            minute = start_time[2:4]
            if int(hour)<12:
                use_date = Time(date_str) + TimeDelta(86400*u.s)
                use_date_str = use_date.datetime.strftime('%Y-%m-%d')
            else:
                use_date_str = date_str
            self.utc_begin_night = parse_datetime_ymd_hm(
                f"{use_date_str} {hour}:{minute}"
            )
        elif end_time is not None:
            hour = end_time[0:2]
            minute = end_time[2:4]
            self.utc_end_night = parse_datetime_ymd_hm(f"{date_str} {hour}:{minute}")

        utc_tz = timezone.utc
        self.local_begin_night = self.utc_begin_night.replace(
            tzinfo=utc_tz
        ).astimezone(UTC_Offset(self.utc_offset, self.utc_offset_name))
        self.local_end_night = self.utc_end_night.replace(
            tzinfo=utc_tz
        ).astimezone(UTC_Offset(self.utc_offset, self.utc_offset_name))
        self.local_middle_night = self.utc_middle_night.replace(
            tzinfo=utc_tz
        ).astimezone(UTC_Offset(self.utc_offset, self.utc_offset_name))

        timeDiff = self.local_end_night - self.local_begin_night
        self.length_of_night = int(round(timeDiff.total_seconds()/C.round_to))

        self.utc_time_array = np.asarray([self.utc_begin_night +\
            timedelta(seconds=(i * C.round_to))
            for i in range(self.length_of_night)])

        self.local_time_array = np.asarray([self.local_begin_night +\
            timedelta(seconds=(i * C.round_to))
            for i in range(self.length_of_night)])

        self.ephemeris.date = self.utc_begin_night
        lst = self.ephemeris.sidereal_time()
        vunits = (u.hour, u.deg)
        c = SkyCoord(str(lst), self.ephemeris.lat, unit=vunits)

        try:
            table = Vizier.query_region(c, radius = 15 * u.deg, catalog='IV/22')[0]
            mask = table['Vmag']>4.0
            table = table[mask]
            separation = Column([c.separation(SkyCoord(r['RAJ2000'],
                r['DEJ2000'], unit=vunits)).degree for r in table],
                name='separation')
            table.add_column(separation)
            table.sort('separation')
            if len(table)>5:
                table = table[0:5]
                pointing = [SkyCoord(r['RAJ2000'], r['DEJ2000'], unit=vunits)
                    for r in table[0:5]]
            else:
                pointing = [SkyCoord(r['RAJ2000'], r['DEJ2000'], unit=vunits)
                    for r in table]
            for row in table:
                coord = SkyCoord(row['RAJ2000'],row['DEJ2000'], unit=vunits)
                self.pointing.append(
                        {'name': 'HD'+str(row['HD']),
                        'ra': coord.to_string(sep=':', style='hmsdms',
                            precision=3).split()[0],
                        'dec': coord.to_string(sep=':', style='hmsdms',
                            precision=3).split()[1],
                        'mag': row['Vmag']}
                        )
        except (requests.exceptions.ReadTimeout, requests.exceptions.ConnectTimeout,
                requests.exceptions.ProxyError, OSError, Exception) as e:
            logger.warning("Skipping pointing catalog (network or query failed): %s", e)

        self.sidereal_string_array = []
        self.sidereal_radian_array = []
        for utc_time in self.utc_time_array:
            self.ephemeris.date = utc_time
            st = self.ephemeris.sidereal_time()
            tokens = str(st).split(":")
            float_tokens = [float(t) for t in tokens]
            st_string = "%02d:%02d:%02d" % (float_tokens[0],float_tokens[1],
                float_tokens[2])
            self.sidereal_string_array.append(st_string)
            self.sidereal_radian_array.append(st)

    def set_airmass_threshold(self, threshold):
        """Set maximum allowed airmass for scheduling."""
        self.airmass_limit = threshold

    def set_halimit_threshold(self, halimit):
        """Set hour-angle limit for telescopes that use it."""
        self.halimit = halimit

    def swap(self, targets, idx):
        """
        Swap two adjacent targets in the list (for slew minimization).

        Parameters
        ----------
        targets : list of Target
            Current schedule.
        idx : int
            Index of first target to swap (swaps with idx+1).

        Returns
        -------
        list of Target
            New list with targets at idx and idx+1 swapped.
        """
        tgt1 = copy.deepcopy(targets[idx])
        tgt2 = copy.deepcopy(targets[idx + 1])
        tgt2.starting_index = tgt1.starting_index
        tgt1.starting_index = tgt1.starting_index + int(tgt2.total_minutes * 60.0 / C.round_to)
        newtargets = copy.deepcopy(targets)
        newtargets[idx] = tgt2
        newtargets[idx + 1] = tgt1
        return newtargets

    def compute_integrated_am(self, targets):
        """
        Sum integrated airmass over the scheduled window for all targets.

        Parameters
        ----------
        targets : list of Target
            Scheduled targets with starting_index and total_minutes set.

        Returns
        -------
        float
            Total integrated airmass.
        """
        integrated_am = 0.0
        for tgt in targets:
            start = int(tgt.starting_index)
            end = int(tgt.starting_index+tgt.total_minutes*60./C.round_to)
            integrated_am += np.sum(tgt.raw_airmass_array[start:end])
        return integrated_am

    def compute_night_fill_fraction(self, targets):
        """
        Log what fraction of the night is filled by the scheduled targets.

        Parameters
        ----------
        targets : list of Target
            Scheduled targets.

        Returns
        -------
        None
            Only logs; returns None if no targets.
        """
        if len(targets) == 0:
            logger.warning("No targets scheduled")
            return None
        num_elements = len(targets[0].raw_airmass_array)
        n = 0
        for tgt in targets:
            start = int(tgt.starting_index)
            end = int(tgt.starting_index+tgt.total_minutes*60./C.round_to)
            n += (end - start)
        pct = 100.0 * float(n) / num_elements
        logger.info("Filled %s/%s of scheduled (%.2f%%)", n, num_elements, pct)

    def is_allowed(self, targets):
        """
        Check whether all scheduled targets stay within airmass limit.

        Parameters
        ----------
        targets : list of Target
            Scheduled targets.

        Returns
        -------
        bool
            True if every time slot for every target is within airmass_limit.
        """
        for tgt in targets:
            goodtime = np.where(tgt.raw_airmass_array <= self.airmass_limit)[0]
            start = tgt.starting_index
            end = int(tgt.starting_index+tgt.total_minutes*60./C.round_to)
            for i in np.arange(start, end):
                if i not in goodtime:
                    return False
        return True

    def is_contiguous(self, int_array):
        """
        Return whether the integer array is contiguous (no gaps).

        Parameters
        ----------
        int_array : array-like
            Integer indices.

        Returns
        -------
        bool
            True if sorted array equals range(min, max+1).
        """
        i = iter(int_array)
        first = next(i)
        contiguous = all(a == b for a, b in enumerate(i, first + 1))
        return contiguous

    def get_start_of_observing_target(self, tgt):
        """
        Return the UTC time at which observing of this target starts.

        Parameters
        ----------
        tgt : Target
            Scheduled target with starting_index set.

        Returns
        -------
        astropy.time.Time
            Start time of the target's window.
        """
        date = Time(self.utc_begin_night)
        idx = tgt.starting_index
        newdate = date + 60.0 * u.s / C.round_to * idx
        return newdate

    def get_sun_twilight_local_times(self):
        """
        Return sunset, evening 12° twilight, morning 12° twilight, and sunrise in local time.
        Sunset is the evening before the schedule starts (before 12° twilight).

        Returns
        -------
        dict
            Keys: 'sunset', 'evening_twilight_12', 'morning_twilight_12', 'sunrise'.
            Values are timezone-aware datetime in local time.
        """
        utc_tz = timezone.utc
        evening_12_utc = self.utc_begin_night
        morning_12_utc = self.utc_end_night
        evening_12_local = evening_12_utc.replace(tzinfo=utc_tz).astimezone(
            UTC_Offset(self.utc_offset, self.utc_offset_name))
        morning_12_local = morning_12_utc.replace(tzinfo=utc_tz).astimezone(
            UTC_Offset(self.utc_offset, self.utc_offset_name))
        # Sunset and sunrise (horizon 0): use a time before evening twilight so next_setting/next_rising
        # give sunset/sunrise for the night we're plotting, not the next calendar day.
        saved_date = self.ephemeris.date
        self.ephemeris.date = (self.obs_date - timedelta(hours=self.utc_offset))
        self.ephemeris.horizon = '0'
        sunset_utc = self.ephemeris.next_setting(ephem.Sun(), use_center=True).datetime()
        sunrise_utc = self.ephemeris.next_rising(ephem.Sun(), use_center=True).datetime()
        self.ephemeris.horizon = self.horizon
        self.ephemeris.date = saved_date
        sunset_local = sunset_utc.replace(tzinfo=utc_tz).astimezone(
            UTC_Offset(self.utc_offset, self.utc_offset_name))
        sunrise_local = sunrise_utc.replace(tzinfo=utc_tz).astimezone(
            UTC_Offset(self.utc_offset, self.utc_offset_name))
        return {
            'sunset': sunset_local,
            'evening_twilight_12': evening_12_local,
            'morning_twilight_12': morning_12_local,
            'sunrise': sunrise_local,
        }

    def filter_targets(self, targets, cat_params=None, max_time=None):
        """
        Remove targets that cannot be observed or fail catalog criteria.

        Drops targets that exceed night length, have no exposure time, or
        (if cat_params given) have too few PS1 sources in the field.

        Parameters
        ----------
        targets : list of Target
            Candidate targets.
        cat_params : dict, optional
            If provided and non-empty, use catalog (cat_directory, cat_nsources, etc.)
            to filter by number of sources.
        max_time : float, optional
            Unused; reserved for future use.

        Returns
        -------
        list of Target
            Filtered list.
        """
        remove_targets = []
        if cat_params:
            for i,target in enumerate(targets):
                Mmax = cat_params['cat_maxmag']
                radius = cat_params['cat_radius']
                catdir = cat_params['cat_directory']
                filename = catdir + '/' + target.name + '.cat'
                if os.path.exists(filename):
                    cat = Table.read(filename, format='ascii')
                else:
                    logger.info(
                        "Downloading catalog for: %s (ra=%.4f, dec=%.4f)",
                        target.name, target.coord.ra.degree, target.coord.dec.degree,
                    )
                    cat = download_ps1_catalog(target, Mmax=22.0, radius=radius)
                    cat.write(filename, overwrite=True, format='ascii')
                logger.info("Got catalog: %s with %s sources", filename, len(cat))
                if len(cat) < cat_params['cat_nsources']:
                    logger.info("Removing %s (insufficient sources: %s)", target.name, len(cat))
                    remove_targets.append(target.name)
        for i, target in enumerate(targets):
            logger.debug(
                "%s: %s; %s min; Priority: %.2f",
                target.name, target.exposures, target.total_minutes, target.priority,
            )
            if ((self.length_of_night * C.round_to/60. < target.total_minutes) or
                (target.total_minutes_only_exposures <= 0)):
                remove_targets.append(target.name)
        targets = [t for t in targets if t.name not in remove_targets]
        return targets

    def get_slews(self, idxs, slews):
        """Sum total slew (degrees) for the given target index order and slew matrix."""
        total_slew = 0
        idx0 = idxs[0]
        for i,idx1 in enumerate(idxs):
            if i==0: continue
            total_slew += slews[idx0, idx1]
            idx0 = idx1
        return(total_slew)

    def swapPositions(self, l, pos1, pos2):
        """Swap elements at pos1 and pos2 in list l in place; return l."""
        l[pos1], l[pos2] = l[pos2], l[pos1]
        return l

    def create_slews(self, o):
        """Compute pairwise angular separation (degrees) matrix for targets o."""
        coords = SkyCoord(np.array([t.coord for t in o]))
        slews = []
        for c in coords:
            slews.append(list(c.separation(coords).to(u.deg).value))
        all_slews = np.array(slews)
        return(all_slews)

    def minimize_slews(self, o):
        """Greedily reorder targets to reduce total slew (in place). Returns reordered list."""
        all_slews = self.create_slews(o)
        curr_targets = copy.copy(o)
        idx = list(np.arange(len(o)))
        i=0
        curr_idx = list(copy.deepcopy(idx))
        iter_idx = list(copy.deepcopy(idx))
        curr_slews = self.get_slews(curr_idx, all_slews)
        for i in np.arange(len(idx)):
            if i==0: continue
            curr_i = i
            while curr_i > 0:
                iter_targets = copy.deepcopy(curr_targets)
                iter_idx = copy.deepcopy(curr_idx)
                iter_targets = self.swapPositions(iter_targets, curr_i-1, curr_i)
                iter_idx = self.swapPositions(iter_idx, curr_i-1, curr_i)
                iter_slews = self.get_slews(iter_idx, all_slews)
                if iter_slews < curr_slews:
                    curr_slews = iter_slews
                    curr_idx = copy.deepcopy(iter_idx)
                    curr_targets = copy.deepcopy(iter_targets)
                curr_i -= 1
        o = copy.copy(curr_targets)
        return o

    def _slew_delta_swap(self, slews, idxs, pos):
        """
        Return change in total slew if we swap positions pos and pos+1 in idxs.
        Positive means total slew would increase (worse).
        """
        n = len(idxs)
        if pos < 0 or pos + 1 >= n:
            return 0.0
        a, b = idxs[pos], idxs[pos + 1]
        if pos == 0:
            old_seg = slews[a, b] + (slews[b, idxs[2]] if n > 2 else 0.0)
            new_seg = slews[b, a] + (slews[a, idxs[2]] if n > 2 else 0.0)
        else:
            left = idxs[pos - 1]
            old_seg = slews[left, a] + slews[a, b] + (slews[b, idxs[pos + 2]] if pos + 2 < n else 0.0)
            new_seg = slews[left, b] + slews[b, a] + (slews[a, idxs[pos + 2]] if pos + 2 < n else 0.0)
        return new_seg - old_seg

    def minimize_slews(self, o):
        """
        Greedily reorder targets to reduce total slew (nearest-neighbor style).
        Works on indices only; reorders the list once at the end. O(N^2) with O(1) per swap.
        """
        n = len(o)
        if n <= 1:
            return o
        slews = self.create_slews(o)
        curr_idx = np.arange(n, dtype=int).tolist()
        for i in range(1, n):
            curr_i = i
            while curr_i > 0:
                delta = self._slew_delta_swap(slews, curr_idx, curr_i - 1)
                if delta < 0:
                    curr_idx[curr_i - 1], curr_idx[curr_i] = curr_idx[curr_i], curr_idx[curr_i - 1]
                curr_i -= 1
        return [o[j] for j in curr_idx]

    def run_slew_minimization(self, o):
        """Run slew minimization and log the improvement. Returns the reordered list."""
        logger.info("Minimizing slews...")
        if len(o) <= 1:
            return o
        slews = self.create_slews(o)
        idx = list(np.arange(len(o)))
        current_slew = self.get_slews(idx, slews)
        current_slew = float('%.4f' % current_slew)
        o = self.minimize_slews(o)
        idx = list(np.arange(len(o)))
        new_slew = self.get_slews(idx, slews)
        new_slew = float('%.4f' % new_slew)
        pct = float('%.2f' % ((1 - new_slew / current_slew) * 100))
        logger.info("Improved slews from %s deg -> %s deg (%.2f%%)", current_slew, new_slew, pct)
        return o

    def _assign_slots_in_order(self, o, length_of_night, C):
        """
        Re-assign time slots so that targets are observed in list order (slew order).
        Each target gets the earliest contiguous block that fits airmass, so schedule
        file order matches observation time order. Returns list with updated targets.
        """
        time_slots = np.zeros(length_of_night)
        out = []
        for tgt in o:
            gam1 = np.asarray(tgt.raw_airmass_array, dtype=float)
            gam2 = np.copy(gam1)
            gam2[np.where(time_slots == 1)] = 8888
            goodtime = np.where(gam2 <= self.airmass_limit)
            n = len(goodtime[0])
            need = int(tgt.total_minutes * 60.0 / C.round_to)
            first_indices = None
            current_start = -1
            for i in range(n):
                current_start += 1
                end = current_start + need
                if end > n:
                    break
                candidate = goodtime[0][current_start:end]
                if len(candidate) != need or not self.is_contiguous(candidate):
                    continue
                first_indices = candidate
                break
            if first_indices is None:
                logger.warning("Cannot re-fit %s in slew order; dropping.", tgt.name)
                continue
            time_slots[first_indices] = 1
            tgt.starting_index = int(first_indices[0])
            tgt.scheduled_airmass_array = np.asarray(tgt.raw_airmass_array)[first_indices]
            tgt.scheduled_time_array = np.asarray(self.local_time_array)[first_indices]
            out.append(tgt)
        return out

    def print_target_list(self, o):
        """Log each target's name, exposures, total minutes, priority, and coordinates."""
        for tgt in o:
            ra, dec = tgt.coord.to_string(style='hmsdms', sep=':', precision=3).split()
            logger.info(
                "%s: %s; %s min; Priority: %7.4f, RA: %s, Dec: %s",
                tgt.name, tgt.exposures, tgt.total_minutes, tgt.priority, ra, dec,
            )

    def schedule_targets(self, telescope_name, preview_plot=False,
                        output_files=None, fieldcenters=None, cat_params=None,
                        obs_date=None, start_time=None, end_time=None,
                        minimize_slew=False, outdir=None, first=False, second=False,
                        one_off=False):
        """
        Run the scheduling algorithm and write telescope output.

        Computes exposures and net priorities, assigns time slots by priority
        and airmass, optionally minimizes slew, then calls the telescope's
        write_schedule.

        Parameters
        ----------
        telescope_name : str
            Key in self.telescopes (e.g. 'Swope', 'Nickel').
        preview_plot : bool, optional
            Unused; reserved for plotting. Default False.
        output_files : str, optional
            Base name for schedule output files.
        fieldcenters : str, optional
            Path for field centers file.
        cat_params : dict, optional
            Catalog selection parameters (usecat, cat_directory, etc.).
        obs_date : datetime or Time, optional
            Override observation date.
        start_time, end_time : str, optional
            Override night boundaries as 'HHMM'.
        minimize_slew : bool, optional
            If True, run slew minimization on the schedule. Default False.
        outdir : str, optional
            Output directory for schedule files.
        first, second : bool, optional
            Restrict to first or second half of night.
        one_off : bool, optional
            Telescope-specific (e.g. NEWFIRM one-off mode). Default False.
        """
        if cat_params is None:
            cat_params = {}
        telescope = self.telescopes[telescope_name]
        telescope.compute_exposures()
        telescope.compute_net_priorities()
        targets = telescope.get_targets()

        if (obs_date is not None or start_time is not None or
            end_time is not None or first or second):

            self.__init__(self.name, self.lon, self.lat, self.elevation,
                self.horizon, self.telescopes, self.utc_offset,
                self.utc_offset_name, obs_date=obs_date,
                start_time=start_time, end_time=end_time,
                first=first, second=second)

            targets_copy = []
            for tgt in targets:
                tgt.initialize_airmass(self.ephemeris.lat,
                    self.sidereal_radian_array)
                targets_copy.append(copy.deepcopy(tgt))
            telescope.set_targets(copy.deepcopy(targets_copy))
            telescope.compute_exposures()
            telescope.compute_net_priorities()
            targets = telescope.get_targets()

        if obs_date is None:
            obs_date = self.obs_date

        logger.info("Start time %s", self.utc_begin_night)
        logger.info("End time %s", self.utc_end_night)

        def packable(goodtime, tgt):
            m_start = None
            goodtime_len = len(goodtime[0])*C.round_to/60.
            if goodtime_len < tgt.total_minutes:
                return None
            end=int(tgt.total_minutes*60./C.round_to)
            for m in range(len(goodtime[0][:-end])):
                start=int(m+tgt.total_minutes*60./C.round_to)
                if (goodtime[0][start] - goodtime[0][m])*C.round_to/60. < 45:
                    m_start = m
            return m_start

        def squeeze(tgt_i, targets, m_start):
            for tt in targets[:tgt_i]:
                affected_ind_up = goodtime[0][int(m_start + (tgt.total_minutes*60./C.round_to))]
                affected_ind_low = goodtime[0][m_start]
                if (tt.starting_index < affected_ind_up) and (tt.starting_index > affected_ind_low):
                    ind_increment = len(np.where(goodtime[0][m_start : m_start + int(tgt.total_minutes*60./C.round_to)] > tt.starting_index)[0])
                    tt.starting_index += ind_increment
            time_slots[affected_ind_low:affected_ind_up + 1] = 1
            return affected_ind_low

        targets.sort(key = operator.attrgetter('net_priority'))
        length_of_night = len(self.utc_time_array)
        targets = self.filter_targets(targets)
        logger.info("Scheduling targets...")
        time_slots = np.zeros(length_of_night)
        o = []
        bad_o = []

        if cat_params and 'usecat' in cat_params.keys() and cat_params['usecat']:
            targets = self.filter_targets(targets, cat_params=cat_params)

        for tgt_i, tgt in enumerate(targets):
            gam1 = copy.deepcopy(tgt.raw_airmass_array)
            gam2 = copy.deepcopy(tgt.raw_airmass_array)
            found = False
            while not found:
                gam2[np.where(time_slots == 1)] = 8888
                goodtime = np.where(gam2 <= self.airmass_limit)
                n = len(goodtime[0])
                current_start = -1
                best_indices = []
                largest_airmass = 1e+6
                for i in range(n):
                    current_start += 1
                    end = int(current_start + tgt.total_minutes*60./C.round_to)
                    candidate_indices = goodtime[0][current_start:end]
                    contiguous = self.is_contiguous(candidate_indices)
                    if len(candidate_indices) != int(tgt.total_minutes*60./C.round_to) or not contiguous:
                        continue
                    else:
                        integrated_am = np.sum(gam1[candidate_indices])
                        if integrated_am < largest_airmass and contiguous:
                            largest_airmass = integrated_am
                            best_indices = candidate_indices

                if largest_airmass < 1e+6:
                    found = True
                    time_slots[best_indices] = 1
                    tgt.scheduled_airmass_array = np.asarray(tgt.raw_airmass_array)[best_indices]
                    tgt.scheduled_time_array = np.asarray(self.local_time_array)[best_indices]
                    tgt.starting_index = best_indices[0]
                    o.append(tgt)
                    if not self.is_allowed(o):
                        o.remove(tgt)
                        found = False
                    logger.info("Successfully included field %i: %s %.3E", tgt_i, tgt.name, 1./tgt.priority/1000)
                else:
                    m_start = packable(goodtime, tgt)
                    logger.warning("Can't fit %s. Skipping!", tgt.name)
                    bad_o.append(tgt)
                    break

        logger.info("Current schedule:")
        o = sorted(o, key=lambda x: x.coord.ra.degree)
        self.print_target_list(o)
        if minimize_slew:
            o = self.run_slew_minimization(o)
            # Re-assign time slots in slew order so observation order = slew order
            o = self._assign_slots_in_order(o, length_of_night, C)
        logger.info("Final schedule:")
        self.print_target_list(o)
        if len(o) == 0:
            logger.warning("No targets scheduled; skipping write.")
            return []
        self.compute_night_fill_fraction(o)

        tgt0 = o[0]
        tgt1 = o[-1]
        date0 = self.get_start_of_observing_target(tgt0)
        date1 = self.get_start_of_observing_target(tgt1)
        d0 = date0.datetime.strftime('%Y-%m-%d')
        d1 = date1.datetime.strftime('%Y-%m-%d')
        logger.info("Writing out schedule: %s to %s", d0, d1)
        date = Time(obs_date)
        telescope.write_schedule(self.name, date, o, outdir=outdir,
            output_files=output_files, fieldcenters=fieldcenters,
            pointing=self.pointing, one_off=one_off)
        return o


# ---------------------------------------------------------------------------
# Telescope base class
# ---------------------------------------------------------------------------

from abc import ABCMeta, abstractmethod
import csv

from scheduler.common import constants as _C
from scheduler.common.target import TargetType, Target
from scheduler.common import logs  # noqa: F401


class Telescope(metaclass=ABCMeta):
    """
    Abstract base class for telescope-specific exposure and schedule output.

    Subclasses must implement set_targets, get_targets, compute_exposures,
    and write_schedule. They typically define filters (zeropoints), exp_funcs
    (per-target-type exposure logic), and override filter_row or write_schedule
    for telescope-specific file formats.
    """

    @abstractmethod
    def set_targets(self, targets):
        """Store the list of Target instances for this telescope."""
        pass

    @abstractmethod
    def get_targets(self):
        """Return the list of Target instances."""
        pass

    @abstractmethod
    def compute_exposures(self):
        """Compute exposure times for each target from airmass and type."""
        pass

    @abstractmethod
    def write_schedule(self, observatory_name, obs_date, good_targets,
                      output_files=None, fieldcenters=None, pointing=None):
        """Write schedule files (CSV, phot, etc.) for the good targets."""
        pass

    def set_base_overhead(self, overhead):
        """Set base overhead in seconds added per target."""
        self.base_overhead = overhead

    def set_fixed_exptime(self, exptime):
        """Override computed exposure time with a fixed value (seconds)."""
        self.fixed_exptime = exptime

    def round_to_num(self, round_to_num, input_to_round):
        """
        Round a value to the nearest multiple of round_to_num.

        Parameters
        ----------
        round_to_num : float
            Rounding step (e.g. 20 for 20 s).
        input_to_round : float
            Value to round.

        Returns
        -------
        int
            Rounded value.
        """
        return int(round_to_num * round(float(input_to_round) / round_to_num))

    def time_to_S_N(self, desired_s_n, apparent_mag, zeropoint, min_exp=45.0):
        """
        Compute exposure time (seconds) to reach desired S/N for a given magnitude.

        Parameters
        ----------
        desired_s_n : float
            Target signal-to-noise ratio.
        apparent_mag : float
            Apparent magnitude of the source.
        zeropoint : float
            Zero point for the filter.
        min_exp : float, optional
            Minimum exposure time in seconds. Default 45.

        Returns
        -------
        float
            Exposure time in seconds.
        """
        term1 = (desired_s_n / 3.0) ** 2
        term2 = 0.4 * (apparent_mag - zeropoint)
        exp_time = term1 * 10 ** term2
        if exp_time < min_exp:
            exp_time = min_exp
        return exp_time

    def filter_row(self, exp_name, exp_time):
        """Default filter row for CSV: [None, None, None, None, exp_name, exp_time]."""
        return [None, None, None, None, exp_name, exp_time]

    def _compute_exposures_loop(self, fudge_factor_sec=None):
        if fudge_factor_sec is None:
            fudge_factor_sec = getattr(self, "base_overhead", 60.0)
        for tgt in self.targets:
            if getattr(self, "_skip_target", lambda t: False)(tgt):
                continue
            total_possible_time = np.sum(
                np.where(tgt.raw_airmass_array <= _C.airmass_threshold)[0]
            )
            if total_possible_time > 0:
                tgt.total_observable_min = total_possible_time * _C.round_to / 60.0
                self.exp_funcs[tgt.type](tgt)
                total_minutes = (sum(tgt.exposures.values()) + fudge_factor_sec) / 60.0
                tgt.total_minutes = round(total_minutes * 60.0 / _C.round_to) * _C.round_to / 60.0
                only_exposures = sum(tgt.exposures.values()) / 60.0
                tgt.total_minutes_only_exposures = (
                    round(only_exposures * 60.0 / _C.round_to) * _C.round_to / 60.0
                )
            good_airmass = tgt.raw_airmass_array[
                np.asarray(np.where(tgt.raw_airmass_array <= _C.airmass_threshold))
            ]
            integrated_good_am = np.sum(good_airmass)
            if integrated_good_am > 0:
                tgt.total_good_air_mass = integrated_good_am
            if tgt.total_observable_min > 0:
                tgt.fraction_time_obs = float(tgt.total_minutes) / float(tgt.total_observable_min)

    def serialize_output_rows(self, targets, pointing=None, focus=False,
        do_acquisition=False):
        output_rows = []
        if pointing:
            for k,point in enumerate(pointing):
                pointing_row = []
                pointing_row.append('pointing'+str(k+1))
                pointing_row.append('\''+point['ra'])
                pointing_row.append('\''+point['dec'])
                pointing_row.append('\''+('%.1f'%point['mag']))
                pointing_row.append('r\'')
                pointing_row.append('1')
                output_rows.append(pointing_row)
        if focus:
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
            if len(t.exposures)==0:
                continue
            hmsdms = t.coord.to_string(sep=':', style='hmsdms', precision=3)
            ra = hmsdms.split()[0]
            dec = hmsdms.split()[1]
            mag=''
            if t.apparent_mag:
                mag = '%.1f'%t.apparent_mag
            if t.type is not TargetType.GW and do_acquisition:
                tgt_row = []
                tgt_row.append('\''+t.name.lower())
                tgt_row.append('\''+ra)
                tgt_row.append('\''+dec)
                tgt_row.append('\''+mag)
                tgt_row.append('\''+_C.r_prime)
                tgt_row.append(10)
                output_rows.append(tgt_row)
            else:
                tgt_row = []
                tgt_row.append('\''+t.name.lower())
                tgt_row.append('\''+ra)
                tgt_row.append('\''+dec)
                tgt_row.append('\''+mag)
                tgt_row.append('\'')
                tgt_row.append('\'')
                output_rows.append(tgt_row)
            for filt in sorted(t.exposures.keys()):
                output_rows.append(self.filter_row(filt, t.exposures[filt]))
        return output_rows

    def write_csv_output(self, output_rows, file_to_write):
        with open(file_to_write,"w") as csvoutput:
            writer = csv.writer(csvoutput, lineterminator="\n")
            writer.writerows(output_rows)

    def write_fieldcenters(self, targets, fieldcenters_file, ampl='1'):
        with open(fieldcenters_file, 'w') as fc:
            fc.write('# field ampl ra dec epoch raD decD RAoffset DecOffset \n')
            for t in targets:
                name = t.name
                ra = t.coord.ra.degree
                dec = t.coord.dec.degree
                hmsdms = t.coord.to_string(sep=':', style='hmsdms', precision=3)
                ra_hms = hmsdms.split()[0]
                dec_dms = hmsdms.split()[1]
                fc_line = f'{name:<40} {ampl} {ra_hms} {dec_dms} J2000  '+\
                          f'{ra:>11}  {dec:>11}    0.0000000    0.0000000 \n'
                fc.write(fc_line)

    def write_phot_file(self, targets, out_phot_file):
        with open(out_phot_file, 'w') as phot:
            phot.write('name ra dec filter exptime m3sigma priority \n')
            for t in targets:
                name = t.name
                ra = t.coord.ra.degree
                dec = t.coord.dec.degree
                for filt in t.exposures.keys():
                    zeropoint = self.filters[filt]
                    exptime = t.exposures[filt]
                    m3sigma = self.limiting_magnitude(zeropoint, exptime, 3)
                    priority = t.orig_priority
                    phot_line = f'{name} {ra} {dec} {filt} {exptime} {m3sigma} {priority} \n'
                    phot.write(phot_line)

    def limiting_magnitude(self, zeropoint, exptime, s_n):
        """
        Limiting magnitude for given zeropoint, exposure time, and S/N.

        Parameters
        ----------
        zeropoint : float
            Filter zero point.
        exptime : float
            Exposure time in seconds.
        s_n : float
            Signal-to-noise ratio.

        Returns
        -------
        float
            Limiting magnitude.
        """
        return 2.5 * np.log10(exptime / (s_n / 3.0) ** 2) + zeropoint

    def compute_net_priorities(self):
        """
        Compute net priority for each target from fraction of time and exposure.

        Rescales priority so that targets that use less of the night get a
        higher net priority. Logs total priority, good time, and exposure time.
        """
        targets = self.get_targets()
        total_p = np.sum([t.priority for t in targets])
        logger.info("Total Priority: %s", total_p)
        total_good_time = np.sum([t.total_observable_min for t in targets])
        logger.info("Total Good Time: %s", total_good_time)
        total_exp_time = np.sum([t.total_minutes for t in targets])
        logger.info("Total Exposure Time: %s", total_exp_time)
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
                logger.debug("Nat: %s; Net: %0.5f", t.priority, t.net_priority)
        else:
            logger.warning("No valid targets")
