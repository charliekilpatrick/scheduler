import common.Constants as C
import common.Telescope
from common.Utilities import UTC_Offset
from common.Utilities import download_ps1_catalog

import ephem
from dateutil.parser import parse
from datetime import tzinfo, timedelta, datetime
import pytz as pytz
import numpy as np
import operator
import copy
import sys
import os
import progressbar
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.dates as md
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Column, Table
from astropy.time import Time, TimeDelta

class Observatory():
    def __init__(self, name, lon, lat, elevation, horizon, telescopes,
        utc_offset, utc_offset_name, obs_date=None, start_time=None,
        end_time=None, first=False, second=False):

        self.name = name

        # Default is to assume we're schedule for now/tonight.  Otherwise, 
        # date can be passed when the schedule is created.
        if obs_date is None:
            self.obs_date = datetime.utcnow()
        else:
            self.obs_date = Time(obs_date).datetime

        self.ephemeris = None
        self.telescopes = telescopes
        self.horizon = horizon
        self.lon = lon
        self.lat = lat
        self.elevation = elevation

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

        self.ephemeris = ephem.Observer()
        self.ephemeris.lon = self.lon
        self.ephemeris.lat = self.lat
        self.ephemeris.elevation = self.elevation
        self.ephemeris.horizon = self.horizon

        date_str = date.strftime('%Y-%m-%d')

        obs_date = parse("%s 12:00" % date_str) # UTC Noon
        self.obs_date = obs_date

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
            self.utc_begin_night = parse("%s %s:%s" % (use_date_str, hour, minute))
        elif end_time is not None:
            hour = end_time[0:2]
            minute = end_time[2:4]
            self.utc_end_night = parse("%s %s:%s" % (date_str, hour, minute))

        self.local_begin_night = pytz.utc.localize(self.utc_begin_night) \
            .astimezone(UTC_Offset(self.utc_offset,self.utc_offset_name))
        self.local_end_night = pytz.utc.localize(self.utc_end_night) \
            .astimezone(UTC_Offset(self.utc_offset,self.utc_offset_name))
        self.local_middle_night = pytz.utc.localize(self.utc_middle_night) \
            .astimezone(UTC_Offset(self.utc_offset,self.utc_offset_name))

        # Total time in night
        timeDiff = self.local_end_night - self.local_begin_night

        # In units of C.round_to seconds, this is the length of the night
        self.length_of_night = int(round(timeDiff.total_seconds()/C.round_to))

        self.utc_time_array = np.asarray([self.utc_begin_night +\
            timedelta(seconds=(i * C.round_to))
            for i in range(self.length_of_night)])

        self.local_time_array = np.asarray([self.local_begin_night +\
            timedelta(seconds=(i * C.round_to))
            for i in range(self.length_of_night)])

        # Get a bright source for pointing
        self.ephemeris.date = self.utc_begin_night
        lst = self.ephemeris.sidereal_time()
        vunits = (u.hour, u.deg)
        c = SkyCoord(str(lst), self.ephemeris.lat, unit=vunits)

        # Get FK5 catalog and find closest object to zenith at start of night
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


        for utc_time in self.utc_time_array:
            self.ephemeris.date = utc_time
            st = self.ephemeris.sidereal_time()

            tokens = str(st).split(":")
            float_tokens = [float(t) for t in tokens]
            st_string = "%02d:%02d:%02d" % (float_tokens[0],float_tokens[1],
                float_tokens[2])

            self.sidereal_string_array.append(st_string)
            self.sidereal_radian_array.append(st)

    def swap(self, targets, idx):
        tgt1 = copy.deepcopy(targets[idx])
        tgt2 = copy.deepcopy(targets[idx+1])

        # Need to flip their starting indices around.  tgt2 is obvious since we
        # simply set the starting index to whatever tgt1 currently is
        tgt2.starting_index = tgt1.starting_index

        # For tgt1, need to account for when tgt2 will end
        tgt1.starting_index = tgt1.starting_index + int(tgt2.total_minutes * 60./C.round_to)

        newtargets = copy.deepcopy(targets)
        newtargets[idx] = tgt2
        newtargets[idx+1] = tgt1

        return(newtargets)

    def compute_integrated_am(self, targets):
        integrated_am = 0.

        for tgt in targets:
            start = int(tgt.starting_index)
            end = int(tgt.starting_index+tgt.total_minutes*60./C.round_to)
            integrated_am += np.sum(tgt.raw_airmass_array[start:end])

        return(integrated_am)

    def compute_night_fill_fraction(self, targets):
        integrated_am = 0.

        if len(targets)==0:
            print('\n\nNo targets scheduled')
            return(None)

        num_elements = len(targets[0].raw_airmass_array)

        n = 0
        for tgt in targets:
            start = int(tgt.starting_index)
            end = int(tgt.starting_index+tgt.total_minutes*60./C.round_to)
            n += (end - start)

        m = '\n\nFilled {n}/{m} of scheduled ({x}%)'.format(n=n, m=num_elements,
            x='%3.2f'%(100.0*float(n)/num_elements))
        print(m)

    # Check if the current schedule is allowed given constraints
    def is_allowed(self, targets):

        for tgt in targets:
            goodtime = np.where(tgt.raw_airmass_array <= C.airmass_threshold)[0]

            start = tgt.starting_index
            end = int(tgt.starting_index+tgt.total_minutes*60./C.round_to)
            for i in np.arange(start, end):
                if i not in goodtime:
                    return(False)

        return(True)

    def is_contiguous(self, int_array):
        i = iter(int_array)
        first = next(i)
        contiguous = all(a == b for a, b in enumerate(i, first + 1))
        return contiguous

    def get_start_of_observing_target(self, tgt):

        date = Time(self.utc_begin_night)
        # Schedule is segmented into seconds/C.round_to from start of the night
        idx = tgt.starting_index

        newdate = date + 60.0 * u.s  / C.round_to * idx

        return(newdate)

    # Filter targets based on certain criteria
    def filter_targets(self, targets, cat_params={}, max_time=None):

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
                    print('Downloading catalog for:',target.name,
                        target.coord.ra.degree, target.coord.dec.degree)
                    cat = download_ps1_catalog(target, Mmax=22.0, radius=radius)
                    cat.write(filename, overwrite=True, format='ascii')

                m='Got catalog: {0} with {1} sources'.format(filename, len(cat))
                print(m)

                if len(cat) < cat_params['cat_nsources']:
                    print('Removing:',target.name,len(cat))
                    remove_targets.append(target.name)

        for i,target in enumerate(targets):
            fmt = '{name}: {exposures}; {total} min; Priority: {priority}'
            fmt = fmt.format(name=target.name, exposures=target.exposures,
                total=target.total_minutes, priority='%.2f'%target.priority)

            print(fmt)

            # Get rid of targets whose observations can't fit in observable time
            # and with some exposure time
            if ((self.length_of_night * C.round_to/60. < target.total_minutes) or
                (target.total_minutes_only_exposures <= 0)):
                remove_targets.append(target.name)

        targets = [t for t in targets if t.name not in remove_targets]

        return(targets)

    def get_slews(self, idxs, slews):
        total_slew = 0
        idx0 = idxs[0]
        for i,idx1 in enumerate(idxs):
            if i==0: continue
            total_slew += slews[idx0, idx1]
            idx0 = idx1
        return(total_slew)

    def swapPositions(self, l, pos1, pos2):
        l[pos1], l[pos2] = l[pos2], l[pos1]
        return l

    def create_slews(self, o):

        coords = SkyCoord(np.array([t.coord for t in o]))
        slews = []
        for c in coords:
            slews.append(list(c.separation(coords).to(u.deg).value))
        all_slews = np.array(slews)

        return(all_slews)


    def minimize_slews(self, o):

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

        return(o)

    def run_slew_minimization(self, o):

        print('Minimizing slews...')
        all_slews = self.create_slews(o)
        idx = list(np.arange(len(o)))
        current_slew = self.get_slews(idx, all_slews)
        current_slew = float('%.4f'%current_slew)
        for i in np.arange(1):
            print(f'Iteration #{i+1}')
            o = self.minimize_slews(o)

        all_slews = self.create_slews(o)
        idx = list(np.arange(len(o)))
        new_slew = self.get_slews(list(np.arange(len(o))), all_slews)
        new_slew = float('%.4f'%new_slew)
        x = float('%.4f'%((1 - new_slew/current_slew)*100))

        print(f'Improved slews from {current_slew} deg->{new_slew} deg ({x}%)')

    def print_target_list(self, o):

        for tgt in o:
            fmt = '{name}: {exposures}; {total} min; Priority: {priority}, RA: {ra}, Dec: {dec}'
            ra, dec = tgt.coord.to_string(style='hmsdms', sep=':', precision=3).split()
            fmt = fmt.format(name=tgt.name, exposures=tgt.exposures,
                total=tgt.total_minutes, priority='%7.4f'%tgt.priority,
                ra=ra, dec=dec)
            print(fmt)

    def schedule_targets(self, telescope_name, preview_plot=False,
        output_files=None, fieldcenters=None, cat_params={}, obs_date=None,
        start_time=None, end_time=None, minimize_slew=False, outdir=None,
        first=False, second=False, one_off=False):

        # Update internal Target list with priorities and exposures
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

            # Recompute exposures and priorities for new telescope night
            targets_copy = []
            for tgt in targets:
                tgt.initialize_airmass(self.ephemeris.lat, 
                    self.sidereal_radian_array, halimit=tgt.halimit)
                targets_copy.append(copy.deepcopy(tgt))

            telescope.set_targets(copy.deepcopy(targets_copy))
            telescope.compute_exposures()
            telescope.compute_net_priorities()
            targets = telescope.get_targets()

        if obs_date is None:
            obs_date = self.obs_date

        print('Start time',self.utc_begin_night)
        print('End time',self.utc_end_night)

        def packable(goodtime, tgt):
            """Find nearby time intervals that amount to the observing time of a tile"""
            m_start = None
            goodtime_len = len(goodtime[0])*C.round_to/60.
            if goodtime_len < tgt.total_minutes:
                return(None)
            end=int(tgt.total_minutes*60./C.round_to)
            for m in range(len(goodtime[0][:-end])):
                # Exclude sparse intervals t_i - t_0 > 45 min
                start=int(m+tgt.total_minutes*60./C.round_to)
                if (goodtime[0][start] - goodtime[0][m])*C.round_to/60. < 45:
                    m_start = m
            return m_start

        def squeeze(tgt_i, targets, m_start):
            """ Move the scheduled time of the affected tiles to make room for a new tile """
            for tt in targets[:tgt_i]:
                affected_ind_up = goodtime[0][int(m_start + (tgt.total_minutes*60./C.round_to))]
                affected_ind_low = goodtime[0][m_start]
                if (tt.starting_index < affected_ind_up) and (tt.starting_index > affected_ind_low):
                    ind_increment = len(np.where(goodtime[0][m_start : m_start + int(tgt.total_minutes*60./C.round_to)] > tt.starting_index)[0])
                    tt.starting_index += ind_increment
                    # print(tt.starting_index)
            time_slots[affected_ind_low:affected_ind_up + 1] = 1
            return affected_ind_low

        # Sorted by priority and closeness to discovery
        targets.sort(key = operator.attrgetter('net_priority'))

        length_of_night = len(self.utc_time_array)
        targets = self.filter_targets(targets)

        print('\n\nScheduling targets:\n\n')
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

                gam2[np.where(time_slots == 1)] = 8888 # make different than airmass cutoff flag
                goodtime = np.where(gam2 <= C.airmass_threshold)
                n = len(goodtime[0])

                current_start = -1 # So that iterator below starts at 0, the first index
                best_indices = []
                largest_airmass = 1e+6

                # We are crawling forward along the array, grabbing segments of length "total_min",
                # and incrementing in starting index

                for i in range(n):
                    current_start += 1 # start index
                    # how many
                    end = int(current_start + tgt.total_minutes*60./C.round_to) 
                    # array of selected indices
                    candidate_indices = goodtime[0][current_start:end]
                    # Check if this associated integrated airmass corresponds 
                    # to a range of time that's contiguous
                    contiguous = self.is_contiguous(candidate_indices)
                    if len(candidate_indices) != int(tgt.total_minutes*60./C.round_to) or not contiguous: # If this is at the end of the array, it won't be the size we need
                        continue

                    else:
                        # Compute the integrated airmass. We're looking for the 
                        # smallest # => the best conditions
                        integrated_am = np.sum(gam1[candidate_indices])
                        if integrated_am < largest_airmass and contiguous:
                        # if this is the smallest, and is for a contiguous span 
                        # of time, it's the new one to beat
                            largest_airmass = integrated_am
                            best_indices = candidate_indices


                if largest_airmass < 1e+6:

                    found = True
                    time_slots[best_indices] = 1 # reserve these slots

                    # grab the corresponding
                    tgt.scheduled_airmass_array = np.asarray(tgt.raw_airmass_array)[best_indices]
                    tgt.scheduled_time_array = np.asarray(self.local_time_array)[best_indices]
                    tgt.starting_index = best_indices[0]

                    o.append(tgt)
                    if not self.is_allowed(o):
                        o.remove(tgt)
                        found = False
                    print('Successfully included field %i: %s %.3E' %(tgt_i, tgt.name, 1./tgt.priority/1000))
                else:
                    m_start = packable(goodtime, tgt)
                    print("Can't fit %s. Skipping!" % tgt.name)
                    bad_o.append(tgt)
                    break

        print('\n\nCurrent schedule:')
        o = sorted(o, key=lambda x: x.coord.ra.degree)
        self.print_target_list(o)

        # If we're minimizing slew, then run minimization algorithm
        if minimize_slew: self.run_slew_minimization(o)

        print('\n\nFinal schedule:')
        self.print_target_list(o)

        self.compute_night_fill_fraction(o)

        # Get date based on first and last target in schedule
        tgt0 = o[0]
        tgt1 = o[-1]

        date0 = self.get_start_of_observing_target(tgt0)
        date1 = self.get_start_of_observing_target(tgt1)

        d0 = date0.datetime.strftime('%Y-%m-%d')
        d1 = date1.datetime.strftime('%Y-%m-%d')

        print(f'\n\nWriting out schedule: {d0} to {d1}\n\n')

        date = Time(obs_date)

        telescope.write_schedule(self.name, date, o, outdir=outdir,
            output_files=output_files, fieldcenters=fieldcenters, 
            pointing=self.pointing, one_off=one_off)
