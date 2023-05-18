import Constants as C
import Telescope
from Utilities import UTC_Offset

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
from astroquery.mast import Catalogs
from astropy.time import Time

def download_ps1catalog(target, Mmax=18.0, radius=0.05):
        coord = target.coord
        cat = Catalogs.query_region(str(coord.ra.degree)+' '+\
            str(coord.dec.degree), radius=2*radius,
            catalog="Panstarrs", data_release="dr2", table="mean")

        cat = cat.filled()
        mask = (cat['raStack'] < 360.0) & (cat['decStack'] < 90.0)
        cat = cat[mask]

        coords = SkyCoord(cat['raStack'], cat['decStack'], unit='deg')
        mask = coords.separation(coord) < radius * u.degree
        cat = cat[mask]

        mask = (cat['gMeanPSFMag'] < Mmax) & (cat['gMeanPSFMag'] > 0) &\
               (cat['rMeanPSFMag'] < Mmax) & (cat['rMeanPSFMag'] > 0) &\
               (cat['iMeanPSFMag'] < Mmax) & (cat['iMeanPSFMag'] > 0) &\
               (cat['zMeanPSFMag'] < Mmax) & (cat['zMeanPSFMag'] > 0) &\
               (cat['yMeanPSFMag'] < Mmax) & (cat['yMeanPSFMag'] > 0)
        cat = cat[mask]

        return(cat)

class Observatory():
    def __init__(self, name, lon, lat, elevation, horizon, telescopes,
        utc_offset, utc_offset_name, obs_date=None):

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
        self.local_begin_night = None
        self.local_end_night = None

        self.length_of_night = 0.

        self.update_obs_date(self.obs_date)

    def update_obs_date(self, date):

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

        self.local_begin_night = pytz.utc.localize(self.utc_begin_night) \
            .astimezone(UTC_Offset(self.utc_offset,self.utc_offset_name))
        self.local_end_night = pytz.utc.localize(self.utc_end_night) \
            .astimezone(UTC_Offset(self.utc_offset,self.utc_offset_name))

        # Total time in night
        timeDiff = self.local_end_night - self.local_begin_night
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
        table = Vizier.query_region(c, radius = 10 * u.deg, catalog='IV/22')[0]
        separation = Column([c.separation(SkyCoord(r['RAJ2000'],
            r['DEJ2000'], unit=vunits)).degree for r in table],
            name='separation')

        table.add_column(separation)
        table.sort('separation')

        if len(table)>5:
            pointing = [SkyCoord(r['RAJ2000'], r['DEJ2000'], unit=vunits)
                for r in table[0:5]]
        else:
            pointing = [SkyCoord(r['RAJ2000'], r['DEJ2000'], unit=vunits)
                for r in table]

        self.pointing = [
                {'ra': p.to_string(sep=':', style='hmsdms',
                    precision=3).split()[0],
                'dec': p.to_string(sep=':', style='hmsdms',
                    precision=3).split()[1]} for p in pointing]


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
        tgt2.starting_index = targets[idx].starting_index

        # For tgt1, need to account for when tgt2 will end
        tgt1.starting_index = targets[idx].starting_index + int(tgt2.total_minutes * 60./C.round_to)

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

    def schedule_targets(self, telescope_name, preview_plot=False, asap=False,
        output_files=None, fieldcenters=None, cat_params={}, obs_date=None):

        # Update internal Target list with priorities and exposures
        telescope = self.telescopes[telescope_name]
        telescope.compute_exposures()
        telescope.compute_net_priorities()
        targets = telescope.get_targets()

        if obs_date is not None:

            self.__init__(self.name, self.lon, self.lat, self.elevation, 
                self.horizon, self.telescopes, self.utc_offset, 
                self.utc_offset_name, obs_date=obs_date)

            # Recompute exposures and priorities for new telescope night
            targets_copy = []
            for tgt in targets:
                tgt.initialize_airmass(self.lat, self.sidereal_radian_array,
                    halimit=tgt.halimit)
                targets_copy.append(tgt)

            telescope.targets = targets_copy
            telescope.compute_exposures()
            telescope.compute_net_priorities()
            targets = telescope.get_targets()

        def packable(goodtime, tgt):
            """Find nearby time intervals that amount to the observing time of a tile"""
            m_start = None
            goodtime_len = len(goodtime[0])*C.round_to/60.
            if goodtime_len < tgt.total_minutes:
                return(None)
            end=int(tgt.total_minutes*60./C.round_to-1.)
            for m in range(len(goodtime[0][:-end])):
                # Exclude sparse intervals t_i - t_0 > 45 min
                start=int(m+tgt.total_minutes*60./C.round_to-1)
                if (goodtime[0][start] - goodtime[0][m])*C.round_to/60. < 45:
                    m_start = m
            return m_start

        def squeeze(tgt_i, targets, m_start):
            """ Move the scheduled time of the affected tiles to make room for a new tile """
            for tt in targets[:tgt_i]:
                affected_ind_up = goodtime[0][int(m_start + (tgt.total_minutes*60./C.round_to - 1))]
                affected_ind_low = goodtime[0][m_start]
                if (tt.starting_index < affected_ind_up) and (tt.starting_index > affected_ind_low):
                    ind_increment = len(np.where(goodtime[0][m_start : m_start + int(tgt.total_minutes*60./C.round_to)] > tt.starting_index)[0])
                    tt.starting_index += ind_increment
                    # print(tt.starting_index)
            time_slots[affected_ind_low:affected_ind_up + 1] = 1
            return affected_ind_low

        # Sorted by priority and closeness to discovery
        if asap:
            targets.sort(key = operator.attrgetter('priority'))
        else:
            targets.sort(key = operator.attrgetter('net_priority'))

        length_of_night = len(self.utc_time_array) # In minutes

        targets_copy = []
        for tgt in targets:
            fmt = '{name}: {exposures}; {total} min; Priority: {priority}'
            fmt = fmt.format(name=tgt.name, exposures=tgt.exposures,
                total=tgt.total_minutes, priority='%7.4f'%tgt.priority)
            print(fmt)

            # Get rid of targets whose observations can't fit in observable time
            # and with some exposure time
            if (tgt.total_observable_min >= tgt.total_minutes and
                tgt.total_minutes_only_exposures>0):
                targets_copy.append(tgt)

        print('\n\nScheduling targets:\n\n')

        targets = targets_copy
        targets.sort(key = operator.attrgetter('priority'))

        time_slots = np.zeros(length_of_night)
        o = []
        bad_o = []

        if cat_params:
            if cat_params['usecat']:
                print('\n\nUsing catalogs...')
                remove_targets = []
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
                        cat = download_ps1catalog(target, Mmax=22.0,
                            radius=radius)
                        cat.write(filename, overwrite=True, format='ascii')

                    mask = cat['gMeanPSFMag'] < Mmax
                    cat = cat[mask]

                    m='Got catalog: {0} with {1} sources'.format(filename,
                        len(cat))
                    print(m)

                    if len(cat) < cat_params['cat_nsources']:
                        print('Removing:',target.name,len(cat))
                        remove_targets.append(target)

                targnames = [t.name for t in remove_targets]
                targets = [t for t in targets if t not in remove_targets]

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
                    end = int(current_start + tgt.total_minutes*60./C.round_to) # how many
                    candidate_indices = goodtime[0][current_start:end] # array of selected indices
                    # Check if this associated integrated airmass corresponds to a range of time that's contiguous
                    contiguous = self.is_contiguous(candidate_indices)
                    if len(candidate_indices) != int(tgt.total_minutes*60./C.round_to) or not contiguous: # If this is at the end of the array, it won't be the size we need
                        continue

                    else:
                        # Compute the integrated airmass. We're looking for the smallest # => the best conditions
                        integrated_am = np.sum(gam1[candidate_indices])
                        if asap:
                            if (integrated_am <= (C.airmass_threshold * tgt.total_minutes *60./C.round_to)):
                                largest_airmass = integrated_am
                                best_indices = candidate_indices
                                break

                        else:
                        # if this is the smallest, and is for a contiguous span of time, it's the new one to beat
                            if integrated_am < largest_airmass and contiguous:
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
                else:
                    m_start = packable(goodtime, tgt)
                    if asap and (m_start is not None):
                        tgt.starting_index = squeeze(tgt_i, targets, m_start)
                        o.append(tgt)
                        print('Successfully included field %i: %s %.3E' %(tgt_i, tgt.name, 1./tgt.priority/1000))
                        found = True
                        best_indices = np.arange(tgt.starting_index, tgt.starting_index + 4)
                        tgt.scheduled_airmass_array = np.asarray(tgt.raw_airmass_array)[best_indices]
                        tgt.scheduled_time_array = np.asarray(self.local_time_array)[best_indices]
                        break
                    else:
                        print("Can't fit %s. Skipping!" % tgt.name)
                        bad_o.append(tgt)
                        break

        print('\n\nCurrent schedule:')
        o = sorted(o, key = operator.attrgetter('starting_index'))

        for tgt in o:
            fmt = '{name}: {exposures}; {total} min; Priority: {priority}'
            fmt = fmt.format(name=tgt.name, exposures=tgt.exposures,
                total=tgt.total_minutes, priority='%7.4f'%tgt.priority)
            print(fmt)

        self.ephemeris.date = self.utc_begin_night
        lst = self.ephemeris.sidereal_time()
        vunits = (u.hour, u.deg)
        zenith = SkyCoord(str(lst), self.ephemeris.lat, unit=vunits)

        # Find the coordinate with the largest value that's within ~5 hour of
        # zenith at the start of the night
        print('\n\nSorting by right ascension...\n\n')

        ra_limit = zenith.ra.degree - 75.0
        if ra_limit < 0: ra_limit = ra_limit + 360.0

        # Now get target with smallest RA that is larger than ra_limit
        o = sorted(o, key = lambda t: t.coord.ra.degree)
        target = None
        idx = 0
        for i,tgt in enumerate(o):
            if tgt.coord.ra.degree > ra_limit:
                target = tgt
                idx = i
                break

        if idx-1 > 0:
            o = o[idx:]+o[0:idx-1]

        if len(o)>1:
            # Check if we can get a more efficient ordering of targets
            print('\n\nAttempting to optimize schedule ordering...')
            message = 'Iterating'
            bar = progressbar.ProgressBar(maxval=len(o)-1)
            bar.start()
            for i,tgt in enumerate(o):
                if i==len(o)-1:
                    bar.update(i)
                    continue

                bar.update(i)
                swap = self.swap(o, i)

                if (self.compute_integrated_am(self.swap(o, i)) <
                    self.compute_integrated_am(o)):

                    o = copy.deepcopy(self.swap(o, i))
                    idx = i-1

                    while ((self.compute_integrated_am(self.swap(o, idx)) <
                        self.compute_integrated_am(o)) and idx!=0):

                        o = copy.deepcopy(self.swap(o, idx))
                        idx = i-1

            bar.finish()

        # Resort targets
        o = sorted(o, key=lambda tgt: tgt.starting_index)

        print('\n\nOptimized schedule:')
        for tgt in o:
            fmt = '{name:<22}: {exposures:<32}; {total:<5} min; Priority: {priority}'
            fmt = fmt.format(name=tgt.name,
                exposures=", ".join("=".join([_[0],str(int(_[1]))])
                    for _ in tgt.exposures.items()),
                total=str(tgt.total_minutes), priority='%7.4f'%tgt.priority)
            print(fmt)

            # Also get UTC start time
            #tgt.utc_start_time = self.utc_time_array[tgt.starting_index]


        self.compute_night_fill_fraction(o)

        # Get date based on first and last target in schedule
        tgt0 = o[0]
        tgt1 = o[-1]

        date0 = self.get_start_of_observing_target(tgt0)
        date1 = self.get_start_of_observing_target(tgt1)

        d0 = date0.datetime.strftime('%Y-%m-%d')
        d1 = date1.datetime.strftime('%Y-%m-%d')

        print(f'\n\nWriting out schedule: {d0} to {d1}\n\n')

        telescope.write_schedule(self.name, date0, o, output_files=output_files,
            fieldcenters=fieldcenters, pointing=self.pointing)
