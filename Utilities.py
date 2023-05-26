from datetime import tzinfo, timedelta, datetime
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.time import Time
from astropy.table import unique, Column
import csv, requests, sys, numpy as np, copy

from astroquery.mast import Catalogs

from requests.auth import HTTPBasicAuth

target_table_names = ('name', 'ra', 'dec', 'priority', 'date', 'mag', 'type')
target_table_row   = [['X' * 40], [0.], [0.], [0.0],
    [Time('2019-01-01T00:00:00')], [0.0], ['X' * 40]]
max_length = 20000

class UTC_Offset(tzinfo):

    # Offset assumed to be hours
    def __init__(self, offset=0, name=None):
        self.offset = timedelta(seconds=offset*3600)
        self.name = name or self.__class__.__name__

    def utcoffset(self, dt):
        return self.offset

    def tzname(self, dt):
        return self.name

    def dst(self, dt):
        return timedelta(0)

# Chile observes Chile Summer Time (CLST) from 1/1/2017 - 5/13/2017 => UTC-3
# Chile observes Chile Standard Time (CLT) from 5/13/2017 - 8/12/2017 => UTC-4
# Chile observes Chile Summer Time (CLST) from 8/13/2017 - 12/31/2017 => UTC-3
lco_clst_utc_offset = -3 # hours
lco_clt_utc_offset = -4 # hours

# California observes Pacific Standard Time (PST) from 1/1/2017 - 3/12/2017 => UTC-8
# California observes Pacific Daylight Time (PDT) from 3/12/2017 - 11/5/2017 => UTC-7
# California observes Pacific Standard Time (PST) from 11/5/2017 - 12/31/2017 => UTC-8
lick_pst_utc_offset = -8 # hours
lick_pdt_utc_offset = -7 # hours

keck_offset=-10

def is_number(val):
    try:
        val = float(val)
        return(True)
    except ValueError:
        return(False)

def parse_coord(ra, dec):
    if (not (is_number(ra) and is_number(dec)) and
        (':' not in str(ra) and ':' not in str(dec))):
        return(None)

    if (':' in str(ra) and ':' in str(dec)):
        # Input RA/DEC are sexagesimal
        unit = (u.hourangle, u.deg)
    else:
        unit = (u.deg, u.deg)

    try:
        coord = SkyCoord(ra, dec, frame='icrs', unit=unit)
        return(coord)
    except ValueError:
        return(None)


# Get a blank target table
def blank_target_table():
    tab = Table(target_table_row, names=target_table_names)
    return(tab[:0].copy())

# file_name assumed to be CSV with headers...
# input is file name and an optional value gw.  gw represents the preferred
# distance modulus for a gravitational wave event (DISTMEAN from healpix),
# which will be used to calculate the exposure time required
def get_targets(file_name, gw=None, target_mag=-17.0, obstype='',
    priority=1.0, username='', password=''):

    if 'https://' in file_name or 'http://' in file_name:
        # Assume URL
        url = copy.copy(file_name)
        auth = HTTPBasicAuth(username, password)
        response = requests.get(url, auth=auth)
        if response.status_code==200:
            file_name = response.text
        else:
            raise Exception(f'ERROR: could not download targets from {url}.\n'+\
                'Did you remember to pass authentication (--username/--password)?')

    data_table = ascii.read(file_name)

    # Sanitize columns
    for key in data_table.keys():
        newkey = key.lower().replace(' ','_')
        data_table.rename_column(key, newkey)

    for key in data_table.keys():
        if key in ['fieldra','r.a.','right_ascension']:
            data_table.rename_column(key, 'ra')
        if key in ['fielddec','declination','dec.']:
            data_table.rename_column(key, 'dec')
        if key in ['fieldname','object','field_name']:
            data_table.rename_column(key, 'name')
        if key in ['prob']:
            data_table.rename_column(key, 'priority')
        if key in ['lum_dist','distance']
            # Assume distance in Mpc
            distance = data_table[key]
            if 'dm' not in data_table.keys():
                dm = 5. * np.log10(data_table[key].data.astype('float')) + 25.
                data_table.add_column(Column(dm, name='dm'))

    if len(data_table) > max_length:
        data_table.sort('priority')
        data_table.reverse()
        data_table = data_table[:max_length]


    table = blank_target_table()

    priority_data = np.log10(data_table['priority'])
    for row in data_table:
        new_row = {}
        if gw: 
            new_row['type'] = 'GW'
            if 'dm' in row.colnames:
                new_row['mag'] = target_mag+row['dm']
            else:
                new_row['mag'] = target_mag+float(gw)
            # Check if A_lambda is in table
            if 'a_lambda' in data_table.keys():
                # Adjust the magnitudes to account for a_lambda
                new_row['mag'] = new_row['mag'] + row['a_lambda']
        elif 'type' in row.colnames:
            new_row['type'] = row['type']
        else:
            new_row['type'] = 'SN'

        if 'mag' not in new_row.keys():
            new_row['mag'] = 19.0

        new_row['date'] = Time(datetime.now())

        coord = parse_coord(row['ra'], row['dec'])
        new_row['ra'] = coord.ra.degree
        new_row['dec'] = coord.dec.degree
        new_row['name'] = row['name']
        new_row['priority'] = row['priority']

        # Reprioritize by approximate priority
        new_row['priority'] = np.max(priority_data)+1.0-np.log10(row['priority'])

        table.add_row(new_row)

    return(table)

# Download a list of targets from a website that can provides a parseable table
# by ascii.read.  Create priorities based on info from each target.
def download_targets(url):

        # Use requests to get list of targets
        try:
            data = requests.get(url, timeout=30)
        except:
            error = f'ERROR: could not get a response from {url}.  Exiting...'
            print(error)
            sys.exit()

        # Format into a table with the same names as standard target file
        table = ascii.read(data.text)

        for key in table.keys():
            if 'name' in key:
                table.rename_column(key, 'name')
                break

        # Only want to compare to BVgri magnitudes
        match = [(f.startswith('B') or f.startswith('V')
            or f.startswith('g') or f.startswith('r') or f.startswith('i'))
            for f in table['filter']]
        table = table[match]

        # We can get duplicte targets from the queries, weed these by iterating
        # through targets by name and checking the name
        newtable = table[:0].copy()
        for name in list(set(table['name'])):
            # Find the row with the most recent
            subtable = table[table['name']==name]
            idx = np.argmax([Time(t).mjd for t in subtable['obs_date']])

            newtable.add_row(subtable[idx])

        # Start by generating a blank table
        targets = blank_target_table()

        # Now
        now = Time(datetime.now())

        # Add targets one by one from
        for row in newtable:
            # Get rough approximation of target priority
            pri = 1.0
            delta_t_days = Time(datetime.now()).mjd - Time(row['obs_date']).mjd
            effective_mag = row['Recent mag'] + 0.03 * delta_t_days

            # Prioritize interesting classes of transients
            if row['spec_class'] in ['SN IIn','SN Ib','SN Ic','SN IIb','LBV',
                'SN Ibn','SN Ia-91T-like','SN Ia-91bg-like']:
                pri += 1.0

            # Prioritize bright things
            delta_t_days = now.mjd - Time(row['obs_date']).mjd
            effective_mag = row['Recent mag'] + 0.03 * delta_t_days
            pri += 17.0-effective_mag

            # Prioritize follow up requested
            if row['status_id'] == 2:
                pri += 40.0
            if row['status_id'] == 4:
                pri += 20.0

            try:
                if (Time(datetime.utcnow()).mjd - Time(row['disc_date']).mjd) < 50:
                    pri += 2.0
            except AttributeError:
                continue

            pri = round(pri, 5)

            time    = Time(row['obs_date'])
            add_row = [row['name'], row['ra'],row['dec'], pri, now,
                effective_mag, 'SN']
            targets.add_row(add_row)

        # Rearrange priority for descending order
        max_priority = np.max(targets['priority'])+1
        targets['priority'] = max_priority - targets['priority']

        return(targets)

def download_ps1_catalog(target, Mmax=18.0, radius=0.05):
        
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
