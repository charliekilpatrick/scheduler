from datetime import tzinfo, timedelta, datetime
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.time import Time
from astropy.table import unique, Column
import csv, requests, sys, numpy as np, copy

from astroquery.vizier import Vizier
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
        if 'name' in key and 'name' not in data_table.keys():
            data_table.rename_column(key, 'name')
            break

    for key in data_table.keys():
        if key in ['fieldra','r.a.','right_ascension']:
            data_table.rename_column(key, 'ra')
        if key in ['fielddec','declination','dec.']:
            data_table.rename_column(key, 'dec')
        if key in ['fieldname','object','field_name']:
            data_table.rename_column(key, 'name')
        if key in ['prob']:
            data_table.rename_column(key, 'priority')
        if key in ['recent_mag']:
            data_table.rename_column(key, 'mag')


    if 'filter' in data_table.keys():
        mask = np.array([r['filter'] in ['r-ZTF','g-ZTF','r','i','B','V','g',
                'z'] for r in data_table])
        data_table = data_table[mask]

    # We can get duplicte targets from the queries, weed these by iterating
    # through targets by name and checking the name
    if 'obs_date' in data_table.keys():
        newtable = data_table[:0].copy()
        for name in list(set(data_table['name'])):
            # Find the row with the most recent
            subtable = data_table[data_table['name']==name]
            idx = np.argmax([Time(t).mjd for t in subtable['obs_date']])

            newtable.add_row(subtable[idx])

        data_table = copy.copy(newtable)

    if 'priority' not in data_table.keys():
        # Add the original priority to the table as a new column
        data_table.add_column(Column(data_table['priority'], 
            name='orig_priority'))
        # Add targets one by one from
        priority = []
        now = Time(datetime.now())
        for row in data_table:
            # Get rough approximation of target priority
            pri = 1.0

            # Prioritize interesting classes of transients
            if row['spec_class'] in ['SN IIn','SN Ib','SN Ic','SN IIb','LBV',
                'SN Ibn','SN Ia-91T-like','SN Ia-91bg-like']:
                pri += 200.0

            # Prioritize bright things
            delta_t_days = now.mjd - Time(row['obs_date']).mjd
            effective_mag = row['mag'] + 0.03 * delta_t_days
            pri += 20.*(17.0-effective_mag)

            # Prioritize follow up requested
            if row['status_id'] == 2:
                pri += 400.0
            if row['status_id'] == 4:
                pri += 200.0

            pri = round(pri, 5)

            priority.append(pri)

        priority = np.array(priority)
        priority = priority - np.min(priority) + 1.0

        data_table.add_column(Column(priority, name='priority'))

        mask = ~np.isnan(data_table['priority'])
        data_table = data_table[mask]

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
            new_row['mag'] = target_mag+float(gw)
            # Check if A_lambda is in table
            if 'a_lambda' in data_table.keys():
                # Adjust the magnitudes to account for a_lambda
                new_row['mag'] = new_row['mag'] + row['a_lambda']
        elif 'type' in row.colnames:
            new_row['type'] = row['type']
        else:
            new_row['type'] = 'SN'

        if 'mag' in row.colnames:
            new_row['mag'] = row['mag']

        if 'mag' not in new_row.keys():
            new_row['mag'] = 19.0

        new_row['date'] = Time(datetime.now())

        coord = parse_coord(row['ra'], row['dec'])
        new_row['ra'] = coord.ra.degree
        new_row['dec'] = coord.dec.degree
        new_row['name'] = row['name']
        new_row['priority'] = row['priority']

        # Reprioritize by approximate priority
        new_row['priority'] = np.max(priority_data) + 1.0 - np.log10(row['priority'])

        table.add_row(new_row)

    return(table)

def download_ps1_catalog(target, Mmax=18.0, radius=0.05):

    coord = target.coord

    vquery = Vizier(columns=['RAJ2000', 'DEJ2000','gmag','rmag','imag',
            'zmag','ymag'],column_filters={'gmag':('<%f' % Mmax)},
            row_limit=100000)

    w = ('%fd' % (4*radius))
    cat = 'II/349/ps1'
    resp = vquery.query_region(coord, width=w, catalog=cat)

    if len(resp)==0:
        table = Table([[0.],[0.],[0.],[0.],[0.],[0.],[0.]],
            names=('raStack','decStack','gMeanPSFMag','rMeanPSFMag',
                'iMeanPSFMag','zMeanPSFMag','yMeanPSFMag')).copy()[:0]
        return(table)
    else:
        tbdata = resp[0]

    tbdata.rename_column('RAJ2000', 'raStack')
    tbdata.rename_column('DEJ2000', 'decStack')
    tbdata.rename_column('gmag', 'gMeanPSFMag')
    tbdata.rename_column('rmag', 'rMeanPSFMag')
    tbdata.rename_column('imag', 'iMeanPSFMag')
    tbdata.rename_column('zmag', 'zMeanPSFMag')
    tbdata.rename_column('ymag', 'yMeanPSFMag')

    coords = SkyCoord(tbdata['raStack'], tbdata['decStack'], unit='deg')
    mask = coords.separation(coord) < radius * u.degree

    tbdata = tbdata[mask]

    return(tbdata)
