"""Target loading, coordinate parsing, PS1 catalog, timezone helpers, RTML XML."""

from datetime import tzinfo, timedelta, datetime
import copy
import itertools
import numpy as np

from astropy import units as u
from astropy.table import Table, unique, Column
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.time import Time
from astroquery.vizier import Vizier
from requests.auth import HTTPBasicAuth
import requests

from .logging import get_logger

logger = get_logger(__name__)

target_table_names = ('name', 'ra', 'dec', 'priority', 'date', 'mag', 'type', 'orig_priority')
target_table_row = [['X' * 40], [0.], [0.], [0.0], [Time('2019-01-01T00:00:00')], [0.0], ['X' * 40], [0.]]
max_length = 20000


class UTC_Offset(tzinfo):
    """Fixed UTC offset timezone (no DST)."""

    def __init__(self, offset=0, name=None):
        self.offset = timedelta(seconds=offset * 3600)
        self.name = name or self.__class__.__name__

    def utcoffset(self, dt):
        return self.offset

    def tzname(self, dt):
        return self.name

    def dst(self, dt):
        return timedelta(0)


def is_number(val):
    """True if float(val) succeeds."""
    try:
        float(val)
        return True
    except (ValueError, TypeError):
        return False


def parse_coord(ra, dec):
    """Parse RA/Dec (sexagesimal or degrees) to SkyCoord; None on failure."""
    if (not (is_number(ra) and is_number(dec)) and
            (':' not in str(ra) and ':' not in str(dec))):
        return None
    unit = (u.hourangle, u.deg) if (':' in str(ra) and ':' in str(dec)) else (u.deg, u.deg)
    try:
        return SkyCoord(ra, dec, frame='icrs', unit=unit)
    except ValueError:
        return None


def blank_target_table():
    """Empty table with standard target column names."""
    return Table(target_table_row, names=target_table_names)[:0].copy()


def get_targets(file_name, gw=None, target_mag=-17.0, obstype='', priority=1.0,
                username='', password='', newfirm=False):
    """Load target list from file or URL; returns Table."""
    if 'https://' in file_name or 'http://' in file_name:
        auth = HTTPBasicAuth(username, password)
        response = requests.get(file_name, auth=auth)
        if response.status_code != 200:
            logger.error("Could not download targets from %s", file_name)
            raise RuntimeError("Could not download. Check --username/--password?")
        file_name = response.text
    data_table = ascii.read(file_name)
    for key in list(data_table.keys()):
        newkey = key.lower().replace(' ', '_')
        if key != newkey and newkey not in data_table.keys():
            data_table.rename_column(key, newkey)
    for key in data_table.keys():
        if 'name' in key and 'name' not in data_table.keys():
            data_table.rename_column(key, 'name')
            break
        if 'metric' in key and 'priority' not in data_table.keys():
            data_table['priority'] = data_table['metric']
    data_table = unique(data_table, keys='name')
    for key in data_table.keys():
        if key in ['fieldra', 'r.a.', 'right_ascension']:
            data_table.rename_column(key, 'ra')
        if key in ['fielddec', 'declination', 'dec.', '_dec']:
            data_table.rename_column(key, 'dec')
        if key in ['fieldname', 'object', 'field_name', 'objname']:
            data_table.rename_column(key, 'name')
        if key in ['percentile', 'p_3d', 'galaxyprob']:
            data_table.rename_column(key, 'priority')
        if key in ['recent_mag']:
            data_table.rename_column(key, 'mag')
    if 'filter' in data_table.keys():
        mask = np.array([r['filter'] in ['r-ZTF', 'g-ZTF', 'r', 'i', 'B', 'V', 'g', 'z', 'J', 'H', 'Ks'] for r in data_table])
        data_table = data_table[mask]
    if 'obs_date' in data_table.keys():
        newtable = data_table[:0].copy()
        for name in list(set(data_table['name'])):
            subtable = data_table[data_table['name'] == name]
            idx = np.argmax([Time(t).mjd for t in subtable['obs_date']])
            newtable.add_row(subtable[idx])
        data_table = copy.copy(newtable)
    if 'priority' in data_table.keys():
        data_table.add_column(Column(data_table['priority'], name='orig_priority'))
    if 'priority' not in data_table.keys():
        priority_list = []
        now = Time(datetime.now())
        for row in data_table:
            pri = 1.0
            if row.get('spec_class') in ['SN IIn', 'SN Ib', 'SN Ic', 'SN IIb', 'LBV', 'SN Ibn', 'SN Ia-91T-like', 'SN Ia-91bg-like']:
                pri += 200.0
            delta_t_days = now.mjd - Time(row['obs_date']).mjd
            effective_mag = row['mag'] + 0.03 * delta_t_days
            pri += 20.0 * (17.0 - effective_mag)
            if row.get('status_id') == 2:
                pri += 400.0
            if row.get('status_id') == 4:
                pri += 200.0
            priority_list.append(round(pri, 5))
        priority_arr = np.array(priority_list) - np.min(priority_list) + 1.0
        data_table.add_column(Column(priority_arr, name='priority'))
        data_table = data_table[~np.isnan(data_table['priority'])]
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
            new_row['mag'] = target_mag + float(gw)
            if 'a_lambda' in data_table.keys():
                new_row['mag'] = new_row['mag'] + row['a_lambda']
        elif 'type' in row.colnames:
            new_row['type'] = row['type']
        elif newfirm:
            new_row['type'] = 'NEWFIRM'
        else:
            new_row['type'] = 'SN'
        new_row['mag'] = row['mag'] if 'mag' in row.colnames else 19.0
        new_row['date'] = Time(datetime.now())
        coord = parse_coord(row['ra'], row['dec'])
        new_row['ra'] = coord.ra.degree
        new_row['dec'] = coord.dec.degree
        new_row['name'] = row['name']
        new_row['orig_priority'] = row['orig_priority']
        new_row['priority'] = np.max(priority_data) + 1.0 - np.log10(row['priority'])
        table.add_row(new_row)
    if 'orig_priority' not in table.keys():
        table.add_column(Column([None] * len(table), name='orig_priority'))
    return table


def download_ps1_catalog(target, Mmax=18.0, radius=0.05):
    """Query Pan-STARRS1 around target; returns Table."""
    coord = target.coord
    vquery = Vizier(columns=['RAJ2000', 'DEJ2000', 'gmag', 'rmag', 'imag', 'zmag', 'ymag'],
                    column_filters={'gmag': ('<%f' % Mmax)}, row_limit=100000)
    resp = vquery.query_region(coord, width='%fd' % (4 * radius), catalog='II/349/ps1')
    if len(resp) == 0:
        return Table([[0.], [0.], [0.], [0.], [0.], [0.], [0.]],
                     names=('raStack', 'decStack', 'gMeanPSFMag', 'rMeanPSFMag',
                            'iMeanPSFMag', 'zMeanPSFMag', 'yMeanPSFMag')).copy()[:0]
    tbdata = resp[0]
    tbdata.rename_column('RAJ2000', 'raStack')
    tbdata.rename_column('DEJ2000', 'decStack')
    tbdata.rename_column('gmag', 'gMeanPSFMag')
    tbdata.rename_column('rmag', 'rMeanPSFMag')
    tbdata.rename_column('imag', 'iMeanPSFMag')
    tbdata.rename_column('zmag', 'zMeanPSFMag')
    tbdata.rename_column('ymag', 'yMeanPSFMag')
    coords = SkyCoord(tbdata['raStack'], tbdata['decStack'], unit='deg')
    return tbdata[coords.separation(coord) < radius * u.degree]


def dicttoxml(file, dictionary, depth=0, name_name='name', order=None):
    """Write nested dict to XML. order: key order; 'Request' expands to all RequestN keys."""
    if order is None:
        order = []
    if not order:
        it = dictionary
    else:
        reqs = [k for k in dictionary if 'Request' in k]
        it = list(itertools.chain.from_iterable([reqs if x == 'Request' else [x] for x in order]))
    for item in it:
        if item not in dictionary:
            continue
        if item == 'xml':
            ver, enc = dictionary[item]['version'], dictionary[item]['encoding']
            file.write('<?xml version="%s" encoding="%s"?>\n' % (ver, enc))
        else:
            try:
                file.write('\t' * depth + '<%s>%s</%s>\n' % (item, dictionary[item][name_name], item))
            except Exception:
                item_name = 'Request' if 'Request' in item else item
                append = ''
                if isinstance(dictionary[item], dict):
                    for k, v in dictionary[item].items():
                        if not isinstance(v, dict):
                            append += ' %s="%s"' % (k, v) if isinstance(v, str) else ' %s=%s' % (k, v)
                file.write('\t' * depth + '<%s%s>\n' % (item_name, append))
                dicttoxml(file, dictionary[item], depth=depth + 1, name_name=name_name, order=order)
                file.write('\t' * depth + '</%s>\n' % item_name)
