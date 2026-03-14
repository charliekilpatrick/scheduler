"""
Observatory site data: geographic and timezone information.

Each site is a dict with: name, lon, lat, elevation, horizon, utc_offset, utc_offset_name.
lon/lat are strings for ephem; elevation in meters; horizon in degrees (sun);
utc_offset in hours (e.g. -8 for PST); utc_offset_name is the timezone abbreviation.

Sources: observatory web sites, Wikipedia, ESO obs_sites, IAU.
"""

# Default horizon (degrees below horizon for twilight) used when not specified
DEFAULT_HORIZON = "-12"

# Site key -> site info. Keys are used in --obstele (e.g. LCO:Swope).
SITE_DATA = {
    # --- Existing scheduler sites (with telescope implementations) ---
    "LCO": {
        "name": "Las Campanas Observatory",
        "lon": "-70.6915",
        "lat": "-29.0182",
        "elevation": 2402,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -4,  # Chile Standard (CLT); Chile Summer is -3
        "utc_offset_name": "CLT",
    },
    "Lick": {
        "name": "Lick Observatory",
        "lon": "-121.6429",
        "lat": "37.3414",
        "elevation": 1283,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -8,  # PST; PDT is -7
        "utc_offset_name": "PST",
    },
    "Thacher": {
        "name": "Thacher Observatory",
        "lon": "-121.6431",
        "lat": "34.46479",
        "elevation": 630.0,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -8,
        "utc_offset_name": "PST",
    },
    "CTIO": {
        "name": "Cerro Tololo Inter-American Observatory",
        "lon": "-70.8035",
        "lat": "-30.1732",
        "elevation": 2207.0,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -4,
        "utc_offset_name": "CLT",
    },
    # --- Additional public observatories (site data only; no telescope drivers yet) ---
    "McDonald": {
        "name": "McDonald Observatory",
        "lon": "-104.022",
        "lat": "30.6714",
        "elevation": 2077,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -6,  # CST; CDT is -5
        "utc_offset_name": "CST",
    },
    "SidingSpring": {
        "name": "Siding Spring Observatory",
        "lon": "149.0587",
        "lat": "-31.2707",
        "elevation": 1165,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": 10,  # AEST; AEDT is 11
        "utc_offset_name": "AEST",
    },
    "SAAO": {
        "name": "South African Astronomical Observatory (Sutherland)",
        "lon": "20.8107",
        "lat": "-32.3794",
        "elevation": 1798,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": 2,
        "utc_offset_name": "SAST",
    },
    "KittPeak": {
        "name": "Kitt Peak National Observatory",
        "lon": "-111.5967",
        "lat": "31.9583",
        "elevation": 2096,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -7,  # MST, no DST
        "utc_offset_name": "MST",
    },
    "MaunaKea": {
        "name": "Mauna Kea Observatories",
        "lon": "-155.47",
        "lat": "19.82",
        "elevation": 4205,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -10,  # HST
        "utc_offset_name": "HST",
    },
    "Paranal": {
        "name": "ESO Paranal Observatory",
        "lon": "-70.4030",
        "lat": "-24.6275",
        "elevation": 2635,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -4,
        "utc_offset_name": "CLT",
    },
    "LaSilla": {
        "name": "ESO La Silla Observatory",
        "lon": "-70.7333",
        "lat": "-29.2567",
        "elevation": 2400,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -4,
        "utc_offset_name": "CLT",
    },
    "Palomar": {
        "name": "Palomar Observatory",
        "lon": "-116.865",
        "lat": "33.356",
        "elevation": 1706,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": -8,
        "utc_offset_name": "PST",
    },
    "RoqueMuchachos": {
        "name": "Roque de los Muchachos Observatory",
        "lon": "-17.894",
        "lat": "28.758",
        "elevation": 2396,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": 0,  # WET; WEST is 1
        "utc_offset_name": "WET",
    },
    "Teide": {
        "name": "Teide Observatory",
        "lon": "-16.5094",
        "lat": "28.3000",
        "elevation": 2390,
        "horizon": DEFAULT_HORIZON,
        "utc_offset": 0,
        "utc_offset_name": "WET",
    },
}


def get_site(site_key: str) -> dict:
    """
    Return site info dict for a given site key.

    Parameters
    ----------
    site_key : str
        Key from SITE_DATA (e.g. 'LCO', 'McDonald').

    Returns
    -------
    dict
        Site info with name, lon, lat, elevation, horizon, utc_offset, utc_offset_name.

    Raises
    ------
    KeyError
        If site_key is not in SITE_DATA.
    """
    return dict(SITE_DATA[site_key])


def list_site_keys():
    """
    Return sorted list of all site keys in SITE_DATA.

    Returns
    -------
    list of str
        All keys (e.g. 'LCO', 'Lick', 'McDonald').
    """
    return sorted(SITE_DATA.keys())
