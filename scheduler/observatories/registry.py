"""
Registry of observatories built from site data and telescope implementations.

Sites with telescope drivers (LCO, Lick, Thacher, CTIO) get full Observatory
instances. Other sites are in sites.SITE_DATA for reference; get_observatory()
can build an Observatory with no telescopes for any known site.
"""

from .base import Observatory
from .sites import SITE_DATA, get_site
from .telescopes import Swope, Nickel, Thacher, T80S, Blanco

# Site key -> { telescope_key: TelescopeClass }
TELESCOPE_REGISTRY = {
    "LCO": {"Swope": Swope},
    "Lick": {"Nickel": Nickel},
    "Thacher": {"Thacher": Thacher},
    "CTIO": {"T80S": T80S, "Blanco": Blanco},
}


def _build_observatory(site_key, site_info, telescope_classes):
    """
    Build one Observatory instance from site info and telescope classes.

    Parameters
    ----------
    site_key : str
        Observatory key (e.g. 'LCO').
    site_info : dict
        Must have lon, lat, elevation, horizon, utc_offset, utc_offset_name.
    telescope_classes : dict
        Mapping telescope_name -> Telescope subclass (instances are created).

    Returns
    -------
    Observatory
        Configured observatory.
    """
    telescopes = {name: cls() for name, cls in telescope_classes.items()}
    return Observatory(
        name=site_key,
        lon=site_info["lon"],
        lat=site_info["lat"],
        elevation=site_info["elevation"],
        horizon=site_info["horizon"],
        telescopes=telescopes,
        utc_offset=site_info["utc_offset"],
        utc_offset_name=site_info["utc_offset_name"],
    )


def _build_all():
    """Build observatories dict for sites that have telescope implementations."""
    result = {}
    for site_key, site_info in SITE_DATA.items():
        classes = TELESCOPE_REGISTRY.get(site_key)
        if not classes:
            continue
        result[site_key] = _build_observatory(site_key, site_info, classes)
    return result


# Main dict used by run.py and options: LCO, Lick, Thacher, CTIO
observatories = _build_all()


def get_observatory(site_key):
    """
    Return the Observatory instance for a site key.

    Pre-built observatories (LCO, Lick, Thacher, CTIO) are returned from the
    registry. For any other key in SITE_DATA, builds an Observatory with
    no telescopes.

    Parameters
    ----------
    site_key : str
        Key from SITE_DATA (e.g. 'LCO', 'McDonald').

    Returns
    -------
    Observatory
        Observatory instance.

    Raises
    ------
    KeyError
        If site_key is not in SITE_DATA or observatories.
    """
    if site_key in observatories:
        return observatories[site_key]
    if site_key in SITE_DATA:
        return _build_observatory(site_key, get_site(site_key), {})
    raise KeyError("Unknown observatory: %r" % site_key)


def list_observatories(include_sites_without_telescopes=False):
    """
    Return list of observatory (site) keys.

    Parameters
    ----------
    include_sites_without_telescopes : bool, optional
        If True, return all keys from SITE_DATA. Otherwise only keys that
        have at least one telescope implementation. Default False.

    Returns
    -------
    list of str
        Sorted site keys.
    """
    if include_sites_without_telescopes:
        return sorted(SITE_DATA.keys())
    return sorted(observatories.keys())
