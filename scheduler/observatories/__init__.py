"""
Observatories package: site data, base classes, telescope implementations, and registry.

Use the registry for runtime:
  - observatories: dict of site_key -> Observatory (LCO, Lick, Thacher, CTIO)
  - get_observatory(site_key): return Observatory, building from SITE_DATA if needed
  - list_observatories(): list of site keys with telescope implementations

To add a new observatory: add an entry to sites.SITE_DATA.
To add a new telescope driver: add a class in telescopes/ and register it in registry.TELESCOPE_REGISTRY.
"""

from .base import Observatory, Telescope
from .sites import SITE_DATA, get_site, list_site_keys
from .registry import (
    observatories,
    get_observatory,
    list_observatories,
    TELESCOPE_REGISTRY,
)
from .telescopes import Swope, Nickel, Thacher, T80S, Blanco

__all__ = [
    "Observatory",
    "Telescope",
    "SITE_DATA",
    "get_site",
    "list_site_keys",
    "observatories",
    "get_observatory",
    "list_observatories",
    "TELESCOPE_REGISTRY",
    "Swope",
    "Nickel",
    "Thacher",
    "T80S",
    "Blanco",
]
