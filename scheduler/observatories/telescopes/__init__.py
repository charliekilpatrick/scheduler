"""Telescope implementations for scheduler-supported observatories."""

from .swope import Swope
from .nickel import Nickel
from .thacher import Thacher
from .t80s import T80S
from .blanco import Blanco

__all__ = ["Swope", "Nickel", "Thacher", "T80S", "Blanco"]
