"""
Telescope scheduler package.

Main entry point: scheduler.run.main() or the create-schedule console script.
"""
__version__ = "0.3.0"

from scheduler.common import (
    constants,
    Target,
    TargetType,
    handle_target_type,
    get_targets,
    parse_coord,
)
from scheduler.observatories import Observatory, Telescope

__all__ = [
    "__version__",
    "constants",
    "Observatory",
    "Telescope",
    "Target",
    "TargetType",
    "handle_target_type",
    "get_targets",
    "parse_coord",
]
