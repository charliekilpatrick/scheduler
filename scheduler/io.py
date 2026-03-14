"""
I/O helpers for scheduler CSV inputs and outputs.

Used by GW conversion, tile filtering, and schedule comparison tools.
"""
from scheduler.common.logging import get_logger

logger = get_logger(__name__)

# Standard header for scheduler input CSV (GW-converted and similar).
SCHEDULER_CSV_HEADER = "Name,RA,DEC,Priority,DiscDates,DiscMags,Type"


def read_scheduled_names_from_csv(csv_path):
    """
    Read object names from the first column of a scheduler output CSV.

    Skips the header row and returns non-empty names in order. Works for any
    observatory/telescope schedule CSV that uses this format.

    Parameters
    ----------
    csv_path : str or path-like
        Path to the scheduler output CSV file.

    Returns
    -------
    list of str
        Non-empty object names from the first column, in order.
    """
    names = []
    with open(csv_path) as f:
        next(f)  # skip header
        for line in f:
            name = line.split(",")[0].strip()
            if name:
                names.append(name)
    return names
