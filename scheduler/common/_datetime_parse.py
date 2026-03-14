"""
Parse date/datetime strings with common formats (stdlib only).

Used in place of python-dateutil for known formats.
"""
from datetime import datetime


def parse_datetime_ymd_hm(s: str) -> datetime:
    """
    Parse a datetime string 'YYYY-MM-DD HH:MM' or 'YYYY-MM-DD 12:00'.

    Parameters
    ----------
    s : str
        String like '2024-06-15 12:00' or '2024-06-15 22:30'.

    Returns
    -------
    datetime
        Naive datetime.

    Raises
    ------
    ValueError
        If string does not match format.
    """
    return datetime.strptime(s.strip(), "%Y-%m-%d %H:%M")


def parse_date_flexible(s: str) -> datetime:
    """
    Parse a date string (YYYY-MM-DD or YYYYMMDD) to datetime at midnight.

    Parameters
    ----------
    s : str
        Date string in YYYY-MM-DD or YYYYMMDD format.

    Returns
    -------
    datetime
        Naive datetime at 00:00:00.

    Raises
    ------
    ValueError
        If string does not match either format.
    """
    s = s.strip().replace("-", "")
    if len(s) == 8 and s.isdigit():
        return datetime.strptime(s[:4] + "-" + s[4:6] + "-" + s[6:8], "%Y-%m-%d")
    return datetime.strptime(s, "%Y-%m-%d")


def parse_month_day_year(month: str, day: str, year: str) -> datetime:
    """
    Parse 'Month DD, YYYY' or 'Mon DD, YYYY' to datetime at midnight.

    Parameters
    ----------
    month : str
        Full or abbreviated month name (e.g. 'January', 'Jan').
    day : str
        Day of month.
    year : str
        Four-digit year.

    Returns
    -------
    datetime
        Naive datetime at 00:00:00.

    Raises
    ------
    ValueError
        If format is not recognized.
    """
    s = f"{month} {day}, {year}"
    for fmt in ("%B %d, %Y", "%b %d, %Y"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    raise ValueError(f"Cannot parse date: {s}")
