"""
Check whether a schedule should be generated for a given observatory/telescope and date.

Uses a small registry: Lick:Nickel uses Google Sheets (optional env vars);
other obstele values return True (schedule should be generated) by default.
"""
import os

from scheduler.common._datetime_parse import parse_date_flexible
from scheduler.common.logging import get_logger

logger = get_logger(__name__)

# Observatory:telescope keys that use Sheets-based checks (e.g. Nickel).
SHEETS_OBSTELE = frozenset({"Lick:Nickel"})


def need_to_generate_schedule_for_obstele(obstele, date, **env):
    """
    Return True if a schedule should be generated for this observatory/telescope on date.

    For Lick:Nickel, uses Google Sheets (NICKEL_SCHEDULE_SHEET, NICKEL_OBSERVERS_SHEET,
    NICKEL_TEMPLATE_SHEET, GSHEETS_TOKEN). If env vars are missing or Sheets check
    fails, returns True. For any other obstele, returns True.

    Parameters
    ----------
    obstele : str
        Observatory:telescope key (e.g. "Lick:Nickel", "CTIO:Blanco").
    date : datetime or str
        Observation date. If str, parsed with parse_date_flexible (YYYYMMDD or YYYY-MM-DD).
    **env : optional
        Override env vars: current_sheet, observers_sheet, template_sheet, gsheets_token.
        If not provided, read from NICKEL_SCHEDULE_SHEET, NICKEL_OBSERVERS_SHEET,
        NICKEL_TEMPLATE_SHEET, GSHEETS_TOKEN.

    Returns
    -------
    bool
        True if a schedule should be generated, False if one already exists (Nickel only).
    """
    if obstele not in SHEETS_OBSTELE:
        return True
    if hasattr(date, "strftime"):
        dt = date
    else:
        dt = parse_date_flexible(str(date))
    current_sheet = env.get("current_sheet") or os.environ.get("NICKEL_SCHEDULE_SHEET")
    observers_sheet = env.get("observers_sheet") or os.environ.get("NICKEL_OBSERVERS_SHEET")
    template_sheet = env.get("template_sheet") or os.environ.get("NICKEL_TEMPLATE_SHEET")
    gsheets_token = env.get("gsheets_token") or os.environ.get("GSHEETS_TOKEN")
    for name, val in (
        ("NICKEL_SCHEDULE_SHEET", current_sheet),
        ("NICKEL_OBSERVERS_SHEET", observers_sheet),
        ("NICKEL_TEMPLATE_SHEET", template_sheet),
        ("GSHEETS_TOKEN", gsheets_token),
    ):
        if not val:
            logger.warning("Sheets check skipped: %s not set", name)
            return True
    try:
        from scheduler.common.logs import (
            gsheets_params,
            initiate_gsheet,
            check_if_tel_on_date,
            need_to_generate_schedule,
        )
    except ImportError as e:
        logger.warning("Sheets check skipped: %s", e)
        return True
    params = gsheets_params(
        current_sheet,
        template_sheet,
        gsheets_token,
        observers_sheet=observers_sheet,
    )
    sheet = initiate_gsheet(params["gsheets_token"])
    observer, _ = check_if_tel_on_date(sheet, params["OBSERVERS_SHEET"], dt)
    if not observer:
        return False
    sheet_name = dt.strftime("%Y%m%d")
    if need_to_generate_schedule(sheet, params["CURRENT_SHEET"], sheet_name):
        return True
    return False


def main():
    """
    CLI: print 1 if a schedule should be generated for obstele on date, 0 otherwise.

    For Lick:Nickel requires NICKEL_SCHEDULE_SHEET, NICKEL_OBSERVERS_SHEET,
    NICKEL_TEMPLATE_SHEET, GSHEETS_TOKEN. For other observatory:telescope
    always prints 1.
    """
    import sys
    import argparse
    ap = argparse.ArgumentParser(description="Check if schedule should be generated for obstele on date.")
    ap.add_argument("date", help="Date (YYYYMMDD or YYYY-MM-DD)")
    ap.add_argument("-ot", "--obstele", default="Lick:Nickel", help="Observatory:telescope (default Lick:Nickel)")
    args = ap.parse_args()
    need = need_to_generate_schedule_for_obstele(args.obstele, args.date)
    print(1 if need else 0)
    sys.exit(0)
