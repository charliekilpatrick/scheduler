"""
GW (gravitational-wave) tiles conversion and schedule comparison.

Generic functions for any observatory/telescope: convert GW tiles to
scheduler input CSV, filter tiles by schedule output, and compute
schedule-vs-tiles coincidence statistics.
"""
from astropy.io import ascii

from scheduler.common.logging import get_logger
from scheduler.io import SCHEDULER_CSV_HEADER, read_scheduled_names_from_csv

logger = get_logger(__name__)


def _iter_tile_lines(path):
    """Yield (stripped_line, parts) for non-empty, non-comment lines."""
    with open(path) as f:
        for line in f:
            line = line.replace("\n", "").replace(" ", "")
            if not line or line.startswith("F") or line.startswith("#"):
                continue
            yield line, line.split(",")


def convert_gw_to_scheduler(
    input_path,
    output_path,
    *,
    invert_priority=False,
    add_string_regular="1,2019-04-12,12,GW",
    add_string_inverted="2019-04-12,12,GW",
):
    """
    Convert a GW tiles file to scheduler input CSV format.

    Works for any observatory/telescope; output is the standard scheduler
    input format.

    Parameters
    ----------
    input_path : str or path-like
        Path to GW tiles input file.
    output_path : str or path-like
        Path for output CSV (regular mode appends .csv if missing).
    invert_priority : bool, optional
        If True, use inverted-priority format (5 columns, priority from
        input column 7). Default False.
    add_string_regular : str, optional
        String for Priority column in regular (4-column) mode.
    add_string_inverted : str, optional
        String for column 4 in inverted mode; column 3 = 1/input_priority/1000.

    Returns
    -------
    None
    """
    if invert_priority:
        add_string = add_string_inverted
        with open(output_path, "w") as fout:
            print(SCHEDULER_CSV_HEADER, file=fout)
            for line, parts in _iter_tile_lines(input_path):
                if len(parts) < 7:
                    logger.warning(
                        "Line has %s columns, need 7 for invert: %s...",
                        len(parts),
                        line[:60],
                    )
                    continue
                try:
                    inv_priority = 1.0 / float(parts[6]) / 1000.0
                except (ValueError, ZeroDivisionError):
                    inv_priority = 0.0
                parts[3] = str(inv_priority)
                parts[4] = add_string
                print(",".join(parts[:5]), file=fout)
    else:
        if not str(output_path).endswith(".csv"):
            output_path = str(output_path) + ".csv"
        add_string = add_string_regular
        with open(output_path, "w") as fout:
            print(SCHEDULER_CSV_HEADER, file=fout)
            for line, parts in _iter_tile_lines(input_path):
                if len(parts) < 4:
                    parts.extend([""] * (4 - len(parts)))
                parts[3] = add_string
                print(",".join(parts[:4]), file=fout)


def filter_tiles_by_schedule(scheduler_csv_path, tiles_path, output_path):
    """
    Filter a GW tiles file to only lines matching scheduled object names.

    Reads scheduled names from the scheduler CSV (first column), then
    keeps only lines in the tiles file that contain one of those names.
    Generic for any observatory/telescope schedule CSV.

    Parameters
    ----------
    scheduler_csv_path : str or path-like
        Path to scheduler output CSV (any telescope).
    tiles_path : str or path-like
        Path to full GW tiles file.
    output_path : str or path-like
        Path for filtered output file.

    Returns
    -------
    int
        Number of scheduled names (for logging).
    """
    scheduled = read_scheduled_names_from_csv(scheduler_csv_path)
    logger.info("Number of tiles scheduled: %s", len(scheduled))
    with open(tiles_path) as f:
        lines = f.readlines()
    to_keep = []
    for line in lines:
        for name in scheduled:
            if name in line:
                to_keep.append(line)
                break
    with open(output_path, "w") as f:
        for line in to_keep:
            f.write(line)
    return len(scheduled)


def schedule_coincidence_stats(schedule_csv_path, tiles_path):
    """
    Compare a schedule CSV to a tiles file; return coincidence and priority sum.

    Tiles file must be CSV with columns FieldName and Priority. Generic for
    any observatory/telescope schedule CSV.

    Parameters
    ----------
    schedule_csv_path : str or path-like
        Path to scheduler output CSV.
    tiles_path : str or path-like
        Path to GW tiles CSV with FieldName and Priority columns.

    Returns
    -------
    dict
        Keys: n_scheduled, n_coincidence, sum_priority, probability_100 (100 * sum_priority).

    Raises
    ------
    ValueError
        If tiles file lacks FieldName or Priority columns.
    """
    scheduled = read_scheduled_names_from_csv(schedule_csv_path)
    scheduled = [n for n in scheduled if n]
    tiles = ascii.read(tiles_path, data_start=0, delimiter=",")
    if "FieldName" not in tiles.colnames or "Priority" not in tiles.colnames:
        raise ValueError(
            f"Tiles file must have 'FieldName' and 'Priority' columns. "
            f"Found: {tiles.colnames}"
        )
    coincidence = []
    ps = []
    for name in scheduled:
        for j, val in enumerate(tiles["FieldName"]):
            if name == val:
                coincidence.append(name)
                ps.append(tiles["Priority"][j])
                break
    sum_p = sum(ps)
    return {
        "n_scheduled": len(scheduled),
        "n_coincidence": len(coincidence),
        "sum_priority": sum_p,
        "probability_100": 100.0 * sum_p,
    }


def main_convert():
    """CLI: convert GW tiles file to scheduler input CSV."""
    import argparse
    ap = argparse.ArgumentParser(
        description="Convert GW tiles file to scheduler input CSV.",
        epilog="Default: 4 columns. With --invert-priority: 5 columns, priority from input column 7.",
    )
    ap.add_argument("input_file", help="GW tiles input file")
    ap.add_argument("output_file", help="Output CSV path (regular mode appends .csv if missing)")
    ap.add_argument("--invert-priority", action="store_true",
                    help="Use inverted-priority format (5 cols, priority from input col 7)")
    ap.add_argument("--add-string", default=None,
                    help="Override default addstring for Priority/DiscDates column(s)")
    args = ap.parse_args()
    add_regular = args.add_string or "1,2019-04-12,12,GW"
    add_inverted = args.add_string or "2019-04-12,12,GW"
    convert_gw_to_scheduler(
        args.input_file,
        args.output_file,
        invert_priority=args.invert_priority,
        add_string_regular=add_regular,
        add_string_inverted=add_inverted,
    )


def main_filter():
    """CLI: filter tiles file to scheduled tiles using scheduler CSV."""
    import argparse
    ap = argparse.ArgumentParser(description="Filter tiles file to scheduled tiles using scheduler CSV.")
    ap.add_argument("-s", "--scheduler-file", required=True, help="Scheduler output CSV")
    ap.add_argument("-t", "--tiles-file", required=True, help="Full GW tiles file")
    ap.add_argument("-c", "--fc-file", "--output", dest="output", required=True, help="Output filtered file")
    args = ap.parse_args()
    filter_tiles_by_schedule(args.scheduler_file, args.tiles_file, args.output)
    logger.info("Scheduled names: %s", read_scheduled_names_from_csv(args.scheduler_file))


def main_probs():
    """CLI: compare schedule CSV to tiles file; print coincidence and priority sum."""
    import argparse
    ap = argparse.ArgumentParser(description="Compare schedule CSV to tiles file; report coincidence and priority sum.")
    ap.add_argument("-s", "--schedule-file", required=True, help="Scheduler output CSV")
    ap.add_argument("-f", "--tiles-file", required=True, dest="tiles_file", help="Tiles CSV with FieldName and Priority")
    args = ap.parse_args()
    stats = schedule_coincidence_stats(args.schedule_file, args.tiles_file)
    print("Number of tiles scheduled:", stats["n_scheduled"])
    print("Coincidence (overlap with tiles file):", stats["n_coincidence"])
    print("Sum of Priorities:", stats["sum_priority"])
    print("Probability (100 * sum):", stats["probability_100"])
