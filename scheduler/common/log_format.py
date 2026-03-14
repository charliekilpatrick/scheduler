"""Schedule log format: CSV, pickle, and Google Sheets (same schema)."""
import csv
import gzip
import pickle
from typing import Any

from scheduler.common.logging import get_logger

logger = get_logger(__name__)

# Column order and header names for CSV / Google Sheets display
LOG_HEADER = [
    "date",
    "observatory",
    "telescope",
    "observer",
    "target_name",
    "ra",
    "dec",
    "filter",
    "exptime_sec",
    "priority",
    "notes",
]


def validate_log_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """Return entry with all LOG_HEADER keys present; fill missing with ''."""
    out = {}
    for key in LOG_HEADER:
        out[key] = entry.get(key, "")
    return out


def log_entries_to_rows(entries: list[dict[str, Any]]) -> list[list[Any]]:
    """Return [LOG_HEADER] + data rows for CSV or Sheets."""
    rows = [LOG_HEADER]
    for e in entries:
        validated = validate_log_entry(e)
        rows.append([validated[k] for k in LOG_HEADER])
    return rows


def write_log_csv(entries: list[dict[str, Any]], path: str) -> None:
    """Write log entries to CSV (Excel/Sheets-friendly)."""
    rows = log_entries_to_rows(entries)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerows(rows)
    logger.info("Wrote %s log rows to %s", len(entries), path)


def write_log_pickle(entries: list[dict[str, Any]], path: str, compress: bool = True) -> None:
    """Write log entries to pickle (optionally .pkl.gz)."""
    if compress and not path.endswith(".gz"):
        path = path + ".gz"
    payload = {"header": LOG_HEADER, "entries": [validate_log_entry(e) for e in entries]}
    if path.endswith(".gz"):
        with gzip.open(path, "wb") as f:
            pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(path, "wb") as f:
            pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    logger.info("Wrote %s log entries to %s", len(entries), path)


def read_log_pickle(path: str) -> list[dict[str, Any]]:
    """Read log entries from .pkl or .pkl.gz."""
    open_fn = gzip.open if path.endswith(".gz") else open
    mode = "rb"
    with open_fn(path, mode) as f:
        payload = pickle.load(f)
    return payload.get("entries", payload) if isinstance(payload, dict) else payload


def log_entries_to_sheets_values(entries: list[dict[str, Any]]) -> list[list[Any]]:
    """Return rows (header + data) for Sheets API values().update."""
    return log_entries_to_rows(entries)


def make_sample_log_entries() -> list[dict[str, Any]]:
    """Return simulated schedule log entries for inspection and testing."""
    return [
        {
            "date": "2025-03-14",
            "observatory": "LCO",
            "telescope": "Swope",
            "observer": "Observer1",
            "target_name": "SN2025A",
            "ra": "14:32:15.12",
            "dec": "-12:05:22.3",
            "filter": "r",
            "exptime_sec": 300,
            "priority": 0.0015,
            "notes": "First epoch",
        },
        {
            "date": "2025-03-14",
            "observatory": "LCO",
            "telescope": "Swope",
            "observer": "Observer1",
            "target_name": "SN2025A",
            "ra": "14:32:15.12",
            "dec": "-12:05:22.3",
            "filter": "g",
            "exptime_sec": 400,
            "priority": 0.0015,
            "notes": "",
        },
        {
            "date": "2025-03-14",
            "observatory": "Lick",
            "telescope": "Nickel",
            "observer": "Observer2",
            "target_name": "GW240315_field_001",
            "ra": "12:10:05.00",
            "dec": "+05:30:00.0",
            "filter": "i",
            "exptime_sec": 600,
            "priority": 0.002,
            "notes": "GW follow-up",
        },
        {
            "date": "2025-03-14",
            "observatory": "CTIO",
            "telescope": "Blanco",
            "observer": "",
            "target_name": "newfirm_field_042",
            "ra": "04:22:10.50",
            "dec": "-18:15:30.0",
            "filter": "J",
            "exptime_sec": 120,
            "priority": 1.0,
            "notes": "NEWFIRM survey",
        },
        {
            "date": "2025-03-15",
            "observatory": "LCO",
            "telescope": "Swope",
            "observer": "Observer1",
            "target_name": "SN2025B",
            "ra": "20:05:33.00",
            "dec": "-33:22:11.0",
            "filter": "r",
            "exptime_sec": 300,
            "priority": 0.001,
            "notes": "",
        },
    ]


def main():
    """CLI: generate sample schedule logs (CSV + compressed pickle) for inspection."""
    import argparse
    import os
    ap = argparse.ArgumentParser(
        description="Generate sample schedule log files (CSV and .pkl.gz) for inspection."
    )
    ap.add_argument(
        "-o", "--outdir",
        default=".",
        help="Output directory (default: current directory)",
    )
    ap.add_argument(
        "--no-pickle",
        action="store_true",
        help="Do not write the .pkl.gz file",
    )
    args = ap.parse_args()
    entries = make_sample_log_entries()
    csv_path = os.path.join(args.outdir, "schedule_log.csv")
    write_log_csv(entries, csv_path)
    if not args.no_pickle:
        pkl_path = os.path.join(args.outdir, "schedule_log.pkl.gz")
        write_log_pickle(entries, pkl_path, compress=True)
    print("Generated sample logs in", os.path.abspath(args.outdir))
    print("  ", csv_path)
    if not args.no_pickle:
        print("  ", pkl_path)
    print("CSV can be opened in Excel or imported into Google Sheets (File > Import).")
    print("Load pickle in Python: from scheduler.common.log_format import read_log_pickle")
