"""
Programmatic schedule build for a given date and observatory/telescope.

Calls the main scheduler with options.build_args so any obstele can be
driven without hard-coded scripts.
"""
from scheduler.common import options
from scheduler.common._datetime_parse import parse_date_flexible
from scheduler.common.logging import get_logger
from scheduler import run
from scheduler.finders import make_finders

logger = get_logger(__name__)


def build_schedule_for_date(
    date,
    obstele,
    target_files,
    outdir=".",
    *,
    first=False,
    second=False,
    extra_target_files=None,
    newfirm=False,
    minimize_slew=False,
    one_off=False,
    output=None,
    compress=False,
    finders=False,
    finder_target_file=None,
    finder_directory=None,
    finder_clobber=False,
):
    """
    Build the scheduler run for one observatory/telescope on a given date.

    Parameters
    ----------
    date : str or datetime
        Observation date (YYYYMMDD or YYYY-MM-DD).
    obstele : str
        Observatory:telescope (e.g. "CTIO:Blanco", "Lick:Nickel").
    target_files : list of str or str
        Target file path(s) to pass to the scheduler.
    outdir : str
        Output directory. Default ".".
    first : bool
        First half of night only.
    second : bool
        Second half of night only.
    extra_target_files : list of str or None
        Additional target files to append to target_files.
    newfirm : bool
        NEWFIRM survey strategy.
    minimize_slew : bool
        Minimize slew time.
    one_off : bool
        One-off NEWFIRM serialization.
    output : str or None
        Base name for output schedule/summary files.
    compress : bool
        If True, create a tarball of outdir (implementation may be
        convention-based; e.g. newfirm_YYYYMMDD subdir).
    finders : bool
        If True, run finder chart generation after scheduling.
    finder_target_file : str or None
        Path to schedule/target file for finders (must have field, ra, dec, epoch).
        If None and finders is True, finders are skipped with a warning.
    finder_directory : str or None
        Output directory for finder PNGs. If None and finders is True, uses "finders".
    finder_clobber : bool
        Overwrite existing finder files.

    Returns
    -------
    None
    """
    if hasattr(date, "strftime"):
        dt = date
    else:
        dt = parse_date_flexible(str(date))
    date_str = dt.strftime("%Y%m%d")
    files = list(target_files) if isinstance(target_files, (list, tuple)) else [target_files]
    if extra_target_files:
        files.extend(
            list(extra_target_files)
            if isinstance(extra_target_files, (list, tuple))
            else [extra_target_files]
        )
    args = options.build_args(
        obstele,
        file=files,
        date=date_str,
        outdir=outdir,
        output=output,
        first=first,
        second=second,
        newfirm=newfirm,
        minimize_slew=minimize_slew,
        one_off=one_off,
    )
    run.main(args)
    if compress:
        import os
        import subprocess
        import tarfile
        tarname = f"schedule_{obstele.replace(':', '_')}_{date_str}.tar.gz"
        with tarfile.open(os.path.join(outdir, tarname), "w:gz") as tf:
            for name in os.listdir(outdir):
                path = os.path.join(outdir, name)
                if os.path.isdir(path) and not name.startswith("."):
                    tf.add(path, arcname=name)
        logger.info("Created %s", tarname)
    if finders:
        if not finder_target_file:
            logger.warning("finders=True but finder_target_file not set; skipping finders.")
        else:
            make_finders(
                finder_target_file,
                directory=finder_directory or "finders",
                clobber=finder_clobber,
            )


def main():
    """CLI: build schedule for a date and observatory/telescope (generic, not hard-coded)."""
    import argparse
    ap = argparse.ArgumentParser(description="Build schedule for a date and observatory:telescope.")
    ap.add_argument("date", help="Date (YYYYMMDD or YYYY-MM-DD)")
    ap.add_argument("-ot", "--obstele", required=True, help="Observatory:telescope (e.g. CTIO:Blanco, Lick:Nickel)")
    ap.add_argument("-f", "--file", nargs="+", required=True, dest="target_files", help="Target file(s)")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--first", action="store_true", help="First half of night only")
    ap.add_argument("--second", action="store_true", help="Second half of night only")
    ap.add_argument("--extra-file", nargs="*", dest="extra_target_files", help="Extra target file(s)")
    ap.add_argument("--newfirm", action="store_true", help="NEWFIRM survey strategy")
    ap.add_argument("--minimize-slew", action="store_true", help="Minimize slew time")
    ap.add_argument("--one-off", action="store_true", help="One-off NEWFIRM serialization")
    ap.add_argument("-o", "--output", default=None, help="Base name for output schedule/summary files")
    ap.add_argument("--compress", action="store_true", help="Create tarball of output directory")
    ap.add_argument("--finders", action="store_true", help="Generate finder charts after scheduling")
    ap.add_argument("--finder-target-file", default=None, help="Schedule/target file for finders (field,ra,dec,epoch)")
    ap.add_argument("--finder-directory", default=None, help="Output directory for finder PNGs")
    ap.add_argument("--finder-clobber", action="store_true", help="Overwrite existing finder files")
    args = ap.parse_args()
    build_schedule_for_date(
        args.date,
        args.obstele,
        args.target_files,
        outdir=args.outdir,
        first=args.first,
        second=args.second,
        extra_target_files=args.extra_target_files or None,
        newfirm=args.newfirm,
        minimize_slew=args.minimize_slew,
        one_off=args.one_off,
        output=args.output,
        compress=args.compress,
        finders=args.finders,
        finder_target_file=args.finder_target_file,
        finder_directory=args.finder_directory,
        finder_clobber=args.finder_clobber,
    )
