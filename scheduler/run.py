"""
Command-line entry point for creating telescope schedules.

This module is used by the create-schedule console script.
"""

import sys
import warnings

warnings.filterwarnings("ignore")

from pathlib import Path

from scheduler.common import options, utilities
from scheduler.common.target import TargetType, Target, handle_target_type
from scheduler.observatories import observatories
from scheduler.common.logging import get_logger
from scheduler.plot_airmass import plot_airmass_curves
from astropy.table import vstack

logger = get_logger(__name__)


def get_target_data(args):
    """Load or download target data from args.file; returns Table or exits."""
    target_data = None
    if args.file is not None:
        for file in args.file:
            targ_data = utilities.get_targets(
                file,
                gw=args.gw,
                target_mag=args.target_mag,
                username=args.username,
                password=args.password,
                newfirm=args.newfirm,
            )
            if target_data is None:
                target_data = targ_data
            else:
                target_data = vstack([target_data, targ_data])
    else:
        logger.info("Downloading target list for %s...", args.tele_keys[0])
        logger.error(
            "Target list download (e.g. from YSE PZ) is not implemented. "
            "Please provide target file(s) via -f / --file."
        )
        raise SystemExit(
            "ERROR: Target list download (e.g. from YSE PZ) is not implemented. "
            "Please provide target file(s) via -f / --file."
        )

    if target_data is None:
        logger.error("Could not load or download any target data")
        sys.exit(1)

    return target_data


def main(args=None):
    """
    Run scheduling for the given or parsed options.

    If args is None, parses CLI via options.handle_options().
    Otherwise uses the provided Namespace (e.g. from options.build_args) so
    the scheduler can be driven programmatically for any observatory/telescope.

    Parameters
    ----------
    args : argparse.Namespace or None
        If None, options are parsed from sys.argv. Otherwise must have
        observatories, tele_keys, file, outdir, etc. (same shape as handle_options).
    """
    if args is None:
        args = options.handle_options()
    target_data = get_target_data(args)

    for i, obs in enumerate(args.observatories):
        obs.set_airmass_threshold(args.airmass_limit)
        obs.set_halimit_threshold(args.halimit)

        for key in obs.telescopes.keys():
            if args.fixed_overhead:
                obs.telescopes[key].set_base_overhead(args.fixed_overhead)
            if args.fixed_exptime:
                obs.telescopes[key].set_fixed_exptime(args.fixed_exptime)

        targets = []
        for target in target_data:
            coord = utilities.parse_coord(target["ra"], target["dec"])
            target_type, disc_date = handle_target_type(target)
            targets.append(
                Target(
                    name=target["name"],
                    coord=coord,
                    priority=target["priority"],
                    target_type=target_type,
                    observatory_lat=obs.ephemeris.lat,
                    sidereal_radian_array=obs.sidereal_radian_array,
                    ref_date=None,
                    apparent_mag=target["mag"],
                    orig_priority=target["orig_priority"],
                )
            )

        obs.telescopes[args.tele_keys[i]].set_targets(targets)

        logger.info("# of %s targets: %s", args.tele_keys[i], len(targets))
        logger.info("First %s target: %s", args.tele_keys[i], targets[0].name)
        logger.info("Last %s target: %s", args.tele_keys[i], targets[-1].name)

        scheduled = obs.schedule_targets(
            args.tele_keys[i],
            args.plot,
            outdir=args.outdir,
            output_files=args.output,
            fieldcenters=args.fieldcenter,
            cat_params=args.cat_params,
            obs_date=args.date,
            start_time=args.start,
            end_time=args.end,
            first=args.first,
            second=args.second,
            minimize_slew=args.minimize_slew,
            one_off=args.one_off,
        )
        if getattr(args, "plot_airmass", False) and scheduled:
            outpath = Path(args.outdir).resolve() / "airmass_curves.html"
            title = f"Airmass — {args.tele_keys[i]} {args.date}"
            fig = plot_airmass_curves(obs, scheduled, outpath=outpath, title=title)
            fig.show()


if __name__ == "__main__":
    main()
