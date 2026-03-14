"""CLI option definitions and parsing for create-schedule."""

from .logging import get_logger

logger = get_logger(__name__)


def add_options():
    """Build parser and parse sys.argv; returns argparse.Namespace."""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", nargs="+", type=str,
        help="Target list file(s) or URL (astropy.io.ascii-readable).")
    parser.add_argument("--username", default='', type=str, help="Username for URL auth.")
    parser.add_argument("--password", default='', type=str, help="Password for URL auth.")
    parser.add_argument("-ot", "--obstele",
        help="Comma-delimited Observatory:Telescope (e.g. LCO:Swope, Lick:Nickel).")
    parser.add_argument("-d", "--date", help="Observation date YYYYMMDD or YYYY-MM-DD.")
    parser.add_argument("-a", "--now", help="Start now True/False.")
    parser.add_argument("-b", "--start", help="Start time HHMM.")
    parser.add_argument("-c", "--end", help="End time HHMM.")
    parser.add_argument("--first", default=False, action='store_true',
        help="First half of night only.")
    parser.add_argument("--second", default=False, action='store_true',
        help="Second half of night only.")
    parser.add_argument("--newfirm", default=False, action='store_true',
        help="NEWFIRM survey strategy.")
    parser.add_argument("--gw", nargs='?', const=1.0, type=float,
        help="GW event flag; optional distance modulus.")
    parser.add_argument("--target", default=-17.0, type=float,
        help="Target limiting mag for GW.")
    parser.add_argument("--fixed-exptime", default=None, type=float, help="Fixed exposure (s).")
    parser.add_argument("--fixed-overhead", default=None, type=float, help="Fixed overhead (s).")
    parser.add_argument("--halimit", default=None, type=float, help="Hour angle limit.")
    parser.add_argument("--airmass-limit", default=3.0, type=float, help="Airmass limit.")
    parser.add_argument("--usecat", default=False, action='store_true', help="Use PS1 catalog.")
    parser.add_argument("--cat-dir", default='', type=str, help="Catalog directory.")
    parser.add_argument("--cat-radius", default=0.05, type=float, help="Catalog radius (deg).")
    parser.add_argument("--cat-minmag", default=18.0, type=float, help="Catalog min mag.")
    parser.add_argument("--cat-maxmag", default=18.0, type=float, help="Catalog max mag.")
    parser.add_argument("--cat-nsources", default=9, type=int, help="Min catalog sources.")
    parser.add_argument("--outdir", type=str, default='.', help="Output directory.")
    parser.add_argument("-o", "--output", help="Output schedule base name.")
    parser.add_argument("--post", help="Post to Google sheet.")
    parser.add_argument("-pp", "--plot", action='store_true', help="Preview plot.")
    parser.add_argument("-fc", "--fieldcenter", help="Field centers output file.")
    parser.add_argument("--minimize-slew", default=False, action="store_true", help="Minimize slew.")
    parser.add_argument("--plot-airmass", default=False, action="store_true",
        help="Plot airmass curves with scheduled segments highlighted (Plotly).")
    parser.add_argument("--newfirm-dir", default='newfirm', type=str, help="NEWFIRM output dir.")
    parser.add_argument("--one-off", default=False, action="store_true", help="NEWFIRM one-off.")
    return parser.parse_args()


def build_args(obstele, *, file=None, date=None, outdir=".", output=None,
               first=False, second=False, newfirm=False, gw=None, target_mag=-17.0,
               username="", password="", fixed_overhead=None, fixed_exptime=None,
               halimit=None, airmass_limit=3.0, usecat=False, cat_dir="",
               cat_radius=0.05, cat_minmag=18.0, cat_maxmag=18.0, cat_nsources=9,
               plot=False, fieldcenter=None, minimize_slew=False, plot_airmass=False,
               newfirm_dir="newfirm", one_off=False, start=None, end=None, post=None):
    """Build an args-like Namespace for programmatic runs (e.g. run.main(args))."""
    import argparse
    if first and second:
        raise ValueError("Cannot pass both first and second.")
    obs_keys = [obstele.split(":")[0]]
    tele_keys = [obstele.split(":")[1]]
    args = argparse.Namespace(
        obstele=obstele, file=file or [], date=date, outdir=outdir, output=output,
        first=first, second=second, newfirm=newfirm, gw=gw, target=target_mag,
        target_mag=float(target_mag), username=username, password=password,
        fixed_overhead=fixed_overhead, fixed_exptime=fixed_exptime, halimit=halimit,
        airmass_limit=airmass_limit, usecat=usecat, cat_dir=cat_dir, cat_radius=cat_radius,
        cat_minmag=cat_minmag, cat_maxmag=cat_maxmag, cat_nsources=cat_nsources,
        plot=plot, fieldcenter=fieldcenter, minimize_slew=minimize_slew,
        plot_airmass=plot_airmass, newfirm_dir=newfirm_dir, one_off=one_off,
        start=start if not (first or second) else None,
        end=end if not (first or second) else None, post=post, now=None,
    )
    args.cat_params = {
        "usecat": usecat, "cat_radius": cat_radius, "cat_directory": cat_dir,
        "cat_maxmag": cat_maxmag, "cat_minmag": cat_minmag, "cat_nsources": cat_nsources,
    }
    from scheduler.observatories import observatories
    args.obs_keys = obs_keys
    args.tele_keys = tele_keys
    args.observatories = [observatories[k] for k in obs_keys]
    return args


def parse_cat_params(args):
    """Return dict of catalog params from parsed args."""
    return {
        "usecat": args.usecat, "cat_radius": args.cat_radius,
        "cat_directory": args.cat_dir, "cat_maxmag": args.cat_maxmag,
        "cat_minmag": args.cat_minmag, "cat_nsources": args.cat_nsources,
    }


def handle_options():
    """Parse and validate CLI; resolve observatories/tele_keys from --obstele."""
    args = add_options()
    if args.obstele is None:
        logger.error("--obstele cannot be None")
        raise ValueError("--obstele cannot be None")
    if args.first and args.second:
        logger.error("Cannot pass both --first and --second")
        raise ValueError("Cannot pass both --first and --second.")
    if args.first or args.second:
        args.start = None
        args.end = None
    args.target_mag = float(args.target) if args.target else -17.0
    args.cat_params = parse_cat_params(args)
    from scheduler.observatories import observatories
    observatory_telescopes = args.obstele.split(",")
    args.obs_keys = [o.split(":")[0] for o in observatory_telescopes]
    args.tele_keys = [t.split(":")[1] for t in observatory_telescopes]
    args.observatories = [observatories[k] for k in args.obs_keys]
    return args
