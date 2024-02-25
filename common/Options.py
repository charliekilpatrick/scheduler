from common.observatory_defs import observatories

def add_options():
    import argparse
    parser = argparse.ArgumentParser()
    # Options for input files
    parser.add_argument("-f", "--file", nargs="+", type=str,
        help="""astropy.io.ascii interpretable files with targets to schedule.  
        Can also be a URL to a astropy.io.ascii interpretable file.""")
    parser.add_argument("--username", default='', type=str,
        help="Username for parsing input target list when it is a URL.")
    parser.add_argument("--password", default='', type=str,
        help="Password for parsing input target list when it is a URL.")
    
    # Option for selecting output based on observatory/telescope
    parser.add_argument("-ot", "--obstele",
        help="Comma-delimited list of <Observatory>:<Telescope>.")
    
    # Options for constraining date and start/end time of schedule
    parser.add_argument("-d", "--date",
        help="YYYYMMDD formatted observation date [default is today].")
    parser.add_argument("-a", "--now",
        help="Start Now -- True or False")
    parser.add_argument("-b", "--start",
        help="Desired Start Time in the format of HHMM")
    parser.add_argument("-c", "--end",
        help="Desired End Time in the format of HHMM")
    parser.add_argument("--first", default=False, action='store_true',
        help="""Only schedule targets for first half (overrides --start,
        --end, and conflicts with --second).""")
    parser.add_argument("--second", default=False, action='store_true',
        help="""Only schedule targets for second half (overrides --start,
        --end, and conflicts with --first).""")

    # Option to indicate that targets should be considered as NEWFIRM survey
    # fields if a target type is not otherwise specified
    parser.add_argument("--newfirm", default=False, action='store_true',
        help="Schedule targets for NEWFIRM survey strategy.")

    # GW-specific options
    parser.add_argument("--gw", nargs='?', const=1.0, type=float,
        help="Flag that input targets are for a gravitational wave event.  Also"+\
        " append the preferred distance modulus for the event.")
    parser.add_argument("--target", default=-17.0, type=float,
        help="Target limiting mag for gravitational wave observations.")

    # Options for constraining target selection during scheduling
    parser.add_argument("--halimit", default=None,
        help="Input hour angle limit for telescopes with limit.")
    
    # Options for using catalogs to select fields based on number of available
    # calibrators above a threshold magnitude (see Nickel).
    parser.add_argument("--usecat", default=False, action='store_true',
        help="Use a PS1 catalog to decide whether to observe a field.")
    parser.add_argument("--cat-dir", default='', type=str,
        help="Directory to store target-specific catalogs for target selection.")
    parser.add_argument("--cat-radius", default=0.05, type=float,
        help="Check sources in catalog within this radius to see if should use.")
    parser.add_argument("--cat-minmag", default=18.0, type=float,
        help="Count number of sources in catalog up to this maximum magnitude.")
    parser.add_argument("--cat-maxmag", default=18.0, type=float,
        help="Count number of sources in catalog up to this maximum magnitude.")
    parser.add_argument("--cat-nsources", default=9, type=int,
        help="Require a minimum of this number of sources to schedule.")
    
    # Output files and plots
    parser.add_argument("--outdir", type=str, default='.',
        help="Base directory for output from the scheduler.")
    parser.add_argument("-o", "--output",
        help="Base name of the output schedule file (txt) and/or summary (png).")
    parser.add_argument("-fc", "--fieldcenter",
        help="Name of the output ordered fieldcenters file with targets.")
    parser.add_argument("--post",
        help="Post the schedule to the Google sheet for Swope or Nickel")
    parser.add_argument("-pp", "--plot",
        help="Preview the plot during command line execution.",
        action='store_true')
    parser.add_argument("--minimize-slew", default=False, action="store_true",
        help="Minimize slew time between targets while creating schedule.")
    parser.add_argument("--newfirm-dir", default='newfirm', type=str,
        help="Output directory for NEWFIRM observing scripts.")
    
    args = parser.parse_args()

    return(args)

def parse_cat_params(args):
    return({'usecat': args.usecat, 'cat_radius': args.cat_radius,
        'cat_directory': args.cat_dir, 'cat_maxmag': args.cat_maxmag, 
        'cat_minmag': args.cat_minmag, 'cat_nsources': args.cat_nsources})

def handle_options():

    args = add_options()

    # Check for conflicts between options
    if args.obstele is None:
        raise Exception(f'ERROR: --obstele cannot be None')
    if args.first and args.second:
        raise Exception(f'ERROR: cannot pass both --first and --second.')

    if args.first or args.second: args.start = None ; args.end = None
    
    if args.target:
        args.target_mag = float(args.target)
    else:
        args.target_mag = -17.0

    args.cat_params = parse_cat_params(args)

    observatory_telescopes = args.obstele.split(",")
    obs_keys = [o.split(":")[0] for o in observatory_telescopes]
    tele_keys = [t.split(":")[1] for t in observatory_telescopes]

    args.obs_keys = obs_keys
    args.tele_keys = tele_keys
    args.observatories = [observatories[obs_key] for obs_key in obs_keys]

    return(args)
