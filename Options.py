def add_options():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
        help="CSV file with targets to schedule.  Can be a URL.")
    parser.add_argument("-d", "--date",
        help="YYYYMMDD formatted observation date [default is today].")
    parser.add_argument("-ot", "--obstele",
        help="Comma-delimited list of <Observatory>:<Telescope>.")
    parser.add_argument("-pp", "--plot",
        help="Preview the plot during command line execution.",
        action='store_true')
    parser.add_argument("-a", "--now",
        help="Start Now -- True or False")
    parser.add_argument("-b", "--start",
        help="Desired Start Time in the format of HHMM")
    parser.add_argument("-c", "--end",
        help="Desired End Time in the format of HHMM")
    parser.add_argument("--post",
        help="Post the schedule to the Google sheet for Swope or Nickel")
    parser.add_argument("--gw",
        help="The input targets are for a gravitational wave event.  Also"+\
        " append the preferred distance modulus for the event.")
    parser.add_argument("-o", "--output",
        help="Base name of the output schedule file (csv) and summary (png).")
    parser.add_argument("-fc", "--fieldcenter",
        help="Name of the output ordered fieldcenters file with targets.")
    parser.add_argument("--target", default=-17.0, type=float,
        help="Target mag for gravitational wave observations.")
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
    parser.add_argument("--halimit", default=None,
        help="Input hour angle limit for telescopes with limit.")
    parser.add_argument("--username", default='', type=str,
        help="Username for parsing input target list when it is a URL.")
    parser.add_argument("--password", default='', type=str,
        help="Password for parsing input target list when it is a URL.")
    
    args = parser.parse_args()

    return(args)

def parse_cat_params(args):
    return({'usecat': args.usecat, 'cat_radius': args.cat_radius,
        'cat_directory': args.cat_dir, 'cat_maxmag': args.cat_maxmag, 
        'cat_minmag': args.cat_minmag, 'cat_nsources': args.cat_nsources})
