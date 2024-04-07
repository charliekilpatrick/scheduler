# Internal dependencies
from common import Constants
from common import Options
from common import Utilities
from common.Target import TargetType, Target, handle_target_type
from common.observatory_defs import observatories

# Other dependencies
from dateutil.parser import parse
from astropy.coordinates import SkyCoord
from astropy.table import vstack
from astropy.time import Time
from astropy import units as unit
import os
import sys
import warnings
warnings.filterwarnings('ignore')

def get_target_data(args):
    # If a target list is provided via the file name then use it
    target_data = None
    if args.file is not None:
        for file in args.file:
            targ_data = Utilities.get_targets(file, gw=args.gw, 
                target_mag=args.target_mag, username=args.username, 
                password=args.password, newfirm=args.newfirm)
            if target_data is None:
                target_data = targ_data
            else:
                target_data = vstack([target_data, targ_data])
    # Otherwise download a target list from YSE PZ
    else:
        message = '\n\nDownloading target list for {tel}...\n\n'
        print(message.format(tel=args.tele_keys[0]))
        target_data = download_targets(args.tele_keys[0])

    # Check that we got some target data
    if target_data is None:
        error = 'ERROR: could not load or download any target data!'
        print(error)
        sys.exit(1)

    return(target_data)

def main():

    args = Options.handle_options()
    target_data = get_target_data(args)

    for i,obs in enumerate(args.observatories):

        targets = []
        for target in target_data:

            coord = Utilities.parse_coord(target['ra'], target['dec'])
            target_type, disc_date = handle_target_type(target)

            targets.append(
                Target(
                    name=target['name'],
                    coord=coord,
                    priority=target['priority'],
                    target_type=target_type,
                    observatory_lat=obs.ephemeris.lat,
                    sidereal_radian_array=obs.sidereal_radian_array,
                    ref_date=None,
                    apparent_mag=target['mag'],
                    halimit=args.halimit,
                    orig_priority=target['orig_priority']
                )
            )

        obs.telescopes[args.tele_keys[i]].set_targets(targets)

        print("# of %s targets: %s" % (args.tele_keys[i], len(targets)))
        print("First %s target: %s" % (args.tele_keys[i], targets[0].name))
        print("Last %s target: %s" % (args.tele_keys[i], targets[-1].name))

        obs.schedule_targets(args.tele_keys[i], args.plot, outdir=args.outdir,
            output_files=args.output, fieldcenters=args.fieldcenter,
            cat_params=args.cat_params, obs_date=args.date,
            start_time=args.start, end_time=args.end,
            first=args.first, second=args.second,
            minimize_slew=args.minimize_slew,
            one_off=args.one_off)

if __name__ == "__main__": main()

