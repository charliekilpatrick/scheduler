# Internal dependencies
from common import Constants
from common import Options
from common import Utilities
from common.Target import TargetType, Target
from common.observatory_defs import observatories

# Other dependencies
from dateutil.parser import parse
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as unit
import os
import sys
import warnings
warnings.filterwarnings('ignore')

def main():

    args = Options.add_options()

    if args.obstele is None:
        raise Exception(f'ERROR: --obstele cannot be None')

    if args.first and args.second:
        raise Exception(f'ERROR: cannot pass both --first and --second.')

    file_name = args.file
    obs_date = args.date
    observatory_telescopes = args.obstele.split(",")
    preview_plot = args.plot
    fieldcenters = args.fieldcenter
    output_files = args.output
    startNow = args.now in ['True']
    startTime = args.start
    endTime = args.end
    first = args.first
    second = args.second
    outdir = args.outdir

    if first or second: start = None ; end = None
    
    if args.target:
        target_mag = float(args.target)
    else:
        target_mag = -17.0

    obs_keys = [o.split(":")[0] for o in observatory_telescopes]
    tele_keys = [t.split(":")[1] for t in observatory_telescopes]

    cat_params = Options.parse_cat_params(args)

    # If a target list is provided via the file name then use it
    if file_name is not None:
        target_data = Utilities.get_targets(file_name, gw=args.gw, 
            target_mag=target_mag, username=args.username, 
            password=args.password, newfirm=args.newfirm)
    # Otherwise download a target list from YSE PZ
    else:
        message = '\n\nDownloading target list for {tel}...\n\n'
        print(message.format(tel=tele_keys[0]))
        target_data = download_targets(tele_keys[0])

    # Check that we got some target data
    if target_data is None:
        error = 'ERROR: could not load or download any target data!'
        print(error)
        sys.exit(1)

    for i in range(len(observatory_telescopes)):

        targets = []
        obs = observatories[obs_keys[i]]

        for target in target_data:

            target_type = None
            disc_date = None

            coord = SkyCoord(target['ra'], target['dec'],
                unit='deg')

            if target['type'] == 'STD':
                target_type = TargetType.Standard
                disc_date = None
            elif target['type'] == 'NEWFIRM':
                target_type = TargetType.NEWFIRM
                disc_date = None
            elif target['type'] == 'TMP':
                target_type = TargetType.Template
            elif target['type'] == 'SN':
                target_type = TargetType.Supernova
            elif target['type'] == 'GW':
                target_type = TargetType.GW
            else:
                raise ValueError('Unrecognized target type!')

            if 'orig_priority' in target.colnames:
                orig_priority=target['orig_priority']
            else:
                orig_priority=None

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
                    orig_priority=orig_priority
                )
            )

        obs.telescopes[tele_keys[i]].set_targets(targets)

        print("# of %s targets: %s" % (tele_keys[i], len(targets)))
        print("First %s target: %s" % (tele_keys[i], targets[0].name))
        print("Last %s target: %s" % (tele_keys[i], targets[-1].name))

        obs.schedule_targets(tele_keys[i], preview_plot, outdir=outdir,
            output_files=output_files, fieldcenters=fieldcenters,
            cat_params=cat_params, obs_date=obs_date,
            start_time=args.start, end_time=args.end,
            first=first, second=second,
            minimize_slew=args.minimize_slew)

if __name__ == "__main__": main()

