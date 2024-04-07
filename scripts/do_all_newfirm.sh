#!/bin/bash
# do_all_newfirm.sh

#datedir=`date +%Y%m%d`
#date=`date +%Y-%m-%d`
datedir="20240227"
date="2024-02-27"
rsync -avzh ckilpatrick@burbidge.northwestern.edu:/data/NEWFIRM/targets/newfirm_good_templates_$datedir.ecsv $NEWFIRM_DIR/targlists/newfirm_good_templates_$datedir.ecsv
/Users/ckilpatrick/Dropbox/scripts/python/gw/scheduler/scripts/newfirm_schedule.py $date --second --compress --finders --supernova-targets $NEWFIRM_DIR/targlists/supernovae.txt
