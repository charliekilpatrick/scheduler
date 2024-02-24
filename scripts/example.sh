#!/bin/bash
# example.sh
# example commands for running the scheduler on Swope, Nickel, Thacher, T80S

python CreateSchedule.py -f data/S230518h_SWOPE_4D_0.9_bayestar.fits.gz.txt --gw 37.19666 --target -14.5 -fc SWOPE_S230518h.fieldcenters --obstele LCO:Swope --outdir /Users/ckilpatrick/scripts/python/scheduler/output/S230518h
python CreateSchedule.py -f data/S230518h_NICKEL_4D_0.9_bayestar.fits.gz.txt --gw 37.19666 --target -14.5 -fc NICKEL_S230518h.fieldcenters --obstele Lick:Nickel --outdir /Users/ckilpatrick/scripts/python/scheduler/output/S230518h
python CreateSchedule.py -f data/S230518h_THACHER_4D_0.9_bayestar.fits.gz.txt --gw 37.19666 --target -14.5 -fc THACHER_S230518h.fieldcenters --obstele Thacher:Thacher --outdir /Users/ckilpatrick/scripts/python/scheduler/output/S230518h
python CreateSchedule.py -f data/S230518h_T80S_4D_0.9_bayestar.fits.gz.txt --gw 37.19666 --target -14.5 -fc T80S_S230518h.fieldcenters --obstele CTIO:T80S --outdir /Users/ckilpatrick/scripts/python/scheduler/output/S230518h
