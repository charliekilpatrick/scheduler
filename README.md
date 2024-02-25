# Scheduler

Telescope scheduling code based on Dave Coulter's scheduler (https://github.com/davecoulter/Scheduler), with help from Tiara Hung, Jay Sarva, Karelle Sielez.

### Installation

Install dependencies via `requirements.txt`.

### Usage

Basic examples are given under `scripts/example.sh`.  Given a `astropy.io.ascii` interpretable target file `targets.txt` with field name, RA, and Dec, the scheduler can create a target list for a specific telescope (e.g., `LCO:Swope`) for the date August 17, 2017 with:

```
python CreateSchedule.py -f targets.txt --obstele LCO:Swope --date 2017-08-17
```

Target-specific options for adjusting exposure times, target limiting magnitude, and observing bands can be explored via the options (`python CreateSchedule.py -h`) or in the specific telescope modules in `common`.

### Contact

For bugs, issues, or questions about this code, contact Charlie Kilpatrick at ckilpatrick@northwestern.edu.
