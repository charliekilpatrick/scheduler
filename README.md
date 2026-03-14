# Scheduler

Telescope scheduling code for observatory target scheduling based on airmass and priority. Derived from [Dave Coulter's Scheduler](https://github.com/davecoulter/Scheduler), with contributions from Tiara Hung, Jay Sarva, and Karelle Sielez.

## Features

- **Observatories and telescopes**: LCO (Swope), Lick (Nickel), Thacher (Thacher), CTIO (T80S, Blanco/NEWFIRM).
- **Input**: Target lists from files or URLs in `astropy.io.ascii`-readable format (name, RA, Dec; optional priority, type, mag).
- **Output**: Schedule files (CSV, `.phot`, field centers) and telescope-specific formats (RTML for Thacher, NEWFIRM observing scripts for Blanco).
- **Options**: Minimize slew, first/second half of night, custom start/end times (HHMM), catalog-based field filtering (PS1), optional Google Sheets logging for Nickel.

## Requirements

- Python 3.12 or later
- Dependencies are declared in `pyproject.toml` (astropy, astroquery, plotly, kaleido, numpy, photutils, pyephem, etc.)

## Installation

### pip (recommended)

From the repository root:

```bash
pip install -e .
```

This installs the `scheduler` package in editable mode and registers the `create-schedule` command. Dependencies are installed automatically.

**Optional extras:**

- Google Sheets (for Nickel schedule logging):  
  `pip install -e ".[google]"`
- Development (pytest, pytest-cov):  
  `pip install -e ".[dev]"`
- Both:  
  `pip install -e ".[google,dev]"`

### Conda

1. Create and activate the environment:

   ```bash
   conda env create -f environment.yml
   conda activate scheduler
   ```

2. Install the package in development mode from the repo root:

   ```bash
   pip install -e .
   ```

The `environment.yml` includes Python and core dependencies (astropy, numpy, plotly, kaleido, scipy, etc.); astroquery, photutils, and pyephem are installed via pip in that environment.

## Usage

### Command line

**With the installed script:**

```bash
create-schedule -f targets.txt -ot LCO:Swope -d 2017-08-17
```

**From the repo root without installing** (add repo root to `PYTHONPATH` if needed):

```bash
python -m scheduler.run -f targets.txt --obstele LCO:Swope --date 2017-08-17
```

**Main options:**

| Option | Description |
|--------|-------------|
| `-f` / `--file` | One or more target list files (or URLs). Format: `astropy.io.ascii`-readable (name, RA, Dec; optional priority, type, mag). |
| `-ot` / `--obstele` | Comma-separated `Observatory:Telescope` (e.g. `LCO:Swope`, `Lick:Nickel`, `Thacher:Thacher`, `CTIO:T80S`, `CTIO:Blanco`). **Required.** |
| `-d` / `--date` | Observation date: `YYYYMMDD` or `YYYY-MM-DD`. Default: today. |
| `--outdir` | Output directory (default: current directory). |
| `-o` / `--output` | Base name for schedule output files. |
| `-fc` / `--fieldcenter` | Filename for the ordered field centers output. |
| `--first` / `--second` | Schedule only first or second half of night. |
| `-b` / `--start`, `-c` / `--end` | Custom start/end time as `HHMM`. |
| `--airmass-limit` | Maximum airmass (default 3.0). |
| `--minimize-slew` | Run slew minimization when building the schedule. |
| `--gw` [value] | Treat targets as gravitational-wave event; optional distance modulus. |
| `--target` | Target limiting magnitude for GW (default -17.0). |
| `--newfirm` | Treat unspecified targets as NEWFIRM survey fields. |
| `--one-off` | For NEWFIRM, output supernova observations individually. |
| `--fixed-exptime`, `--fixed-overhead` | Override exposure time (s) or overhead (s). |
| `--usecat`, `--cat-dir`, `--cat-*` | Use PS1 catalog for field filtering; see `-h`. |

Full list:

```bash
create-schedule -h
```

### Example commands (bash)

**GW scheduling** (edit paths and `--obstele`/`--outdir` for your setup):

```bash
create-schedule -f data/S230518h_SWOPE_4D_0.9_bayestar.fits.gz.txt --gw 37.19666 --target -14.5 \
  -fc SWOPE_S230518h.fieldcenters --obstele LCO:Swope --outdir output/S230518h
create-schedule -f data/S230518h_NICKEL_4D_0.9_bayestar.fits.gz.txt --gw 37.19666 --target -14.5 \
  -fc NICKEL_S230518h.fieldcenters --obstele Lick:Nickel --outdir output/S230518h
create-schedule -f data/S230518h_THACHER_4D_0.9_bayestar.fits.gz.txt --gw 37.19666 --target -14.5 \
  -fc THACHER_S230518h.fieldcenters --obstele Thacher:Thacher --outdir output/S230518h
create-schedule -f data/S230518h_T80S_4D_0.9_bayestar.fits.gz.txt --gw 37.19666 --target -14.5 \
  -fc T80S_S230518h.fieldcenters --obstele CTIO:T80S --outdir output/S230518h
```

**NEWFIRM (CTIO:Blanco)** — set `NEWFIRM_DIR` and use `scheduler-build-date`:

```bash
# Single date (e.g. today)
scheduler-build-date $(date +%Y-%m-%d) -ot CTIO:Blanco \
  -f $NEWFIRM_DIR/targlists/newfirm_good_templates_$(date +%Y%m%d).ecsv \
  --outdir $NEWFIRM_DIR/obsplans --newfirm --minimize-slew --one-off --second

# With finders and optional extra targets
scheduler-build-date 2024-02-23 -ot CTIO:Blanco \
  -f $NEWFIRM_DIR/targlists/newfirm_good_templates_20240223.ecsv \
  --extra-file $NEWFIRM_DIR/targlists/supernovae.txt \
  --outdir $NEWFIRM_DIR/obsplans --newfirm --minimize-slew --one-off --second \
  --finders --finder-target-file $NEWFIRM_DIR/obsplans/newfirm_20240223/newfirm_schedule_20240223.txt \
  --finder-directory $NEWFIRM_DIR/finders/newfirm_20240223
```

**Fetch NEWFIRM target list and build** (run from repo or path where scheduler is installed):

```bash
datedir=$(date +%Y%m%d)
date=$(date +%Y-%m-%d)
rsync -avzh user@host:/data/NEWFIRM/targets/newfirm_good_templates_$datedir.ecsv \
  $NEWFIRM_DIR/targlists/newfirm_good_templates_$datedir.ecsv
scheduler-build-date $date -ot CTIO:Blanco \
  -f $NEWFIRM_DIR/targlists/newfirm_good_templates_$datedir.ecsv \
  --outdir $NEWFIRM_DIR/obsplans --newfirm --minimize-slew --one-off --second
```

### Additional CLI commands (after `pip install -e .`)

| Command | Purpose |
|--------|---------|
| **`scheduler-convert-gw`** | Convert a GW tiles file to scheduler input CSV. `scheduler-convert-gw <input> <output> [--invert-priority] [--add-string STR]` |
| **`scheduler-filter-tiles`** | Filter a tiles file to scheduled tiles using a scheduler CSV. `-s scheduler.csv -t tiles_file -c output_file` |
| **`scheduler-probs`** | Compare schedule CSV to tiles file (FieldName, Priority); print coincidence and priority sum. `-s schedule.csv -f tiles_file` |
| **`scheduler-finders`** | Generate 2MASS-J finder charts for a target list (field, ra, dec, epoch). `scheduler-finders <target_list> [-d dir] [--clobber]` |
| **`scheduler-check-schedule`** | Print 1 if a schedule should be generated for obstele on date, 0 otherwise. `scheduler-check-schedule <date> [--obstele Lick:Nickel]` |
| **`scheduler-build-date`** | Build schedule for a date and observatory:telescope. `scheduler-build-date <date> -ot <obstele> -f <file> [file ...] [--outdir .] [--first \| --second] [--newfirm] [--finders] ...` |
| **`scheduler-sample-logs`** | Generate sample schedule log files (CSV + .pkl.gz) for inspection. `scheduler-sample-logs [-o outdir] [--no-pickle]` |

### Python API

**Programmatic build and helpers (generic for any observatory:telescope):**

```python
from scheduler.build import build_schedule_for_date
from scheduler.schedule_check import need_to_generate_schedule_for_obstele
from scheduler.gw_convert import convert_gw_to_scheduler, filter_tiles_by_schedule, schedule_coincidence_stats
from scheduler.finders import make_finders
from scheduler.io import read_scheduled_names_from_csv, SCHEDULER_CSV_HEADER

# Build schedule for a date and obstele (no hard-coded telescope)
build_schedule_for_date("20250115", "CTIO:Blanco", ["targets.ecsv"], outdir=".", newfirm=True, minimize_slew=True)

# Check if schedule should be generated (Lick:Nickel uses Sheets; others return True)
need_to_generate_schedule_for_obstele("Lick:Nickel", "20250115")

# GW conversion and comparison
convert_gw_to_scheduler("tiles.txt", "scheduler_input.csv", invert_priority=False)
filter_tiles_by_schedule("schedule.csv", "tiles.txt", "filtered_tiles.txt")
stats = schedule_coincidence_stats("schedule.csv", "tiles.csv")  # tiles need FieldName, Priority

# Finders (target list: field, ra, dec, epoch)
make_finders("targets.txt", directory="finders", clobber=False)

# Sample schedule logs (CSV / pickle) and log format
from scheduler.common.log_format import make_sample_log_entries, write_log_csv, write_log_pickle, read_log_pickle
entries = make_sample_log_entries()
write_log_csv(entries, "schedule_log.csv")
write_log_pickle(entries, "schedule_log.pkl.gz")
```

**Full scheduler run (parse CLI or build args):**

```python
from scheduler.common import constants, options, utilities
from scheduler.common.target import Target, TargetType, handle_target_type
from scheduler.observatories import observatories
from scheduler import run

# Option A: parse CLI
args = options.handle_options()

# Option B: build args programmatically for any obstele
args = options.build_args("LCO:Swope", file=["targets.txt"], date="20250115", outdir=".")

run.main(args)  # runs scheduling for args.observatories / args.tele_keys
```

**Loading targets and calling schedule_targets (lower level):**

```python
args = options.handle_options()

# Load targets from file/URL
target_data = utilities.get_targets(
    args.file[0],
    gw=args.gw,
    target_mag=args.target_mag,
    username=args.username,
    password=args.password,
    newfirm=args.newfirm,
)

# Build Target list and run scheduling for an observatory/telescope
obs = observatories["LCO"]
# ... build targets from target_data, then:
obs.telescopes["Swope"].set_targets(targets)
obs.schedule_targets(
    "Swope",
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
```

## Logging

The package uses Python’s `logging` module. Configuration is centralized in `scheduler/common/logging.py`, which provides `get_logger(name)`. All former `print` output is emitted via the `scheduler` logger hierarchy (e.g. `scheduler.observatories.base`).

- **Log level**: Set the environment variable `SCHEDULER_LOG_LEVEL` to `DEBUG`, `INFO`, `WARNING`, `ERROR`, or `CRITICAL`. Default is `INFO`.
- **Format**: Log messages are sent to stderr with a timestamp, level, logger name, and message.

Example:

```bash
SCHEDULER_LOG_LEVEL=DEBUG create-schedule -f targets.txt -ot LCO:Swope -d 2024-01-15
```

## Observatories package

Observatory and telescope code lives under `scheduler/observatories/`:

- **`sites.py`** — Site data (lat, lon, elevation, horizon, UTC offset, timezone name) for LCO, Lick, Thacher, CTIO, and other public observatories (McDonald, Siding Spring, SAAO, Kitt Peak, Mauna Kea, Paranal, La Silla, Palomar, Roque de los Muchachos, Teide, etc.). Use `get_site(site_key)` or `list_site_keys()`.
- **`base.py`** — `Observatory` (night boundaries, scheduling) and abstract `Telescope` base (exposures, CSV/phot/field centers).
- **`telescopes/`** — Concrete telescope implementations: Swope, Nickel, Thacher, T80S, Blanco.
- **`registry.py`** — Builds the `observatories` dict (LCO, Lick, Thacher, CTIO) from site data and `TELESCOPE_REGISTRY`; provides `get_observatory(site_key)` and `list_observatories()`.

**Adding a new observatory (site only):** Add an entry to `scheduler/observatories/sites.py` in `SITE_DATA` with keys: `name`, `lon`, `lat`, `elevation`, `horizon`, `utc_offset`, `utc_offset_name`. The site will be available via `get_observatory(site_key)` (with no telescopes until you add a driver).

**Adding a new telescope driver:** Implement a class in `scheduler/observatories/telescopes/` that subclasses `Telescope` and implements `set_targets`, `get_targets`, `compute_exposures`, and `write_schedule`. Register it in `registry.TELESCOPE_REGISTRY` under the appropriate site key (e.g. `"LCO": {"Swope": Swope}`). Re-export the class from `telescopes/__init__.py`.

## Repository layout

```
scheduler/                    # Main Python package
├── __init__.py               # Re-exports from scheduler.common and scheduler.observatories
├── run.py                    # CLI entry point (create-schedule)
├── observatories/            # Observatory and telescope definitions
│   ├── __init__.py           # Observatory, Telescope, observatories, get_observatory, etc.
│   ├── base.py               # Observatory and Telescope base classes
│   ├── sites.py               # Site data (lat, lon, elevation, UTC, timezone)
│   ├── registry.py            # observatories dict and TELESCOPE_REGISTRY
│   └── telescopes/            # Swope, Nickel, Thacher, T80S, Blanco
└── common/                   # Core subpackage
    ├── __init__.py
    ├── constants.py          # Airmass threshold, time grid, filter names
    ├── options.py            # CLI options and build_args
    ├── utilities.py         # Target I/O, parse_coord, PS1 catalog, dicttoxml
    ├── target.py            # TargetType, Target, handle_target_type
    ├── logging.py           # get_logger()
    ├── log_format.py        # Schedule log CSV/pickle/Sheets, make_sample_log_entries
    └── logs.py              # Google Sheets for Nickel (optional)

pyproject.toml                # Package metadata, dependencies, [project.scripts], pytest config
environment.yml               # Conda environment definition

tests/                        # Pytest suite for scheduler.common
├── __init__.py
├── test_constants.py
├── test_options.py
├── test_target.py
├── test_telescope.py
├── test_utilities.py
└── test_observatory_defs.py

.github/workflows/ci.yml      # CI: build, install, test on Python 3.12 and 3.13
LICENSE
.gitignore
```

## Testing

From the repository root (with the package installed, e.g. `pip install -e .`):

```bash
pytest tests/ -v
```

With coverage (install dev deps: `pip install -e ".[dev]"`):

```bash
pytest tests/ -v --cov=scheduler --cov-report=term-missing
```

**CI**: On push or pull request to `main` or `master`, `.github/workflows/ci.yml` runs on Ubuntu with Python 3.12 and 3.13: it installs the package with `pip install -e ".[dev]"` and runs `pytest tests/ -v --tb=short`.

## License

See `LICENSE` (MIT).

## Contact

For bugs, issues, or questions: **Charlie Kilpatrick**, ckilpatrick@northwestern.edu.
