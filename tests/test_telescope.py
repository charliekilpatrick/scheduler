"""Tests for Telescope base class."""

import pytest
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

from scheduler.observatories import Telescope
from scheduler.common.target import Target, TargetType


class ConcreteTelescope(Telescope):
    def __init__(self):
        self.targets = None
        self.name = "Test"
        self.filters = {}
        self.exp_funcs = {}

    def set_targets(self, targets):
        self.targets = targets

    def get_targets(self):
        return self.targets

    def compute_exposures(self):
        pass

    def write_schedule(self, observatory_name, obs_date, good_targets,
                      output_files=None, fieldcenters=None, pointing=None):
        pass


class TestTelescopeBase:
    def test_filter_row(self):
        t = ConcreteTelescope()
        row = t.filter_row("r", 300)
        assert row == [None, None, None, None, "r", 300]

    def test_round_to_num(self):
        t = ConcreteTelescope()
        assert t.round_to_num(20, 45) == 40
        assert t.round_to_num(20, 55) == 60

    def test_time_to_S_N(self):
        t = ConcreteTelescope()
        exp = t.time_to_S_N(30.0, 18.0, 25.0, min_exp=45.0)
        assert exp >= 45.0
        exp2 = t.time_to_S_N(30.0, 18.0, 25.0, min_exp=10.0)
        assert exp2 >= 10.0

    def test_limiting_magnitude(self):
        t = ConcreteTelescope()
        m = t.limiting_magnitude(25.0, 300.0, 3.0)
        assert 15.0 < m < 35.0

    def test_set_base_overhead(self):
        t = ConcreteTelescope()
        t.set_base_overhead(120.0)
        assert t.base_overhead == 120.0

    def test_set_fixed_exptime(self):
        t = ConcreteTelescope()
        t.set_fixed_exptime(300.0)
        assert t.fixed_exptime == 300.0
