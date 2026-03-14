"""Tests for common.target."""

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy import units as u

from scheduler.common.target import TargetType, handle_target_type, Target


class TestTargetType:
    def test_values(self):
        assert TargetType.Supernova.value == 1
        assert TargetType.GW.value == 4
        assert TargetType.NEWFIRM.value == 6


class TestHandleTargetType:
    def test_sn(self):
        tt, disc = handle_target_type({"type": "SN"})
        assert tt is TargetType.Supernova
        assert disc is None

    def test_gw(self):
        tt, disc = handle_target_type({"type": "GW"})
        assert tt is TargetType.GW
        assert disc is None

    def test_std(self):
        tt, disc = handle_target_type({"type": "STD"})
        assert tt is TargetType.Standard
        assert disc is None

    def test_newfirm(self):
        tt, disc = handle_target_type({"type": "NEWFIRM"})
        assert tt is TargetType.NEWFIRM
        assert disc is None

    def test_tmp(self):
        tt, disc = handle_target_type({"type": "TMP"})
        assert tt is TargetType.Template
        assert disc is None

    def test_invalid_raises(self):
        with pytest.raises(ValueError, match="Unrecognized"):
            handle_target_type({"type": "INVALID"})


class TestTarget:
    @pytest.fixture
    def coord(self):
        return SkyCoord(ra=180.0 * u.deg, dec=-30.0 * u.deg)

    @pytest.fixture
    def sidereal_array(self):
        return list(np.linspace(0, 2 * np.pi, 100))

    def test_init(self, coord, sidereal_array):
        t = Target(
            name="Test",
            coord=coord,
            priority=1.0,
            target_type=TargetType.GW,
            observatory_lat=0.5,
            sidereal_radian_array=sidereal_array,
            apparent_mag=18.0,
            orig_priority=0.5,
        )
        assert t.name == "Test"
        assert t.priority == 1.0
        assert t.type is TargetType.GW
        assert t.apparent_mag == 18.0
        assert t.raw_airmass_array is not None
        assert len(t.raw_airmass_array) == len(sidereal_array)

    def test_compute_airmass(self, coord, sidereal_array):
        t = Target(
            name="T",
            coord=coord,
            priority=1.0,
            target_type=TargetType.GW,
            observatory_lat=0.5,
            sidereal_radian_array=sidereal_array,
        )
        am = t.compute_airmass(0.5, np.array(sidereal_array), halimit=6.0)
        assert len(am) == len(sidereal_array)
        assert np.all((am >= 1.0) | (am > 10.0) | (am == 9999))
