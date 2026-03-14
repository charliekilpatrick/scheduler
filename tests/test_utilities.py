"""Tests for common.utilities."""

import pytest
from astropy.coordinates import SkyCoord
from astropy import units as u

from scheduler.common.utilities import (
    is_number,
    parse_coord,
    blank_target_table,
    UTC_Offset,
)


class TestIsNumber:
    def test_integer(self):
        assert is_number(42) is True
        assert is_number("42") is True

    def test_float(self):
        assert is_number(3.14) is True
        assert is_number("-1.5") is True

    def test_non_number(self):
        assert is_number("abc") is False
        assert is_number("12:30:00") is False


class TestParseCoord:
    def test_degree_input(self):
        c = parse_coord(180.0, -45.0)
        assert c is not None
        assert abs(c.ra.degree - 180.0) < 1e-6
        assert abs(c.dec.degree - (-45.0)) < 1e-6

    def test_sexagesimal_input(self):
        c = parse_coord("12:00:00", "45:00:00")
        assert c is not None
        assert abs(c.ra.hour - 12.0) < 1e-6
        assert abs(c.dec.degree - 45.0) < 1e-6

    def test_invalid_returns_none(self):
        assert parse_coord("x", "y") is None


class TestBlankTargetTable:
    def test_columns(self):
        tab = blank_target_table()
        assert "name" in tab.colnames
        assert "ra" in tab.colnames
        assert "dec" in tab.colnames
        assert "priority" in tab.colnames

    def test_empty(self):
        tab = blank_target_table()
        assert len(tab) == 0


class TestUTCOffset:
    def test_offset(self):
        tz = UTC_Offset(-5, "EST")
        from datetime import datetime
        dt = datetime(2020, 1, 1, 12, 0, 0)
        off = tz.utcoffset(dt)
        assert off.total_seconds() == -5 * 3600
