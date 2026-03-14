"""Tests for common.constants."""

import pytest
import scheduler.common.constants as C


class TestConstants:
    def test_airmass_threshold(self):
        assert C.airmass_threshold == 1.75

    def test_round_to(self):
        assert C.round_to == 20

    def test_filter_names(self):
        assert C.u_band == "u"
        assert C.B_band == "B"
        assert C.r_prime == "r'"
        assert C.J_band == "J"
