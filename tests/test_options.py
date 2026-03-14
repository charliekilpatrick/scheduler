"""Tests for common.options."""

import pytest
from scheduler.common.options import parse_cat_params


class TestParseCatParams:
    def test_parse_cat_params(self):
        class Args:
            usecat = True
            cat_radius = 0.1
            cat_dir = "/tmp/cats"
            cat_maxmag = 19.0
            cat_minmag = 17.0
            cat_nsources = 10
        args = Args()
        d = parse_cat_params(args)
        assert d["usecat"] is True
        assert d["cat_radius"] == 0.1
        assert d["cat_directory"] == "/tmp/cats"
        assert d["cat_maxmag"] == 19.0
        assert d["cat_minmag"] == 17.0
        assert d["cat_nsources"] == 10
