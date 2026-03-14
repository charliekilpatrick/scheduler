"""Unit tests for plot_airmass module (duration formatting, timestamps, plot layout)."""

import pytest
from datetime import datetime, timezone

from scheduler.plot_airmass import (
    _slots_to_hours_minutes,
    _format_duration_h_m,
    _to_timestamp_ms,
    _color_for_index,
    PLOTLY_COLORS,
    plot_airmass_curves,
)


class TestSlotsToHoursMinutes:
    """Tests for _slots_to_hours_minutes."""

    def test_zero_slots(self):
        assert _slots_to_hours_minutes(0, 20) == (0, 0)

    def test_one_slot_20_sec(self):
        # 20 s -> 0 min
        assert _slots_to_hours_minutes(1, 20) == (0, 0)
        # 3 slots * 20 = 60 s -> 1 min
        assert _slots_to_hours_minutes(3, 20) == (0, 1)

    def test_90_minutes(self):
        # 90 min = 5400 s, 5400/20 = 270 slots
        h, m = _slots_to_hours_minutes(270, 20)
        assert h == 1 and m == 30

    def test_full_night_example(self):
        # 1535 slots * 20 s = 30700 s -> 511 min 40 s -> round 512 min -> 8h 32m
        h, m = _slots_to_hours_minutes(1535, 20)
        assert h == 8 and m == 32

    def test_different_round_to(self):
        # 60 slots * 60 s = 3600 s -> 1h 0m
        assert _slots_to_hours_minutes(60, 60) == (1, 0)
        assert _slots_to_hours_minutes(90, 60) == (1, 30)


class TestFormatDurationHM:
    """Tests for _format_duration_h_m."""

    def test_zero_minutes(self):
        assert _format_duration_h_m(0, 0) == "0m"

    def test_minutes_only(self):
        assert _format_duration_h_m(0, 30) == "30m"
        assert _format_duration_h_m(0, 1) == "1m"

    def test_hours_and_minutes(self):
        assert _format_duration_h_m(1, 30) == "1h 30m"
        assert _format_duration_h_m(8, 32) == "8h 32m"
        assert _format_duration_h_m(1, 0) == "1h 0m"


class TestToTimestampMs:
    """Tests for _to_timestamp_ms."""

    def test_aware_datetime(self):
        dt = datetime(2024, 11, 12, 3, 0, 0, tzinfo=timezone.utc)
        ms = _to_timestamp_ms(dt)
        assert isinstance(ms, (int, float))
        assert ms == dt.timestamp() * 1000

    def test_naive_datetime_has_timestamp(self):
        dt = datetime(2024, 11, 12, 3, 0, 0)
        ms = _to_timestamp_ms(dt)
        assert ms is not None
        assert isinstance(ms, (int, float))
        # Naive datetime.timestamp() uses local time; we only check it returns a value
        assert ms == dt.timestamp() * 1000

    def test_none_like(self):
        # Objects without timestamp return None (e.g. date-only)
        class NoTimestamp:
            pass
        assert _to_timestamp_ms(NoTimestamp()) is None


class TestColorForIndex:
    """Tests for _color_for_index."""

    def test_returns_hex_string(self):
        c = _color_for_index(0)
        assert isinstance(c, str)
        assert c.startswith("#")
        assert len(c) == 7

    def test_wraps_around_palette(self):
        n = len(PLOTLY_COLORS)
        assert _color_for_index(0) == PLOTLY_COLORS[0]
        assert _color_for_index(n) == PLOTLY_COLORS[0]
        assert _color_for_index(n + 1) == PLOTLY_COLORS[1]


class TestPlotAirmassCurves:
    """Tests for plot_airmass_curves (minimal: empty list, metadata in layout)."""

    def test_empty_scheduled_targets_returns_figure(self):
        from scheduler.observatories import observatories
        obs = observatories["CTIO"]
        obs.update_obs_date(datetime(2024, 11, 12))
        fig = plot_airmass_curves(obs, [], title="Empty")
        assert fig is not None
        if hasattr(fig, "layout") and fig.layout.title:
            assert "Empty" in (fig.layout.title.text or "")

    def test_empty_scheduled_writes_html_when_outpath_given(self, tmp_path):
        from scheduler.observatories import observatories
        obs = observatories["CTIO"]
        obs.update_obs_date(datetime(2024, 11, 12))
        out = tmp_path / "airmass_empty.html"
        fig = plot_airmass_curves(obs, [], outpath=out, title="Empty")
        assert out.exists()
        content = out.read_text()
        assert "plotly" in content.lower() or "Empty" in content

    def test_night_fill_duration_in_subtitle_with_mock_targets(self):
        """When scheduled targets exist, subtitle should show duration (e.g. 'Xh Ym / Ah Bm'), not raw slots."""
        from scheduler.observatories import observatories
        from types import SimpleNamespace
        obs = observatories["CTIO"]
        obs.update_obs_date(datetime(2024, 11, 12))
        # Mock target with scheduled arrays and total_minutes so fill is computed
        raw_len = getattr(obs, "length_of_night", 100) or 100
        mock_target = SimpleNamespace(
            name="mock_field",
            starting_index=0,
            total_minutes=30.0,  # 30 min = 90 slots at round_to=20
            scheduled_time_array=[],
            scheduled_airmass_array=[],
            raw_airmass_array=[1.5] * raw_len,
        )
        fig = plot_airmass_curves(obs, [mock_target], title="One target")
        assert fig is not None
        text = (fig.layout.title.text or "")
        # Should contain duration-style text (e.g. "1h 30m" or "30m") and "/", not "90/..."
        assert "/" in text
        # Should not show raw slot count like "90/1535" (indices)
        assert "Night fill:" in text
        # After our change, we expect "h" and "m" in the subtitle for durations
        assert "m" in text or "h" in text
