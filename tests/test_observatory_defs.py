"""Unit tests for observatory definitions."""

import pytest
from scheduler.observatories import observatories


class TestObservatoryDefs:
    """Tests for observatories registry."""

    def test_observatories_dict(self):
        assert "LCO" in observatories
        assert "Lick" in observatories
        assert "Thacher" in observatories
        assert "CTIO" in observatories

    def test_lco_has_swope(self):
        lco = observatories["LCO"]
        assert "Swope" in lco.telescopes

    def test_ctio_has_t80s_and_blanco(self):
        ctio = observatories["CTIO"]
        assert "T80S" in ctio.telescopes
        assert "Blanco" in ctio.telescopes

    def test_observatory_has_ephemeris_after_update(self):
        from datetime import datetime
        obs = observatories["LCO"]
        obs.update_obs_date(datetime(2024, 6, 15))
        assert obs.ephemeris is not None
        assert obs.utc_begin_night is not None
        assert obs.utc_end_night is not None
        assert len(obs.sidereal_radian_array) > 0

    def test_get_sun_twilight_local_times_ordering(self):
        """Twilight/sun times should be ordered: sunset < evening_12 < morning_12 < sunrise."""
        from datetime import datetime
        obs = observatories["CTIO"]
        obs.update_obs_date(datetime(2024, 11, 12))
        tw = obs.get_sun_twilight_local_times()
        assert "sunset" in tw
        assert "evening_twilight_12" in tw
        assert "morning_twilight_12" in tw
        assert "sunrise" in tw
        sunset = tw["sunset"]
        eve12 = tw["evening_twilight_12"]
        morn12 = tw["morning_twilight_12"]
        sunrise = tw["sunrise"]
        assert sunset < eve12, "Sunset should be before 12° evening twilight"
        assert eve12 < morn12, "Evening 12° should be before morning 12°"
        assert morn12 < sunrise, "Morning 12° should be before sunrise"
        # Sunset should be same local calendar day as evening twilight (night of obs_date)
        assert sunset.date() == eve12.date() or (eve12.date() - sunset.date()).days <= 1
