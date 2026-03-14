"""
Plot airmass curves for the night with scheduled segments highlighted (Plotly).
"""
import numpy as np
from datetime import timezone, datetime

from scheduler.common import constants as C
from scheduler.common.logging import get_logger
from scheduler.common.utilities import UTC_Offset

logger = get_logger(__name__)

# Long color sequence for many targets (Plotly qualitative + repeated)
PLOTLY_COLORS = (
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
)


def _color_for_index(i: int) -> str:
    return PLOTLY_COLORS[i % len(PLOTLY_COLORS)]


def _to_timestamp_ms(dt):
    """Convert datetime to milliseconds since epoch (UTC) for Plotly."""
    if hasattr(dt, "timestamp"):
        return dt.timestamp() * 1000
    return None


def _slots_to_hours_minutes(slots, round_to_seconds):
    """Convert number of time-array slots to (hours, minutes) for display."""
    total_seconds = slots * round_to_seconds
    total_minutes = int(round(total_seconds / 60.0))
    h, m = divmod(total_minutes, 60)
    return h, m


def _format_duration_h_m(h, m):
    """Format duration as 'Xh Ym' (omit hours if 0)."""
    if h == 0:
        return f"{m}m"
    return f"{h}h {m}m"


def plot_airmass_curves(observatory, scheduled_targets, outpath=None, title=None):
    """
    Build a Plotly figure of airmass vs time for the night.

    Only the segment when each target is observed is plotted (no full-night curve).
    Bottom x-axis: observatory local time. Top x-axis: UTC (parallel, same scale).
    Shows sunset/sunrise/12° twilight lines and metadata: target count, night fill
    as duration (e.g. "4h 30m / 8h 32m (53%)"), and date. Airmass axis fixed 1--3.5.
    """
    try:
        import plotly.graph_objects as go
    except ImportError as e:
        raise ImportError("plotly is required for --plot-airmass; install with: pip install plotly") from e

    fig = go.Figure()
    n = len(scheduled_targets)
    if n == 0:
        logger.warning("No scheduled targets for airmass plot")
        if title:
            fig.update_layout(title=title)
        if outpath:
            fig.write_html(outpath)
        return fig

    # Use observation time order so trace/legend match schedule
    scheduled_targets = sorted(scheduled_targets, key=lambda t: getattr(t, "starting_index", 0))

    AM_MIN, AM_MAX = 1.0, 3.5
    show_legend = n <= 80
    length_of_night = getattr(observatory, "length_of_night", 0) or (
        len(scheduled_targets[0].raw_airmass_array) if scheduled_targets else 0
    )
    filled_slots = 0
    utc_tz = timezone.utc
    local_tz = UTC_Offset(observatory.utc_offset, getattr(observatory, "utc_offset_name", ""))

    # Build traces using UTC timestamps (ms) for x so we can show local + UTC axes
    all_ts = []
    for i, tgt in enumerate(scheduled_targets):
        t_sched = np.asarray(tgt.scheduled_time_array)
        am_sched = np.clip(np.asarray(tgt.scheduled_airmass_array, dtype=float), AM_MIN, AM_MAX)
        if len(t_sched) == 0:
            continue
        start = int(getattr(tgt, "starting_index", 0))
        end = int(start + tgt.total_minutes * 60.0 / C.round_to)
        filled_slots += (end - start)
        # Convert local times to UTC timestamps (ms) for Plotly
        x_ts = []
        for t in t_sched:
            if hasattr(t, "replace") and hasattr(t, "astimezone"):
                if getattr(t, "tzinfo", None) is None:
                    t = t.replace(tzinfo=local_tz)
                t_utc = t.astimezone(utc_tz)
                ts_ms = _to_timestamp_ms(t_utc)
                if ts_ms is not None:
                    x_ts.append(ts_ms)
                    all_ts.append(ts_ms)
            else:
                ts_ms = _to_timestamp_ms(t)
                if ts_ms is not None:
                    x_ts.append(ts_ms)
                    all_ts.append(ts_ms)
        if not x_ts:
            continue
        color = _color_for_index(i)
        fig.add_trace(
            go.Scatter(
                x=x_ts,
                y=am_sched.tolist(),
                mode="lines",
                line=dict(color=color, width=2.5),
                name=tgt.name,
                legendgroup=tgt.name,
                showlegend=show_legend,
            )
        )

    fill_pct = (100.0 * filled_slots / length_of_night) if length_of_night else 0
    filled_h, filled_m = _slots_to_hours_minutes(filled_slots, C.round_to)
    total_h, total_m = _slots_to_hours_minutes(length_of_night, C.round_to)
    fill_duration_str = _format_duration_h_m(filled_h, filled_m)
    total_duration_str = _format_duration_h_m(total_h, total_m)
    obs_date = getattr(observatory, "obs_date", None)
    obs_date_str = obs_date.strftime("%Y-%m-%d") if obs_date and hasattr(obs_date, "strftime") else str(obs_date) if obs_date else ""
    metadata_lines = [
        f"Targets: {n}",
        f"Night fill: {fill_duration_str} / {total_duration_str} ({fill_pct:.1f}%)",
        f"Date: {obs_date_str}",
    ]
    if all_ts:
        ts_min, ts_max = min(all_ts), max(all_ts)
    else:
        ts_min = ts_max = None

    # Vertical lines: sunset, 12° evening/morning twilight, sunrise (use UTC timestamp for x)
    twilight_getter = getattr(observatory, "get_sun_twilight_local_times", None)
    if callable(twilight_getter):
        try:
            tw = twilight_getter()
            for key, label in [
                ("sunset", "Sunset"),
                ("evening_twilight_12", "12° twilight (eve)"),
                ("morning_twilight_12", "12° twilight (morn)"),
                ("sunrise", "Sunrise"),
            ]:
                t = tw.get(key)
                if t is not None:
                    if getattr(t, "tzinfo", None) is None:
                        t = t.replace(tzinfo=local_tz)
                    t_utc = t.astimezone(utc_tz)
                    x_ts = _to_timestamp_ms(t_utc)
                    if x_ts is not None:
                        fig.add_vline(
                            x=x_ts,
                            line_dash="dash",
                            line_color="gray",
                            opacity=0.9,
                            annotation_text=label,
                            annotation_position="top",
                        )
                    if ts_min is not None and x_ts is not None:
                        all_ts.append(x_ts)
        except Exception as e:
            logger.debug("Could not add twilight lines: %s", e)

    if all_ts:
        ts_min, ts_max = min(all_ts), max(all_ts)
    else:
        ts_min = ts_max = None

    # Tick labels: bottom = observatory local, top = UTC
    n_ticks = 8
    if ts_min is not None and ts_max is not None and ts_max > ts_min:
        tickvals = [ts_min + (ts_max - ts_min) * i / (n_ticks - 1) for i in range(n_ticks)]
        tick_local = []
        tick_utc = []
        for ts_ms in tickvals:
            utc_dt = datetime.fromtimestamp(ts_ms / 1000.0, tz=utc_tz)
            local_dt = utc_dt.astimezone(local_tz)
            tick_local.append(local_dt.strftime("%H:%M\n%m-%d"))
            tick_utc.append(utc_dt.strftime("%H:%M\n%m-%d"))
    else:
        tickvals = tick_local = tick_utc = None

    # Subtitle with metadata
    subtitle = "  |  ".join(metadata_lines)
    full_title = title or "Airmass at time of observation"
    if subtitle:
        full_title = f"{full_title}<br><sub>{subtitle}</sub>"

    xaxis_cfg = dict(
        title="Local time",
        type="date",
        range=[ts_min, ts_max] if ts_min is not None else None,
    )
    if tickvals and tick_local:
        xaxis_cfg["tickvals"] = tickvals
        xaxis_cfg["ticktext"] = tick_local
    xaxis2_cfg = dict(
        title="UTC",
        type="date",
        overlaying="x",
        side="top",
        anchor="y",
        range=[ts_min, ts_max] if ts_min is not None else None,
        showgrid=False,
    )
    if tickvals and tick_utc:
        xaxis2_cfg["tickvals"] = tickvals
        xaxis2_cfg["ticktext"] = tick_utc

    fig.update_layout(
        title=dict(text=full_title, font=dict(size=14)),
        xaxis=xaxis_cfg,
        xaxis2=xaxis2_cfg,
        yaxis=dict(
            title="Airmass",
            range=[AM_MAX, AM_MIN],
            fixedrange=True,
            autorange=False,
        ),
        height=600,
        template="plotly_white",
        legend=dict(
            title="Target (observation order)",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
        ),
        margin=dict(r=200, t=120),
    )

    if outpath:
        fig.write_html(outpath)
        logger.info("Wrote airmass plot: %s", outpath)
    return fig
