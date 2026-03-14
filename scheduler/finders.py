"""
Generate 2MASS-J finder charts for target lists.

Generic for any observatory/telescope; target list format: field, ra, dec, epoch.
Uses Plotly for interactive-quality static PNG export.
"""
import os
import warnings

import numpy as np
from astroquery.skyview import SkyView
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.ndimage import gaussian_filter
from photutils.background import Background2D, MedianBackground
from astropy.stats import SigmaClip
import plotly.graph_objects as go

from scheduler.common.logging import get_logger

warnings.filterwarnings("ignore")
logger = get_logger(__name__)


def make_finders(target_list, directory="finders", clobber=False):
    """
    Create 2MASS-J finder chart PNGs for each target in target_list.

    Parameters
    ----------
    target_list : str or path-like
        File with columns: field, ra, dec, epoch (space or comma separated).
    directory : str, optional
        Output directory for PNGs (created if missing). Default "finders".
    clobber : bool, optional
        If False, skip targets that already have a finder file. Default False.

    Returns
    -------
    int
        Number of finder files written.
    """
    sv = SkyView()
    os.makedirs(directory, exist_ok=True)
    params = {
        "radius": 14.0 * u.arcmin,
        "pixels": "1000",
        "deedger": "skyview.process.Deedger",
        "survey": ["2MASS-J"],
    }
    target_table = ascii.read(target_list, names=("field", "ra", "dec", "epoch"))
    n_saved = 0
    for row in target_table:
        coord = SkyCoord(row["ra"], row["dec"], unit=(u.hour, u.deg))
        rah, decd = coord.to_string(style="hmsdms", sep=":", precision=3).split()
        name = row["field"]
        outpath = os.path.join(directory, f"{name}_finder.png")
        if os.path.exists(outpath) and not clobber:
            logger.warning("%s exists and clobber=False. Skipping.", outpath)
            continue
        images = sv.get_images(position=coord, **params)
        image = images[0][0]
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(
            image.data, (9, 9), filter_size=(3, 3),
            sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
        )
        image.data = image.data - bkg.background
        smoothed = gaussian_filter(image.data, 1.3)
        vmin = float(np.percentile(image.data.flatten(), 10))
        vmax = float(np.percentile(image.data.flatten(), 99.0))
        ny, nx = smoothed.shape

        # Plotly Heatmap: first row of z is at bottom (y=0), matching "lower" origin
        fig = go.Figure(
            data=go.Heatmap(
                z=smoothed,
                colorscale="Greys_r",
                zmin=vmin,
                zmax=vmax,
                showscale=False,
            )
        )

        # Compass: N/E lines in pixel coords (lower-right)
        n, m = nx, ny
        fig.add_shape(
            type="line",
            x0=n - 200, y0=10, x1=n - 10, y1=10,
            line=dict(color="red", width=2),
        )
        fig.add_shape(
            type="line",
            x0=n - 10, y0=10, x1=n - 10, y1=200,
            line=dict(color="red", width=2),
        )

        # N and E labels
        fig.add_annotation(x=n - 20, y=200, text="N", showarrow=False, font=dict(color="red", size=14))
        fig.add_annotation(x=n - 200, y=20, text="E", showarrow=False, font=dict(color="red", size=14))

        # Target name and coords (paper coords, right side)
        fig.add_annotation(
            xref="paper", yref="paper", x=0.98, y=0.95, xanchor="left",
            text=name, showarrow=False, font=dict(size=14, color="black"),
        )
        fig.add_annotation(
            xref="paper", yref="paper", x=0.98, y=0.90, xanchor="left",
            text=f"{rah}  {decd}", showarrow=False, font=dict(size=12, color="black"),
        )
        fig.add_annotation(
            xref="paper", yref="paper", x=0.98, y=0.85, xanchor="left",
            text="2MASS J-band", showarrow=False, font=dict(size=12, color="black"),
        )

        fig.update_layout(
            xaxis=dict(title="RA (J2000)", constrain="domain", scaleanchor="y"),
            yaxis=dict(title="Dec (J2000)"),
            margin=dict(l=60, r=120, t=40, b=60),
            width=1200,
            height=900,
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="white",
        )

        fig.write_image(outpath, scale=1)
        logger.info("Saved to %s", outpath)
        n_saved += 1
    return n_saved


def main():
    """CLI: generate 2MASS-J finder charts for a target list (field, ra, dec, epoch)."""
    import argparse
    ap = argparse.ArgumentParser(description="Generate 2MASS-J finder charts for a target list.")
    ap.add_argument("target_list", help="File with columns: field, ra, dec, epoch")
    ap.add_argument("-d", "--directory", default="finders", help="Output directory")
    ap.add_argument("--clobber", action="store_true", help="Overwrite existing finder files")
    args = ap.parse_args()
    make_finders(args.target_list, directory=args.directory, clobber=args.clobber)
