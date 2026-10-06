"""Data loading helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np

DEFAULT_DATA = (
    Path(__file__).resolve().parents[2] / "data" / "RadialVelocities.dat"
)


def load_rotation_curve(
    path: str | Path | None = None, *, sort: bool = True
) -> tuple[np.ndarray, np.ndarray]:
    """Load galactocentric radius (kpc) and circular velocity (km/s).

    The file may start with a comment line. Rows are ``R  V``.
    """
    path = Path(path) if path is not None else DEFAULT_DATA
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected two-column rotation-curve file, got {data.shape}")
    radius = data[:, 0]
    velocity = data[:, 1]
    if sort:
        order = np.argsort(radius)
        radius = radius[order]
        velocity = velocity[order]
    return radius, velocity
