"""Three-component galactic mass model for rotation curves.

The circular velocity is built from a Plummer bulge, a Miyamoto–Nagai disk,
and an Allen–Santillán-like isothermal halo. Scale lengths are held fixed;
the free parameters are the scaled masses ``Mb``, ``Md`` and ``Mh``
(gravitational constant absorbed into the units).

Physically,

.. math::

    v_c(R) = \\sqrt{v_b^2(R) + v_d^2(R) + v_h^2(R)}

with

.. math::

    v_b^2 &= M_b\\, R^2 / (R^2 + B_b^2)^{3/2}, \\\\
    v_d^2 &= M_d\\, R^2 / (R^2 + (B_d + A_d)^2)^{3/2}, \\\\
    v_h^2 &= M_h / (R^2 + A_h^2)^{1/2}.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

# Fixed geometric scales (kpc), matching the original homework constants
# (Irrgang / Allen–Santillán style galactic potential).
SCALE_PARAMS: Mapping[str, float] = {
    "Bb": 0.2497,
    "Bd": 5.16,
    "Ad": 0.3105,
    "Ah": 64.3,
}

PARAM_NAMES = ("Mb", "Md", "Mh")


@dataclass(frozen=True)
class GalaxyPotential:
    """Analytic bulge + disk + halo rotation-curve model.

    Parameters
    ----------
    Bb, Bd, Ad, Ah:
        Geometric scale lengths in kiloparsecs.
    """

    Bb: float = SCALE_PARAMS["Bb"]
    Bd: float = SCALE_PARAMS["Bd"]
    Ad: float = SCALE_PARAMS["Ad"]
    Ah: float = SCALE_PARAMS["Ah"]

    def component_velocities_squared(
        self,
        radius: np.ndarray,
        Mb: float,
        Md: float,
        Mh: float,
        Ah: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(vb², vd², vh²)`` at each radius."""
        r = np.asarray(radius, dtype=float)
        r2 = r * r
        ah = self.Ah if Ah is None else float(Ah)
        vb2 = Mb * r2 / np.power(r2 + self.Bb**2, 1.5)
        vd2 = Md * r2 / np.power(r2 + (self.Bd + self.Ad) ** 2, 1.5)
        vh2 = Mh / np.sqrt(r2 + ah**2)
        return vb2, vd2, vh2

    def circular_velocity(
        self,
        radius: np.ndarray,
        Mb: float,
        Md: float,
        Mh: float,
        Ah: float | None = None,
    ) -> np.ndarray:
        """Circular speed ``vc(R)`` in km/s for scaled masses Mb, Md, Mh."""
        vb2, vd2, vh2 = self.component_velocities_squared(radius, Mb, Md, Mh, Ah=Ah)
        return np.sqrt(np.maximum(vb2 + vd2 + vh2, 0.0))

    def components(
        self,
        radius: np.ndarray,
        Mb: float,
        Md: float,
        Mh: float,
        Ah: float | None = None,
    ) -> dict[str, np.ndarray]:
        """Named circular-speed contributions (not summed in quadrature)."""
        vb2, vd2, vh2 = self.component_velocities_squared(radius, Mb, Md, Mh, Ah=Ah)
        return {
            "bulge": np.sqrt(np.maximum(vb2, 0.0)),
            "disk": np.sqrt(np.maximum(vd2, 0.0)),
            "halo": np.sqrt(np.maximum(vh2, 0.0)),
            "total": np.sqrt(np.maximum(vb2 + vd2 + vh2, 0.0)),
        }

    def from_theta(
        self,
        radius: np.ndarray,
        theta: np.ndarray,
        Ah: float | None = None,
    ) -> np.ndarray:
        """Evaluate ``vc`` from a parameter vector ``[Mb, Md, Mh]``."""
        Mb, Md, Mh = theta
        return self.circular_velocity(radius, Mb, Md, Mh, Ah=Ah)
