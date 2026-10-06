"""Likelihood and prior for Bayesian rotation-curve inference."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import GalaxyPotential, PARAM_NAMES, SCALE_PARAMS


@dataclass
class RotationCurvePosterior:
    """Log-posterior for Gaussian noise on the observed rotation curve.

    Parameters are sampled in log-space. Base vector is
    ``theta = log(Mb, Md, Mh)``. Optional extras, in order:

    - ``fit_ah`` → append ``log(Ah)``
    - ``fit_sigma`` → append ``log(sigma)``
    """

    radius: np.ndarray
    velocity: np.ndarray
    potential: GalaxyPotential
    sigma: float = 2.0
    fit_sigma: bool = False
    fit_ah: bool = False
    log_mass_bounds: tuple[float, float] = (np.log(1e-8), np.log(1e9))
    log_ah_bounds: tuple[float, float] = (np.log(5.0), np.log(300.0))
    log_sigma_bounds: tuple[float, float] = (np.log(0.1), np.log(50.0))

    def __post_init__(self) -> None:
        self.radius = np.asarray(self.radius, dtype=float)
        self.velocity = np.asarray(self.velocity, dtype=float)
        if self.radius.shape != self.velocity.shape:
            raise ValueError("radius and velocity must have the same shape")

    @property
    def n_params(self) -> int:
        return 3 + int(self.fit_ah) + int(self.fit_sigma)

    @property
    def param_names(self) -> tuple[str, ...]:
        names: tuple[str, ...] = PARAM_NAMES
        if self.fit_ah:
            names = names + ("Ah",)
        if self.fit_sigma:
            names = names + ("sigma",)
        return names

    def unpack(self, theta: np.ndarray) -> tuple[np.ndarray, float, float]:
        """Return ``(masses, Ah, sigma)`` in physical units from log-theta."""
        masses = np.exp(theta[:3])
        idx = 3
        if self.fit_ah:
            ah = float(np.exp(theta[idx]))
            idx += 1
        else:
            ah = float(self.potential.Ah)
        sigma = float(np.exp(theta[idx])) if self.fit_sigma else self.sigma
        return masses, ah, sigma

    def predict(self, radius: np.ndarray, theta: np.ndarray) -> np.ndarray:
        masses, ah, _ = self.unpack(theta)
        return self.potential.from_theta(radius, masses, Ah=ah)

    def default_start(self) -> np.ndarray:
        """Reasonable log-space starting point for MH."""
        start = [np.log(1.0), np.log(1.4e4), np.log(2.6e4)]
        if self.fit_ah:
            start.append(np.log(SCALE_PARAMS["Ah"]))
        if self.fit_sigma:
            start.append(np.log(self.sigma))
        return np.asarray(start, dtype=float)

    def proposal_scale(self) -> np.ndarray:
        """Default relative MH step sizes in log-space."""
        scale = [0.8, 0.04, 0.04]
        if self.fit_ah:
            scale.append(0.05)
        if self.fit_sigma:
            scale.append(0.05)
        return np.asarray(scale, dtype=float)

    def log_prior(self, theta: np.ndarray) -> float:
        lo, hi = self.log_mass_bounds
        if np.any(theta[:3] < lo) or np.any(theta[:3] > hi):
            return -np.inf
        idx = 3
        if self.fit_ah:
            alo, ahi = self.log_ah_bounds
            if theta[idx] < alo or theta[idx] > ahi:
                return -np.inf
            idx += 1
        if self.fit_sigma:
            slo, shi = self.log_sigma_bounds
            if theta[idx] < slo or theta[idx] > shi:
                return -np.inf
        return 0.0

    def log_likelihood(self, theta: np.ndarray) -> float:
        masses, ah, sigma = self.unpack(theta)
        pred = self.potential.from_theta(self.radius, masses, Ah=ah)
        resid = (self.velocity - pred) / sigma
        n = self.velocity.size
        return -0.5 * np.sum(resid**2) - n * np.log(sigma) - 0.5 * n * np.log(2 * np.pi)

    def log_posterior(self, theta: np.ndarray) -> float:
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta)

    def chi_squared(self, theta: np.ndarray) -> float:
        masses, ah, sigma = self.unpack(theta)
        pred = self.potential.from_theta(self.radius, masses, Ah=ah)
        return float(np.sum(((self.velocity - pred) / sigma) ** 2))
