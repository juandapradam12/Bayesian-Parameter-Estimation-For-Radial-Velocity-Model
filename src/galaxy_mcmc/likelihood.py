"""Likelihood and prior for Bayesian rotation-curve inference."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import GalaxyPotential, PARAM_NAMES


@dataclass
class RotationCurvePosterior:
    """Log-posterior for Gaussian noise on the observed rotation curve.

    Parameters are sampled in log-space: ``theta = log(Mb, Md, Mh)``.
    An optional free log-noise parameter can be appended as a fourth
    element of ``theta`` when ``fit_sigma`` is True.
    """

    radius: np.ndarray
    velocity: np.ndarray
    potential: GalaxyPotential
    sigma: float = 2.0
    fit_sigma: bool = False
    log_mass_bounds: tuple[float, float] = (np.log(1e-8), np.log(1e9))
    log_sigma_bounds: tuple[float, float] = (np.log(0.1), np.log(50.0))

    def __post_init__(self) -> None:
        self.radius = np.asarray(self.radius, dtype=float)
        self.velocity = np.asarray(self.velocity, dtype=float)
        if self.radius.shape != self.velocity.shape:
            raise ValueError("radius and velocity must have the same shape")

    @property
    def n_params(self) -> int:
        return 4 if self.fit_sigma else 3

    @property
    def param_names(self) -> tuple[str, ...]:
        return PARAM_NAMES + (("sigma",) if self.fit_sigma else ())

    def unpack(self, theta: np.ndarray) -> tuple[np.ndarray, float]:
        masses = np.exp(theta[:3])
        sigma = float(np.exp(theta[3])) if self.fit_sigma else self.sigma
        return masses, sigma

    def log_prior(self, theta: np.ndarray) -> float:
        lo, hi = self.log_mass_bounds
        if np.any(theta[:3] < lo) or np.any(theta[:3] > hi):
            return -np.inf
        if self.fit_sigma:
            slo, shi = self.log_sigma_bounds
            if theta[3] < slo or theta[3] > shi:
                return -np.inf
        return 0.0

    def log_likelihood(self, theta: np.ndarray) -> float:
        masses, sigma = self.unpack(theta)
        pred = self.potential.from_theta(self.radius, masses)
        resid = (self.velocity - pred) / sigma
        n = self.velocity.size
        return -0.5 * np.sum(resid**2) - n * np.log(sigma) - 0.5 * n * np.log(2 * np.pi)

    def log_posterior(self, theta: np.ndarray) -> float:
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta)

    def chi_squared(self, theta: np.ndarray) -> float:
        masses, sigma = self.unpack(theta)
        pred = self.potential.from_theta(self.radius, masses)
        return float(np.sum(((self.velocity - pred) / sigma) ** 2))
