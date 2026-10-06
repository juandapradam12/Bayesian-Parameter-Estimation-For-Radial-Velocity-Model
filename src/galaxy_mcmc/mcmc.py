"""Metropolis–Hastings MCMC with adaptive Gaussian proposals."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .likelihood import RotationCurvePosterior


LogProbFn = Callable[[np.ndarray], float]


@dataclass
class MCMCResult:
    """Container for MCMC output."""

    chain: np.ndarray  # (n_steps, n_params) in sampling space (log masses)
    log_prob: np.ndarray
    acceptance_rate: float
    burn_in: int
    param_names: tuple[str, ...]
    proposal_scale: np.ndarray

    @property
    def samples(self) -> np.ndarray:
        """Posterior samples after burn-in, transformed to physical space."""
        raw = self.chain[self.burn_in :]
        physical = np.exp(raw)
        return physical

    def summary(self, percentiles: tuple[float, ...] = (16, 50, 84)) -> dict[str, dict[str, float]]:
        samples = self.samples
        out: dict[str, dict[str, float]] = {}
        for i, name in enumerate(self.param_names):
            qs = np.percentile(samples[:, i], percentiles)
            out[name] = {f"p{int(p)}": float(q) for p, q in zip(percentiles, qs)}
            out[name]["mean"] = float(np.mean(samples[:, i]))
            out[name]["std"] = float(np.std(samples[:, i]))
        return out


@dataclass
class MetropolisHastings:
    """Random-walk Metropolis–Hastings sampler with optional adaptation.

    Proposals are independent Gaussians in the sampling coordinates
    (log-masses). During an optional warm-up window the proposal scale is
    adapted toward a target acceptance rate (~0.25 for multi-dimensional
    Gaussian targets).
    """

    log_prob: LogProbFn
    n_params: int
    proposal_scale: np.ndarray | float = 0.1
    target_accept: float = 0.25
    adapt_until: int = 0
    rng: np.random.Generator = field(default_factory=lambda: np.random.default_rng())

    def __post_init__(self) -> None:
        self.proposal_scale = np.broadcast_to(
            np.asarray(self.proposal_scale, dtype=float), self.n_params
        ).copy()

    def run(
        self,
        start: np.ndarray,
        n_steps: int,
        burn_in: int = 0,
        param_names: tuple[str, ...] | None = None,
    ) -> MCMCResult:
        start = np.asarray(start, dtype=float)
        if start.shape != (self.n_params,):
            raise ValueError(f"start must have shape ({self.n_params},)")

        chain = np.empty((n_steps, self.n_params))
        log_prob = np.empty(n_steps)
        current = start.copy()
        current_lp = self.log_prob(current)
        if not np.isfinite(current_lp):
            raise ValueError("Starting point has non-finite log-posterior")

        accepted = 0
        scale = self.proposal_scale.copy()

        for i in range(n_steps):
            proposal = current + scale * self.rng.normal(size=self.n_params)
            proposal_lp = self.log_prob(proposal)
            if np.isfinite(proposal_lp):
                log_alpha = proposal_lp - current_lp
                if np.log(self.rng.random()) < log_alpha:
                    current = proposal
                    current_lp = proposal_lp
                    accepted += 1

            chain[i] = current
            log_prob[i] = current_lp

            if self.adapt_until and i > 0 and i < self.adapt_until and (i + 1) % 100 == 0:
                rate = accepted / (i + 1)
                if rate < self.target_accept * 0.7:
                    scale *= 0.9
                elif rate > self.target_accept * 1.3:
                    scale *= 1.1

        names = param_names or tuple(f"p{i}" for i in range(self.n_params))
        return MCMCResult(
            chain=chain,
            log_prob=log_prob,
            acceptance_rate=accepted / n_steps,
            burn_in=burn_in,
            param_names=names,
            proposal_scale=scale,
        )


def run_ensemble(
    posterior: RotationCurvePosterior,
    n_chains: int = 4,
    n_steps: int = 20_000,
    burn_in: int = 5_000,
    seed: int = 42,
    start_guess: np.ndarray | None = None,
) -> list[MCMCResult]:
    """Run several independent Metropolis–Hastings chains.

    Starting points are jittered around ``start_guess`` (default: a rough
    least-squares-inspired guess) so Gelman–Rubin style comparisons are
    meaningful.
    """
    rng = np.random.default_rng(seed)
    if start_guess is None:
        start_guess = np.array([np.log(1.0), np.log(1.4e4), np.log(2.6e4)])
        if posterior.fit_sigma:
            start_guess = np.append(start_guess, np.log(posterior.sigma))

    results: list[MCMCResult] = []
    # Relative proposal scales: bulge mass is poorly constrained → larger steps
    base_scale = np.array([0.8, 0.04, 0.04])
    if posterior.fit_sigma:
        base_scale = np.append(base_scale, 0.05)

    for c in range(n_chains):
        jitter = 0.15 * base_scale * rng.normal(size=posterior.n_params)
        start = start_guess + jitter
        # ensure prior-valid start
        for _ in range(50):
            if np.isfinite(posterior.log_posterior(start)):
                break
            start = start_guess + 0.15 * base_scale * rng.normal(size=posterior.n_params)
        else:
            start = start_guess.copy()

        sampler = MetropolisHastings(
            log_prob=posterior.log_posterior,
            n_params=posterior.n_params,
            proposal_scale=base_scale,
            adapt_until=max(burn_in, 1000),
            rng=np.random.default_rng(seed + c + 1),
        )
        result = sampler.run(
            start=start,
            n_steps=n_steps,
            burn_in=burn_in,
            param_names=posterior.param_names,
        )
        results.append(result)
    return results


def gelman_rubin(chains: list[np.ndarray]) -> np.ndarray:
    """Compute the Gelman–Rubin ``R-hat`` diagnostic per parameter.

    Each chain array should already be burned-in and have shape
    ``(n_samples, n_params)``.
    """
    if len(chains) < 2:
        raise ValueError("Need at least two chains for Gelman–Rubin")
    stacked = np.stack(chains, axis=0)  # (m, n, p)
    m, n, p = stacked.shape
    chain_means = stacked.mean(axis=1)
    chain_vars = stacked.var(axis=1, ddof=1)
    between = n * chain_means.var(axis=0, ddof=1)
    within = chain_vars.mean(axis=0)
    var_hat = ((n - 1) / n) * within + between / n
    return np.sqrt(var_hat / within)
