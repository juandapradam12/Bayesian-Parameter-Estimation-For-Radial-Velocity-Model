"""Unit tests for the galactic potential and MCMC utilities."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from galaxy_mcmc.io import load_rotation_curve
from galaxy_mcmc.likelihood import RotationCurvePosterior
from galaxy_mcmc.mcmc import MetropolisHastings, gelman_rubin, run_ensemble
from galaxy_mcmc.model import GalaxyPotential


def test_load_data():
    r, v = load_rotation_curve()
    assert r.ndim == 1 and v.ndim == 1
    assert len(r) == len(v) == 300
    assert np.all(np.diff(r) >= 0)
    assert r.min() > 0


def test_circular_velocity_positive_and_finite():
    pot = GalaxyPotential()
    r = np.linspace(0.5, 100, 50)
    v = pot.circular_velocity(r, Mb=100.0, Md=1.0e4, Mh=2.0e4)
    assert np.all(np.isfinite(v))
    assert np.all(v >= 0)
    # Halo should dominate at large R for these masses
    comps = pot.components(r, 100.0, 1.0e4, 2.0e4)
    assert comps["halo"][-1] > comps["bulge"][-1]


def test_velocity_increases_with_mass():
    pot = GalaxyPotential()
    r = np.array([10.0])
    v1 = pot.circular_velocity(r, 1e3, 1e4, 1e4)[0]
    v2 = pot.circular_velocity(r, 1e3, 2e4, 1e4)[0]
    assert v2 > v1


def test_log_posterior_finite_at_reasonable_point():
    r, v = load_rotation_curve()
    post = RotationCurvePosterior(r, v, GalaxyPotential(), sigma=2.2)
    theta = np.log([1.0, 1.4e4, 2.6e4])
    lp = post.log_posterior(theta)
    assert np.isfinite(lp)


def test_log_prior_rejects_out_of_bounds():
    r, v = load_rotation_curve()
    post = RotationCurvePosterior(r, v, GalaxyPotential())
    assert post.log_prior(np.array([100.0, 0.0, 0.0])) == -np.inf


def test_metropolis_hastings_runs():
    r, v = load_rotation_curve()
    post = RotationCurvePosterior(r, v, GalaxyPotential(), sigma=2.2)
    sampler = MetropolisHastings(
        log_prob=post.log_posterior,
        n_params=3,
        proposal_scale=np.array([0.5, 0.05, 0.05]),
        adapt_until=200,
        rng=np.random.default_rng(0),
    )
    result = sampler.run(
        start=np.log([1.0, 1.4e4, 2.6e4]),
        n_steps=500,
        burn_in=100,
        param_names=post.param_names,
    )
    assert result.chain.shape == (500, 3)
    assert 0.0 < result.acceptance_rate < 1.0
    assert result.samples.shape[0] == 400


def test_ensemble_recovers_disk_and_halo():
    r, v = load_rotation_curve()
    post = RotationCurvePosterior(r, v, GalaxyPotential(), sigma=2.2)
    results = run_ensemble(post, n_chains=2, n_steps=3000, burn_in=800, seed=1)
    combined = np.concatenate([res.samples for res in results], axis=0)
    md = np.median(combined[:, 1])
    mh = np.median(combined[:, 2])
    # Rough recovery against the least-squares solution ~1.45e4, 2.62e4
    assert 8e3 < md < 2.5e4
    assert 1.5e4 < mh < 4e4


def test_gelman_rubin_near_one_for_identical_chains():
    rng = np.random.default_rng(0)
    base = rng.normal(size=(1000, 3))
    chains = [base + 0.01 * rng.normal(size=base.shape) for _ in range(3)]
    rhat = gelman_rubin(chains)
    assert np.all(rhat < 1.1)
