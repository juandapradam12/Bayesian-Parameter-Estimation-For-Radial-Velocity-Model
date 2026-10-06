#!/usr/bin/env python3
"""Mock-data recovery check for the Metropolis–Hastings sampler.

Generates a synthetic rotation curve from known (Mb, Md, Mh, Ah), runs MCMC,
and prints whether the true values fall inside the 16–84% posterior interval.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from galaxy_mcmc.likelihood import RotationCurvePosterior
from galaxy_mcmc.mcmc import run_ensemble
from galaxy_mcmc.model import GalaxyPotential


TRUE = {"Mb": 50.0, "Md": 1.4e4, "Mh": 2.6e4, "Ah": 64.3}


def make_mock(n: int = 80, sigma: float = 1.5, seed: int = 0):
    rng = np.random.default_rng(seed)
    radius = np.linspace(1.0, 250.0, n)
    pot = GalaxyPotential(Ah=TRUE["Ah"])
    truth = pot.circular_velocity(radius, TRUE["Mb"], TRUE["Md"], TRUE["Mh"])
    velocity = truth + sigma * rng.normal(size=n)
    return radius, velocity, sigma


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-steps", type=int, default=12_000)
    p.add_argument("--burn-in", type=int, default=3_000)
    p.add_argument("--n-chains", type=int, default=2)
    p.add_argument("--seed", type=int, default=7)
    args = p.parse_args()

    radius, velocity, sigma = make_mock(seed=args.seed)
    posterior = RotationCurvePosterior(
        radius=radius,
        velocity=velocity,
        potential=GalaxyPotential(),
        sigma=sigma,
        fit_ah=True,
        fit_sigma=False,
    )
    results = run_ensemble(
        posterior,
        n_chains=args.n_chains,
        n_steps=args.n_steps,
        burn_in=args.burn_in,
        seed=args.seed,
    )
    samples = np.concatenate([r.samples for r in results], axis=0)
    names = results[0].param_names

    print("Mock recovery (true vs posterior 16/50/84%)")
    ok = True
    for i, name in enumerate(names):
        truth = TRUE[name]
        p16, p50, p84 = np.percentile(samples[:, i], [16, 50, 84])
        inside = p16 <= truth <= p84
        ok = ok and inside
        flag = "OK" if inside else "MISS"
        print(f"  {name}: true={truth:.4g}  median={p50:.4g}  [{p16:.4g}, {p84:.4g}]  {flag}")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
