#!/usr/bin/env python3
"""Run Bayesian MCMC inference on the galaxy rotation curve."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Allow running without installation
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from galaxy_mcmc.io import load_rotation_curve
from galaxy_mcmc.likelihood import RotationCurvePosterior
from galaxy_mcmc.mcmc import run_ensemble
from galaxy_mcmc.model import GalaxyPotential
from galaxy_mcmc.plotting import (
    plot_corner,
    plot_enclosed_mass,
    plot_residuals,
    plot_rotation_curve_fit,
    plot_traces,
    print_diagnostics,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Bayesian mass-parameter estimation from a galaxy rotation curve."
    )
    p.add_argument(
        "--data",
        type=Path,
        default=ROOT / "data" / "RadialVelocities.dat",
        help="Path to two-column R [kpc], V [km/s] data file",
    )
    p.add_argument("--outdir", type=Path, default=ROOT / "results", help="Output directory")
    p.add_argument("--n-steps", type=int, default=25_000, help="Steps per chain")
    p.add_argument("--burn-in", type=int, default=5_000, help="Burn-in steps discarded")
    p.add_argument("--n-chains", type=int, default=4, help="Independent MH chains")
    p.add_argument("--seed", type=int, default=42, help="RNG seed")
    p.add_argument(
        "--fit-ah",
        action="store_true",
        help="Also infer the halo scale length Ah [kpc]",
    )
    p.add_argument(
        "--fit-sigma",
        action="store_true",
        help="Also infer the Gaussian noise scale σ",
    )
    p.add_argument(
        "--sigma",
        type=float,
        default=2.2,
        help="Fixed observational noise (km/s) if --fit-sigma is not set",
    )
    p.add_argument(
        "--quick",
        action="store_true",
        help="Short run for smoke tests (2k steps, 2 chains)",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.quick:
        args.n_steps = 2_000
        args.burn_in = 500
        args.n_chains = 2

    radius, velocity = load_rotation_curve(args.data)
    potential = GalaxyPotential()
    posterior = RotationCurvePosterior(
        radius=radius,
        velocity=velocity,
        potential=potential,
        sigma=args.sigma,
        fit_sigma=args.fit_sigma,
        fit_ah=args.fit_ah,
    )

    free = ", ".join(posterior.param_names)
    print(
        f"Loaded {len(radius)} points | R ∈ [{radius.min():.2f}, {radius.max():.2f}] kpc | "
        f"V ∈ [{velocity.min():.2f}, {velocity.max():.2f}] km/s"
    )
    print(f"Free parameters: {free}")
    print(
        f"Running {args.n_chains} Metropolis–Hastings chains × {args.n_steps} steps "
        f"(burn-in {args.burn_in})..."
    )

    results = run_ensemble(
        posterior,
        n_chains=args.n_chains,
        n_steps=args.n_steps,
        burn_in=args.burn_in,
        seed=args.seed,
    )

    report = print_diagnostics(results)
    print(report)

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "posterior_summary.txt").write_text(report + "\n")

    import numpy as np

    combined = np.concatenate([r.samples for r in results], axis=0)
    names = results[0].param_names
    summary = {
        name: {
            "p16": float(np.percentile(combined[:, i], 16)),
            "p50": float(np.percentile(combined[:, i], 50)),
            "p84": float(np.percentile(combined[:, i], 84)),
        }
        for i, name in enumerate(names)
    }
    summary["acceptance_rates"] = [r.acceptance_rate for r in results]
    summary["n_steps"] = args.n_steps
    summary["burn_in"] = args.burn_in
    summary["n_chains"] = args.n_chains
    summary["seed"] = args.seed
    summary["param_names"] = list(names)
    summary["reproducibility"] = (
        f"python scripts/run_inference.py --seed {args.seed} "
        f"--n-steps {args.n_steps} --burn-in {args.burn_in} "
        f"--n-chains {args.n_chains}"
        + (" --fit-ah" if args.fit_ah else "")
        + (" --fit-sigma" if args.fit_sigma else "")
    )
    (args.outdir / "posterior_summary.json").write_text(json.dumps(summary, indent=2))
    np.save(args.outdir / "posterior_samples.npy", combined)

    # Human-readable chain export
    header = ",".join(names)
    np.savetxt(
        args.outdir / "posterior_samples.csv",
        combined,
        delimiter=",",
        header=header,
        comments="",
    )

    plot_rotation_curve_fit(
        radius,
        velocity,
        combined,
        potential,
        args.outdir / "rotation_curve_fit.png",
        param_names=names,
    )
    plot_residuals(
        radius,
        velocity,
        combined,
        potential,
        args.outdir / "residuals.png",
        param_names=names,
    )
    plot_enclosed_mass(
        radius,
        combined,
        potential,
        args.outdir / "enclosed_mass.png",
        param_names=names,
    )
    plot_traces(results, args.outdir / "mcmc_traces.png")
    plot_corner(combined, names, args.outdir / "posterior_corner.png")

    print(f"\nWrote figures and summaries to {args.outdir}/")
    print(f"Reproducibility: {summary['reproducibility']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
