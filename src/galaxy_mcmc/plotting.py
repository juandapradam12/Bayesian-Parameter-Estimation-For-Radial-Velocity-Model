"""Plotting utilities for posterior diagnostics and model fits."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .mcmc import MCMCResult, gelman_rubin
from .model import GalaxyPotential


def _decode_sample(
    sample: np.ndarray,
    param_names: tuple[str, ...],
    potential: GalaxyPotential,
) -> tuple[np.ndarray, float]:
    """Return ``(masses, Ah)`` from a physical-space sample row."""
    masses = sample[:3]
    if "Ah" in param_names:
        ah = float(sample[param_names.index("Ah")])
    else:
        ah = float(potential.Ah)
    return masses, ah


def plot_rotation_curve_fit(
    radius: np.ndarray,
    velocity: np.ndarray,
    samples: np.ndarray,
    potential: GalaxyPotential,
    outfile: str | Path,
    param_names: tuple[str, ...] = ("Mb", "Md", "Mh"),
    n_posterior_draws: int = 200,
    seed: int = 0,
) -> Path:
    """Data + median model + posterior predictive band + components."""
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    r_grid = np.linspace(radius.min(), radius.max(), 400)
    idx = rng.choice(len(samples), size=min(n_posterior_draws, len(samples)), replace=False)
    curves = []
    for i in idx:
        masses, ah = _decode_sample(samples[i], param_names, potential)
        curves.append(potential.from_theta(r_grid, masses, Ah=ah))
    curves = np.asarray(curves)
    lo, med, hi = np.percentile(curves, [16, 50, 84], axis=0)
    med_masses, med_ah = _decode_sample(np.median(samples, axis=0), param_names, potential)
    comps = potential.components(r_grid, *med_masses, Ah=med_ah)

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.fill_between(r_grid, lo, hi, color="#4C78A8", alpha=0.25, label="68% posterior band")
    ax.plot(r_grid, med, color="#4C78A8", lw=2.2, label="Median model")
    ax.plot(r_grid, comps["bulge"], ls="--", color="#F58518", lw=1.4, label="Bulge")
    ax.plot(r_grid, comps["disk"], ls="--", color="#54A24B", lw=1.4, label="Disk")
    ax.plot(r_grid, comps["halo"], ls="--", color="#E45756", lw=1.4, label="Halo")
    ax.scatter(radius, velocity, s=14, c="#333333", alpha=0.75, zorder=5, label="Data")
    ax.set_xlabel("Galactocentric radius $R$ [kpc]")
    ax.set_ylabel("Circular velocity $v_c$ [km/s]")
    ax.set_title("Galaxy rotation curve — Bayesian mass decomposition")
    ax.legend(frameon=False, loc="upper right")
    ax.set_xlim(0, radius.max())
    ax.set_ylim(0, max(velocity.max(), hi.max()) * 1.08)
    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    plt.close(fig)
    return outfile


def plot_residuals(
    radius: np.ndarray,
    velocity: np.ndarray,
    samples: np.ndarray,
    potential: GalaxyPotential,
    outfile: str | Path,
    param_names: tuple[str, ...] = ("Mb", "Md", "Mh"),
) -> Path:
    """Residuals (data − median model) versus radius."""
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    med_masses, med_ah = _decode_sample(np.median(samples, axis=0), param_names, potential)
    model = potential.from_theta(radius, med_masses, Ah=med_ah)
    resid = velocity - model

    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    ax.axhline(0.0, color="k", lw=0.8)
    ax.scatter(radius, resid, s=14, c="#333333", alpha=0.75)
    ax.set_xlabel("Galactocentric radius $R$ [kpc]")
    ax.set_ylabel(r"Residual $V_{\mathrm{data}} - V_{\mathrm{model}}$ [km/s]")
    ax.set_title("Fit residuals")
    ax.set_xlim(0, radius.max())
    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    plt.close(fig)
    return outfile


def plot_traces(
    results: list[MCMCResult],
    outfile: str | Path,
) -> Path:
    """Trace plots for all chains and parameters."""
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    n_params = results[0].chain.shape[1]
    names = results[0].param_names
    fig, axes = plt.subplots(n_params, 1, figsize=(9, 2.2 * n_params), sharex=True)
    if n_params == 1:
        axes = [axes]
    colors = plt.cm.tab10(np.linspace(0, 1, len(results)))
    for p in range(n_params):
        for c, res in enumerate(results):
            axes[p].plot(np.exp(res.chain[:, p]), color=colors[c], alpha=0.7, lw=0.6)
            axes[p].axvline(res.burn_in, color="k", ls=":", lw=0.8)
        axes[p].set_ylabel(names[p])
        axes[p].set_yscale("log")
    axes[-1].set_xlabel("Step")
    fig.suptitle("MCMC traces (log-scale physical parameters)", y=1.01)
    fig.tight_layout()
    fig.savefig(outfile, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return outfile


def plot_corner(
    samples: np.ndarray,
    param_names: tuple[str, ...],
    outfile: str | Path,
) -> Path:
    """Simple corner / pairwise posterior plot without external deps."""
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    n = samples.shape[1]
    fig, axes = plt.subplots(n, n, figsize=(2.4 * n, 2.4 * n))
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if i < j:
                ax.axis("off")
                continue
            if i == j:
                ax.hist(samples[:, i], bins=40, color="#4C78A8", alpha=0.85)
                ax.set_yticks([])
            else:
                ax.scatter(
                    samples[:, j],
                    samples[:, i],
                    s=3,
                    alpha=0.25,
                    c="#4C78A8",
                    rasterized=True,
                )
            if i == n - 1:
                ax.set_xlabel(param_names[j])
            else:
                ax.set_xticklabels([])
            if j == 0 and i != 0:
                ax.set_ylabel(param_names[i])
            elif j != 0:
                ax.set_yticklabels([])
    fig.suptitle("Posterior corner plot", y=1.02)
    fig.tight_layout()
    fig.savefig(outfile, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return outfile


def print_diagnostics(results: list[MCMCResult]) -> str:
    """Format acceptance rates, parameter summaries, and R-hat."""
    lines = []
    for i, res in enumerate(results):
        lines.append(f"Chain {i}: acceptance = {res.acceptance_rate:.3f}")
    burned = [res.samples for res in results]
    combined = np.concatenate(burned, axis=0)
    names = results[0].param_names
    lines.append("")
    lines.append("Posterior summary (combined chains):")
    for i, name in enumerate(names):
        p16, p50, p84 = np.percentile(combined[:, i], [16, 50, 84])
        lines.append(
            f"  {name}: {p50:.4g}  [{p16:.4g}, {p84:.4g}]  (16/50/84%)"
        )
    if len(results) >= 2:
        rhat = gelman_rubin([np.log(s) for s in burned])
        lines.append("")
        lines.append("Gelman–Rubin R-hat (on log parameters):")
        for name, r in zip(names, rhat):
            lines.append(f"  {name}: {r:.4f}")
    return "\n".join(lines)
