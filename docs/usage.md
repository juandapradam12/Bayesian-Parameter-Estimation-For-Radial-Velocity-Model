# Usage guide

How to install, run, and interpret **Galaxy MCMC**. Physics background: [theory.md](theory.md).

## Installation

```bash
python3 -m pip install -r requirements.txt
python3 -m pip install -e .
```

Requires **Python ≥ 3.10**, NumPy, Matplotlib, SciPy. Pytest is used for tests.

## One-command workflows

| Command | What it does |
| --- | --- |
| `make test` | Run the full pytest suite (incl. mock recovery) |
| `make mock` | Standalone synthetic recovery script |
| `make quick` | Short smoke inference (`--quick --fit-ah`) |
| `make run` | Recommended full run (`--fit-ah --fit-sigma`) |
| `cd legacy && make -f Tarea5.mk` | Educational C pipeline |

## Command-line interface

```bash
python3 scripts/run_inference.py --help
```

| Flag | Meaning |
| --- | --- |
| `--data PATH` | Two-column `R [kpc]  V [km/s]` file |
| `--outdir DIR` | Output directory (default `results/`) |
| `--n-steps N` | Steps per chain (default 25000) |
| `--burn-in N` | Discarded warm-up (default 5000) |
| `--n-chains N` | Independent MH chains (default 4) |
| `--fit-ah` | Infer halo scale \(A_h\) |
| `--fit-sigma` | Infer noise \(\sigma\) |
| `--sigma S` | Fixed noise [km/s] if not fitting \(\sigma\) |
| `--quick` | Short smoke run (2k steps, 2 chains) |
| `--seed N` | RNG seed (written into the JSON summary) |

### Recommended first run

```bash
python3 scripts/run_inference.py --fit-ah --fit-sigma
```

## Interpreting outputs

All figures and tables land in `results/` (or `--outdir`).

1. **`posterior_summary.txt` / `.json`**
   - Acceptance rates ~0.15–0.4 are healthy for this proposal.
   - \(\hat{R}\) ≲ 1.05 means chains agree.
   - JSON includes `seed` and a copy-paste `reproducibility` command.
2. **`rotation_curve_fit.png`** — science figure: data, median model, 68% band, bulge/disk/halo.
3. **`residuals.png`** — data − model; look for trends with \(R\).
4. **`enclosed_mass.png`** — scaled \(M(<R)=v_c^2 R\) with posterior band.
5. **`mcmc_traces.png`** — mixing and burn-in.
6. **`posterior_corner.png`** — pairwise degeneracies.
7. **`posterior_samples.npy` / `.csv`** — burned-in draws for downstream analysis.

### Loading samples

```python
import numpy as np

samples = np.load("results/posterior_samples.npy")
# columns follow param_names in posterior_summary.json
# e.g. Mb, Md, Mh [, Ah] [, sigma]
print(samples.shape)
```

Or open `posterior_samples.csv` in any spreadsheet / pandas.

## Mock recovery

```bash
make mock
# or: python3 scripts/mock_recovery.py
```

Builds a synthetic rotation curve from known \((M_b, M_d, M_h, A_h)\), runs MCMC, and checks that each true value falls inside the 16–84% posterior interval. This is the simplest proof that the sampler works.

## Library API

```python
from galaxy_mcmc import GalaxyPotential, load_rotation_curve
from galaxy_mcmc.likelihood import RotationCurvePosterior
from galaxy_mcmc.mcmc import MetropolisHastings, run_ensemble
import numpy as np
from numpy.random import default_rng

R, V = load_rotation_curve("data/RadialVelocities.dat")
pot = GalaxyPotential()
post = RotationCurvePosterior(
    R, V, pot, sigma=2.2, fit_ah=True, fit_sigma=False
)

# Several chains (recommended)
chains = run_ensemble(post, n_chains=4, n_steps=20_000, burn_in=5_000, seed=42)
samples = np.concatenate([c.samples for c in chains], axis=0)

# Or a single adaptive MH chain
mh = MetropolisHastings(
    post.log_posterior,
    n_params=post.n_params,
    proposal_scale=post.proposal_scale(),
    adapt_until=2000,
    rng=default_rng(0),
)
result = mh.run(
    start=post.default_start(),
    n_steps=10_000,
    burn_in=2_000,
    param_names=post.param_names,
)
```

## Your own galaxy

Provide a whitespace-separated file:

```
#radius(kpc) Velocity(km/s)
0.5  180.0
1.0  200.0
...
```

```bash
python3 scripts/run_inference.py \
  --data my_galaxy.dat \
  --outdir results/my_galaxy \
  --fit-ah --fit-sigma
```

If the galaxy has a very different size, retune fixed scales in `GalaxyPotential(...)` (or free \(A_h\)). Only `model.py` needs to change to swap analytic components.

## Tests & CI

```bash
make test
```

Coverage includes data loading, model sanity, prior bounds, short MH runs, approximate recovery of \(M_d\)/\(M_h\), mock recovery with free \(A_h\), enclosed-mass consistency, and Gelman–Rubin on synthetic chains.

GitHub Actions (`.github/workflows/tests.yml`) runs the same suite on every push and pull request.

## Troubleshooting

| Symptom | Likely fix |
| --- | --- |
| Acceptance ≪ 0.1 | Longer adaptation / smaller steps; try `--quick` first to sanity-check |
| Acceptance ≫ 0.5 | Proposals too timid; increase steps or check starting point |
| \(\hat{R} \gg 1\) | More steps / burn-in; bulge \(M_b\) is often poorly constrained — that is expected |
| Poor fit at large \(R\) | Enable `--fit-ah` |
| Crashes on missing deps | `pip install -r requirements.txt && pip install -e .` |
