# Usage guide

## Installation

```bash
python3 -m pip install -r requirements.txt
python3 -m pip install -e .
```

Requires Python ≥ 3.10, NumPy, Matplotlib, SciPy.

## Command-line interface

```bash
python3 scripts/run_inference.py --help
```

Useful flags:

| Flag | Meaning |
| --- | --- |
| `--data PATH` | Custom `R V` file |
| `--outdir DIR` | Where to write figures / JSON (default `results/`) |
| `--n-steps N` | Steps per chain (default 25000) |
| `--burn-in N` | Discarded warm-up steps (default 5000) |
| `--n-chains N` | Independent MH chains (default 4) |
| `--fit-sigma` | Infer noise \(\sigma\) jointly with masses |
| `--sigma S` | Fixed noise in km/s if not fitting \(\sigma\) |
| `--quick` | Short smoke run |
| `--seed N` | Reproducibility |

### Interpreting outputs

1. Open `results/posterior_summary.txt`.
   - Acceptance rates near \(0.15\)–\(0.4\) are healthy for this proposal.
   - \(\hat{R}\) close to \(1\) (≲ 1.05) indicates chain agreement.
2. Check `mcmc_traces.png` for burn-in and mixing.
3. Use `posterior_corner.png` to spot degeneracies.
4. `rotation_curve_fit.png` is the science figure: data, median model, 68% band, and component contributions.

### Loading posterior samples

```python
import numpy as np
samples = np.load("results/posterior_samples.npy")
# columns: Mb, Md, Mh [, sigma]
Mb, Md, Mh = samples[:, 0], samples[:, 1], samples[:, 2]
```

## Library API (minimal)

```python
from galaxy_mcmc import GalaxyPotential, load_rotation_curve
from galaxy_mcmc.likelihood import RotationCurvePosterior
from galaxy_mcmc.mcmc import MetropolisHastings, run_ensemble

R, V = load_rotation_curve("data/RadialVelocities.dat")
pot = GalaxyPotential()
post = RotationCurvePosterior(R, V, pot, sigma=2.2, fit_sigma=False)

# One chain
from numpy.random import default_rng
import numpy as np
mh = MetropolisHastings(post.log_posterior, n_params=3,
                        proposal_scale=[0.8, 0.04, 0.04],
                        adapt_until=2000, rng=default_rng(0))
result = mh.run(start=np.log([1.0, 1.4e4, 2.6e4]),
                n_steps=10000, burn_in=2000,
                param_names=post.param_names)

# Or several chains + R-hat
chains = run_ensemble(post, n_chains=4, n_steps=20000, burn_in=5000)
```

## Bringing your own galaxy

Provide a whitespace-separated file:

```
#radius(kpc) Velocity(km/s)
0.5  180.0
1.0  200.0
...
```

Then:

```bash
python3 scripts/run_inference.py --data my_galaxy.dat --outdir results/my_galaxy --fit-sigma
```

If your galaxy has a different characteristic size, you may need to retune the fixed scale lengths in `GalaxyPotential(...)` or expose them as free parameters — the code is structured so that only `model.py` must change.

## Tests

```bash
make test
```

Tests cover data loading, model sanity, prior bounds, a short MH run, approximate recovery of \(M_d\) and \(M_h\), and Gelman–Rubin on synthetic chains.
