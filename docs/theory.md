# Theory: rotation curves and Bayesian mass estimation

This note explains the physics and inference setup behind **Galaxy MCMC**. For commands and outputs, see [usage.md](usage.md).

## From mass to circular velocity

In an axisymmetric galactic potential $\Phi(R,z)$, stars on circular mid-plane orbits obey

$$
\frac{v_c^2(R)}{R} = \left.\frac{\partial\Phi}{\partial R}\right|_{z=0}.
$$

For a multi-component galaxy, $\Phi=\Phi_b+\Phi_d+\Phi_h$, so

$$
v_c^2(R)=v_b^2(R)+v_d^2(R)+v_h^2(R).
$$

This project uses analytic profiles standard in Milky-Way-style potential models:

| Component | Family | Role in $v_c(R)$ |
| --- | --- | --- |
| Bulge | Plummer sphere | Central rise |
| Disk | Miyamoto–Nagai | Intermediate radii |
| Halo | Isothermal / Allen–Santillán-like | Outer, slowly declining curve |

With $G$ absorbed into the mass units, the explicit formulae are:

$$
\begin{aligned}
v_b^2 &= M_b\, R^2\big/\big(R^2+B_b^2\big)^{3/2},\\
v_d^2 &= M_d\, R^2\big/\big(R^2+(B_d+A_d)^2\big)^{3/2},\\
v_h^2 &= M_h\big/\big(R^2+A_h^2\big)^{1/2}.
\end{aligned}
$$

Default geometric scales (kpc): $B_b=0.2497$, $B_d=5.16$, $A_d=0.3105$, $A_h=64.3$.

### Which parameters are free?

| Parameter | Default | Notes |
| --- | --- | --- |
| $M_b, M_d, M_h$ | **always free** | Scaled masses |
| $A_h$ | fixed, optional free (`fit_ah=True`) | Halo scale; often informative on this dataset |
| $B_b, B_d, A_d$ | fixed | Strongly degenerate with masses on a single noisy curve |
| $\sigma$ | fixed or free (`fit_sigma=True`) | Gaussian noise [km/s] |

Keeping disk/bulge scales fixed yields a well-posed problem that still answers the science question: **how much mass sits in each component?** Freeing $A_h$ is a small, useful extension without exploding dimensionality.

### Enclosed mass

From $v_c^2 = G M(< R)/R$ with $G$ absorbed,

$$
M(< R) = v_c^2(R)\, R
$$

in the same scaled mass units. The pipeline plots a posterior band for this profile.

## Bayesian formulation

Parameters $\theta$ (log-space sampling of the free quantities above) are inferred from data $D=\{(R_i,V_i)\}$:

$$
p(\theta\mid D)\propto p(D\mid\theta)\,p(\theta).
$$

**Likelihood** — independent Gaussian errors:

$$
\log\mathcal{L}(\theta)
= -\tfrac12\sum_i\left(\frac{V_i-v_c(R_i;\theta)}{\sigma}\right)^2
- N\log\sigma + \mathrm{const}.
$$

**Priors** — uniform on $\log M_j$ over a wide box; when enabled, uniform on $\log A_h$ and $\log\sigma$ within broad physical bounds.

**Sampler** — random-walk Metropolis–Hastings with Gaussian proposals in log-space. Proposal widths adapt during warm-up toward ~25% acceptance. Several independent chains are compared with Gelman–Rubin $\hat{R}$.

## Reading the posterior

- **Marginals** on $M_d$, $M_h$ (and $A_h$, $\sigma$): mass / scale constraints with 16–84% intervals.
- **Corner plots**: reveal disk–halo (and mass–scale) degeneracies.
- **Predictive band on $v_c(R)$**: model uncertainty, not just a point fit.
- **Component curves**: show *where* bulge, disk, and halo matter.
- **Residuals**: check for systematic mismodeling vs radius.
- **Enclosed mass**: cumulative mass growth with radius.

On the included dataset, $M_b$ is typically consistent with near-zero: the data do not require a significant bulge.

## Historical bugs fixed in this rewrite

The original homework C code aimed at this model but was unreliable:

1. Integer division (`3/4 == 0` in C) zeroed key exponents.
2. `pow(R, R)` was used instead of $R^2$.
3. The likelihood used a single data point instead of the full curve.
4. Both Metropolis branches wrote the proposal (no true rejection).
5. Masses were drawn in $[0,1]$, far from the data scale.
6. Disk scale used $B_b+A_d$ instead of $B_d+A_d$.

The Python package and cleaned `legacy/CurvaRotacion.c` implement the corrected physics and a proper MH update. A mock-recovery test injects known parameters and checks that truth falls in the 68% posterior interval.
