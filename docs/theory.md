# Theory: rotation curves and Bayesian mass estimation

## From mass to circular velocity

In an axisymmetric galactic potential \(\Phi(R,z)\), stars on circular orbits in the mid-plane obey

\[
\frac{v_c^2(R)}{R} = \left.\frac{\partial\Phi}{\partial R}\right|_{z=0}.
\]

For a multi-component galaxy one decomposes \(\Phi=\Phi_b+\Phi_d+\Phi_h\), so

\[
v_c^2(R)=v_b^2(R)+v_d^2(R)+v_h^2(R).
\]

This project uses analytic profiles common in Milky-Way potential models
(Plummer bulge, Miyamoto–Nagai disk, isothermal / Allen–Santillán-like halo):

| Component | Density / potential family | Role in \(v_c(R)\) |
| --- | --- | --- |
| Bulge | Plummer sphere | Peaked central contribution |
| Disk | Miyamoto–Nagai | Intermediate radii |
| Halo | Extended isothermal-like | Flat / slowly declining outer curve |

With gravitational constant absorbed into the mass units, the explicit formulae are those in the README. Scale lengths \((B_b,B_d,A_d,A_h)\) are held fixed; only the mass normalizations \((M_b,M_d,M_h)\) are free.

### Why not fit scale lengths too?

Geometric scales are strongly degenerate with masses on a single noisy rotation curve. Fixing literature values (Irrgang / Allen–Santillán style constants) yields a well-posed three-parameter inference problem that still answers the scientific question: **how much mass sits in each component?**

## Bayesian formulation

Parameters \(\theta=(M_b,M_d,M_h)\) [and optionally \(\sigma\)] are inferred from data \(D=\{(R_i,V_i)\}\):

\[
p(\theta\mid D)\propto p(D\mid\theta)\,p(\theta).
\]

**Likelihood.** Independent Gaussian measurement errors with scale \(\sigma\):

\[
p(D\mid\theta)=\prod_{i=1}^{N}\mathcal{N}\big(V_i;\,v_c(R_i;\theta),\,\sigma\big).
\]

**Prior.** Uniform on \(\log M_j\) over a wide interval (Jeffreys-like positivity prior). When \(\sigma\) is free, a uniform prior on \(\log\sigma\) is used.

**Posterior exploration.** Random-walk Metropolis–Hastings in \(\log M\) space with Gaussian proposals. Proposal widths are adapted during warm-up toward an acceptance rate near \(0.25\). Several independent chains are compared with the Gelman–Rubin statistic \(\hat{R}\).

## What the posterior tells you

- **Marginals** \(p(M_d\mid D)\), \(p(M_h\mid D)\): mass constraints with credible intervals.
- **Joint structure**: corner plots reveal mass degeneracies (e.g. disk–halo trade-offs).
- **Predictive band**: drawing \(v_c(R)\) from posterior samples visualizes model uncertainty, not just a best-fit curve.
- **Component curves**: even when \(M_b\) is consistent with zero, the disk and halo contributions clarify *where* each mass matters.

## Historical bugs fixed in this rewrite

The original homework C code intended the model above but contained several defects that made results unreliable:

1. Integer division (`3/4 == 0` in C) zeroed key exponents.
2. `pow(R, R)` was used instead of \(R^2\).
3. The likelihood used a single data point instead of the full curve.
4. Both branches of the Metropolis accept/reject step wrote the proposal (no true rejection).
5. Masses were drawn in \([0,1]\), far from the scale needed by the data.
6. Disk scale used \(B_b+A_d\) instead of \(B_d+A_d\).

The Python package and the cleaned `legacy/CurvaRotacion.c` implement the corrected physics and a proper MH update.
