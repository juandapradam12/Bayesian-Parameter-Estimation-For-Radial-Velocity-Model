# Legacy C pipeline

The `legacy/` directory keeps the original homework structure after correcting the algorithmic issues described in [theory.md](theory.md):

- `CurvaRotacion.c` — cleaned Metropolis–Hastings fit
- `Plots.py` — data vs model figure
- `Tarea5.mk` — build/run targets
- `Resultados_hw5.tex` / `.pdf` — original write-up (provenance)

## Build & run

```bash
cd legacy
make -f Tarea5.mk
```

This will:

1. Compile `CurvaRotacion.x`
2. Run MH and write `Datos.dat`, `Ajuste.dat`, `chain.dat`
3. Produce `DatosyModelo.png`

Optional PDF refresh (if you have LaTeX):

```bash
pdflatex Resultados_hw5.tex
```

## When to use it

- Teaching the Metropolis–Hastings algorithm in plain C
- Comparing a minimal implementation to the Python package
- Portfolio / historical context for the original assignment

For analysis, diagnostics, free \(A_h\), mock recovery, and reusable inference, use the Python package in `src/galaxy_mcmc` and the root [README](../README.md).
