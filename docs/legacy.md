# Legacy C pipeline

The `legacy/` directory preserves the original homework structure (`CurvaRotacion.c`, `Plots.py`, `Tarea5.mk`, and the old PDF write-up) after correcting the algorithmic issues documented in [theory.md](theory.md).

## Build & run

```bash
cd legacy
make -f Tarea5.mk
```

This will:

1. Compile `CurvaRotacion.x`
2. Run Metropolis–Hastings, writing `Datos.dat`, `Ajuste.dat`, and `chain.dat`
3. Produce `DatosyModelo.png`

## When to use it

- Teaching the MH algorithm in plain C without Python dependencies
- Comparing a minimal implementation against the Python package
- Historical / portfolio context for the original assignment

For analysis, diagnostics, and reusable inference, prefer the Python package in `src/galaxy_mcmc`.
