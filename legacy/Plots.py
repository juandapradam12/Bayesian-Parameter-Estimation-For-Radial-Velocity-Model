"""Plot observed rotation curve vs. legacy C-code fit."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
datos = np.loadtxt(here / "Datos.dat")
ajuste = np.loadtxt(here / "Ajuste.dat")

# Sort for a clean curve
order = np.argsort(datos[:, 0])
r = datos[order, 0]
v = datos[order, 1]
model = ajuste[order]

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.scatter(r, v, s=12, c="k", alpha=0.7, label="Data")
ax.plot(r, model, color="#4C78A8", lw=2, label="MH posterior mean")
ax.set_xlabel("Radius [kpc]")
ax.set_ylabel("Rotation velocity [km/s]")
ax.set_title("Legacy C Metropolis–Hastings fit")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(here / "DatosyModelo.png", dpi=140)
print(f"Saved {here / 'DatosyModelo.png'}")
