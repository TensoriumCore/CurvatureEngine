
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("constraints_evolution.csv")
time = np.abs(df["time"])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Hamiltonian constraint
ax1.plot(time, df["hamiltonian_L2"], color='royalblue', linewidth=0.8)
ax1.set_xlabel("dt", fontsize=12)
ax1.set_ylabel(r"$\|\mathcal{H}\|_{L^2}$", fontsize=12)
ax1.set_title("Hamiltonian constraint (128$^3$ grid)", fontsize=14)
ax1.set_yscale("log")
ax1.grid(True, linestyle="--", alpha=0.6)

# Momentum constraint
momentum_norm = np.linalg.norm(df[["momentum_x_L2", "momentum_y_L2", "momentum_z_L2"]], axis=1)
ax2.plot(time, momentum_norm, color='crimson', linewidth=0.8)
ax2.set_xlabel("dt", fontsize=12)
ax2.set_ylabel(r"$\|\vec{\mathcal{M}}\|_{L^2}$", fontsize=12)
ax2.set_title("Momentum (128$^3$ grid)", fontsize=14)
ax2.set_yscale("log")
ax2.grid(True, linestyle="--", alpha=0.6)

# Supprimer ticks et labels des axes Y
ax1.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
ax2.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

plt.tight_layout()
plt.show()

