
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from glob import glob

data_folder = "."
pattern = os.path.join(data_folder, "gamma_slice_t*.csv")

files = glob(pattern)
files = [f for f in files if re.search(r"t([\d.]+)\.csv", f)]
files = sorted(files, key=lambda f: float(re.search(r"t([\d.]+)\.csv", f).group(1)))

if not files:
    print("[Erreur] Aucun fichier gamma_slice_t*.csv trouvé.")
    exit(1)

print(f"[INFO] {len(files)} fichiers gamma chargés.")

data0 = np.loadtxt(files[0], delimiter=",", skiprows=1)
N = int(np.sqrt(data0.shape[0]))
X = data0[:, 0].reshape(N, N)
Z = data0[:, 1].reshape(N, N)

comp_names = [
    "gamma_00", "gamma_01", "gamma_02",
    "gamma_10", "gamma_11", "gamma_12",
    "gamma_20", "gamma_21", "gamma_22"
]
component_idx = 3 

fig, ax = plt.subplots()
gamma_data = data0[:, 2 + component_idx].reshape(N, N)
plot = ax.pcolormesh(X, Z, gamma_data, shading='auto', cmap='plasma')
cbar = fig.colorbar(plot, ax=ax)
cbar.set_label(comp_names[component_idx])

ax.set_xlabel("x")
ax.set_ylabel("z")
ax.set_title(f"{comp_names[component_idx]} at t = 0.000")

def update(frame_idx):
    data = np.loadtxt(files[frame_idx], delimiter=",", skiprows=1)
    gamma = data[:, 2 + component_idx].reshape(N, N)
    plot.set_array(gamma.ravel())
    t = re.search(r"t([\d.]+)\.csv", files[frame_idx]).group(1)
    ax.set_title(f"{comp_names[component_idx]} at t = {t}")
    return plot,

ani = animation.FuncAnimation(fig, update, frames=len(files), interval=100, blit=False)
plt.tight_layout()
plt.show()
