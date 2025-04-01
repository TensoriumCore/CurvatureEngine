import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from glob import glob

data_folder = "." 
file_pattern = os.path.join(data_folder, "chi_slice_t*.dat")

files = glob(file_pattern)
files = [f for f in files if re.search(r"t([\d\.]+)\.dat$", f)]
files = sorted(files, key=lambda f: float(re.search(r"t([\d\.]+)\.dat$", f).group(1)))

print(f"[INFO] Found {len(files)} chi slices.")

data0 = np.loadtxt(files[0])
N = int(np.sqrt(data0.shape[0]))
X = data0[:, 0].reshape(N, N)
Y = data0[:, 1].reshape(N, N)

fig, ax = plt.subplots()
cmap = plt.get_cmap("plasma")
chi_data = data0[:, 2].reshape(N, N)
plot = ax.pcolormesh(X, Y, chi_data, cmap=cmap, shading="auto")
cbar = fig.colorbar(plot, ax=ax)
cbar.set_label("chi")

ax.set_title("chi evolution")
ax.set_xlabel("x")
ax.set_ylabel("y")

def update(frame_idx):
    data = np.loadtxt(files[frame_idx])
    chi = data[:, 2].reshape(N, N)
    plot.set_array(chi.ravel())
    time = re.search(r"t([\d\.]+)\.dat$", files[frame_idx]).group(1)
    ax.set_title(f"chi at t = {time}")
    return plot,

# === Lancer l’animation ===
ani = animation.FuncAnimation(
    fig, update, frames=len(files), interval=100, blit=False, repeat=True
)

plt.tight_layout()
plt.show()
