import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("fluid_slice.csv")

x = data["x"].values
z = data["z"].values
rho = data["rho"].values
p = data["p"].values
vx = data["vx"].values
vy = data["vy"].values
vz = data["vz"].values

NX = len(np.unique(x))
NZ = len(np.unique(z))

fig, axes = plt.subplots(1, 5, figsize=(28, 6))

rho_grid = rho.reshape(NX, NZ)
im1 = axes[0].imshow(rho_grid, extent=[x.min(), x.max(), z.min(), z.max()], origin="lower", cmap="inferno")
axes[0].set_title("Densité du fluide (ρ)")
fig.colorbar(im1, ax=axes[0])

p_grid = p.reshape(NX, NZ)
im2 = axes[1].imshow(p_grid, extent=[x.min(), x.max(), z.min(), z.max()], origin="lower", cmap="viridis")
axes[1].set_title("Pression du fluide (p)")
fig.colorbar(im2, ax=axes[1])

vx_grid = vx.reshape(NX, NZ)
im3 = axes[2].imshow(vx_grid, extent=[x.min(), x.max(), z.min(), z.max()], origin="lower", cmap="hot")
axes[2].set_title("Composante $v_x$ du fluide")
fig.colorbar(im3, ax=axes[2])

vy_grid = vy.reshape(NX, NZ)
im4 = axes[3].imshow(vy_grid, extent=[x.min(), x.max(), z.min(), z.max()], origin="lower", cmap="hot")
axes[3].set_title("Composante $v_y$ du fluide")
fig.colorbar(im4, ax=axes[3])

vz_grid = vz.reshape(NX, NZ)
im5 = axes[4].imshow(vz_grid, extent=[x.min(), x.max(), z.min(), z.max()], origin="lower", cmap="viridis")
axes[4].set_title("Composante $v_z$ du fluide")
fig.colorbar(im5, ax=axes[4])

plt.suptitle("Distribution et vitesse du fluide en (x, z)", fontsize=16)
plt.tight_layout()
plt.show()
