
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_path = "gauge_slice.csv"
df = pd.read_csv(file_path)

x = df["x"].values
z = df["z"].values

gauge_labels = [
    "alpha", "beta0", "beta1", "beta2",
    "d_alpha_dt", "d_beta0_dt", "d_beta1_dt", "d_beta2_dt"
]

NX = len(np.unique(x))  
NZ = len(np.unique(z))  

fig, axes = plt.subplots(2, 4, figsize=(18, 8))  
fig.subplots_adjust(hspace=0.3, wspace=0.3)

for idx, col in enumerate(gauge_labels):
    data = df[col].values.reshape(NX, NZ)
    ax = axes[idx // 4, idx % 4]
    im = ax.imshow(data, extent=[x.min(), x.max(), z.min(), z.max()], origin="lower", aspect="auto", cmap="coolwarm")
    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label=col)

plt.suptitle("Champs de jauge et dérivées temporelles", fontsize=16)
plt.show()
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
fig2.subplots_adjust(wspace=0.3)

for ax, col in zip(axes2, ["d_beta1_dt", "d_beta2_dt"]):
    data = df[col].values.reshape(NX, NZ)
    im = ax.imshow(
        data,
        extent=[x.min(), x.max(), z.min(), z.max()],
        origin="lower",
        aspect="auto",
        cmap="hot",
    )
    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig2.colorbar(im, ax=ax, label=col)

plt.suptitle("d_beta1_dt et d_beta2_dt", fontsize=14)
plt.show()
