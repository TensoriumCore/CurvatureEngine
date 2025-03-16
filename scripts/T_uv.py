import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Chargement du CSV
file_path = "T_energy_momentum.csv"
df = pd.read_csv(file_path)

# Extraction des coordonnées
x = df["x"].values
z = df["z"].values

# Définition des composantes du tenseur énergie-impulsion (4x4)
T_labels = [
    "T_00", "T_01", "T_02", "T_03",
    "T_10", "T_11", "T_12", "T_13",
    "T_20", "T_21", "T_22", "T_23",
    "T_30", "T_31", "T_32", "T_33"
]

# Dimensions de la grille
NX = len(np.unique(x))
NZ = len(np.unique(z))

# Tri des données pour assurer un bon ordonnancement de la grille
df_sorted = df.sort_values(by=["z", "x"])

# Extraction et reshape des composantes de T^{\mu\nu}
T_matrices = {col: df_sorted[col].values.reshape(NZ, NX) for col in T_labels}

# Définition des valeurs min/max pour l'affichage (évite d'être trop influencé par les extrêmes)
all_vals = np.concatenate([T_matrices[col].flatten() for col in T_labels])
vmin, vmax = np.percentile(all_vals, [5, 95])  # Échelle adaptative

# Création des sous-plots pour le tenseur 4x4
fig, axes = plt.subplots(4, 4, figsize=(20, 16))
fig.subplots_adjust(hspace=0.3, wspace=0.3)

# Affichage des cartes de chaleur pour chaque composante de T^{\mu\nu}
for idx, col in enumerate(T_labels):
    data = T_matrices[col]
    ax = axes[idx // 4, idx % 4]
    im = ax.imshow(
        data,
        extent=[x.min(), x.max(), z.min(), z.max()],
        origin="lower",
        aspect="auto",
        cmap="inferno",
        vmin=vmin,
        vmax=vmax
    )
    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label=col)

# Titre principal
plt.suptitle("Tenseur énergie-impulsion \( T^{\mu\nu} \)", fontsize=18)
plt.tight_layout()
plt.show()
