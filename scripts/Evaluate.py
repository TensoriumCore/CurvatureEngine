#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import matplotlib.pyplot as plt
import numpy as np

# Nom du fichier texte généré par votre programme C++
filename = "constraints_output.txt"

# -----------------------------------------------------------------------------
# 1) Lecture et extraction des données
# -----------------------------------------------------------------------------

# Regex pour extraire i, j, k
cell_pattern = re.compile(r"Pour la cellule \((\d+),\s*(\d+),\s*(\d+)\):")
# Regex pour extraire Hamiltonian
ham_pattern  = re.compile(r"Hamiltonian = ([\d\.\-eE]+)")
# Regex pour extraire Momentum
mom_pattern  = re.compile(r"Momentum = \[\s*([\d\.\-eE]+),\s*([\d\.\-eE]+),\s*([\d\.\-eE]+)\s*\]")

# Listes pour stocker les données
cells = []          # (i, j, k)
hamiltonians = []   # Valeur de l'Hamiltonian
momentum_x = []
momentum_y = []
momentum_z = []

with open(filename, "r") as f:
    lines = f.readlines()

i_line = 0
while i_line < len(lines):
    line = lines[i_line].strip()
    cell_match = cell_pattern.match(line)
    if cell_match:
        i_val = int(cell_match.group(1))
        j_val = int(cell_match.group(2))
        k_val = int(cell_match.group(3))
        ham_val = None
        mom_x_val = mom_y_val = mom_z_val = None

        if i_line + 1 < len(lines):
            ham_match = ham_pattern.search(lines[i_line + 1])
            if ham_match:
                ham_val = float(ham_match.group(1))

        if i_line + 2 < len(lines):
            mom_match = mom_pattern.search(lines[i_line + 2])
            if mom_match:
                mom_x_val = float(mom_match.group(1))
                mom_y_val = float(mom_match.group(2))
                mom_z_val = float(mom_match.group(3))

        cells.append((i_val, j_val, k_val))
        hamiltonians.append(ham_val)
        momentum_x.append(mom_x_val)
        momentum_y.append(mom_y_val)
        momentum_z.append(mom_z_val)

        i_line += 3  # On saute les trois lignes (cellule, Hamiltonian, Momentum)
    else:
        i_line += 1

cells = np.array(cells)
hamiltonians = np.array(hamiltonians)
momentum_x = np.array(momentum_x)
momentum_y = np.array(momentum_y)
momentum_z = np.array(momentum_z)

# -----------------------------------------------------------------------------
# 2) Déterminer les dimensions (NX, NY, NZ) si nécessaire
#    On suppose que vos indices i, j, k vont de 0 à NX-1, NY-1, NZ-1
# -----------------------------------------------------------------------------
max_i = cells[:,0].max()
max_j = cells[:,1].max()
max_k = cells[:,2].max()

NX = int(max_i + 1)
NY = int(max_j + 1)
NZ = int(max_k + 1)

print(f"Dimensions déduites : NX={NX}, NY={NY}, NZ={NZ}")

# -----------------------------------------------------------------------------
# 3) Visualisation d'un plan 2D (i,j) pour un k fixé (heatmap)
# -----------------------------------------------------------------------------
# Choisissez le k que vous voulez visualiser
k_slice = NZ // 2  # Par défaut, on prend le milieu
# Vous pouvez changer k_slice en fonction de ce que vous voulez inspecter.

# Construire une matrice 2D (NX x NY) pour l'Hamiltonian
H_slice = np.full((NX, NY), np.nan)  # initialisé à NaN

# Idem pour le momentum, si souhaité
Mx_slice = np.full((NX, NY), np.nan)
My_slice = np.full((NX, NY), np.nan)
Mz_slice = np.full((NX, NY), np.nan)

for idx in range(len(cells)):
    i, j, k = cells[idx]
    if k == k_slice:
        H_slice[i, j] = hamiltonians[idx]
        Mx_slice[i, j] = momentum_x[idx]
        My_slice[i, j] = momentum_y[idx]
        Mz_slice[i, j] = momentum_z[idx]

# --- Tracé de l'Hamiltonian en 2D ---
plt.figure(figsize=(7, 5))
plt.title(f"Hamiltonian au plan k={k_slice}")
# Attention : imshow considère l'axe 0 comme vertical (lignes) et l'axe 1 comme horizontal (colonnes).
# On transpose H_slice pour avoir i sur l'axe X et j sur l'axe Y, ou on peut jouer avec origin/extent.
plt.imshow(H_slice.T, origin='lower', 
           extent=[0, NX, 0, NY], aspect='auto', 
           cmap='viridis')
plt.colorbar(label="Hamiltonian")
plt.xlabel("i")
plt.ylabel("j")
plt.tight_layout()
plt.show()

# --- Tracé d'une composante du Momentum en 2D (par ex. Mx) ---
plt.figure(figsize=(7, 5))
plt.title(f"Momentum X au plan k={k_slice}")
plt.imshow(Mx_slice.T, origin='lower', 
           extent=[0, NX, 0, NY], aspect='auto', 
           cmap='coolwarm')
plt.colorbar(label="Momentum X")
plt.xlabel("i")
plt.ylabel("j")
plt.tight_layout()
plt.show()

# -----------------------------------------------------------------------------
# 4) Visualisation 1D : on fixe j et k, et on fait varier i
# -----------------------------------------------------------------------------
j_fixed = NY // 2
k_fixed = NZ // 2

# Récupérer l'Hamiltonian en fonction de i
i_vals = []
H_vals_line = []
for idx in range(len(cells)):
    i, j, k = cells[idx]
    if j == j_fixed and k == k_fixed:
        i_vals.append(i)
        H_vals_line.append(hamiltonians[idx])

# Tri par i croissant
i_vals = np.array(i_vals)
H_vals_line = np.array(H_vals_line)
sort_idx = np.argsort(i_vals)
i_vals = i_vals[sort_idx]
H_vals_line = H_vals_line[sort_idx]


# Définir la taille de la fenêtre de lissage (ajustable)
window_size = 5  # Nombre de points utilisés pour la moyenne mobile

# Calcul de la moyenne mobile (avec gestion des bords)
H_vals_smooth = np.convolve(H_vals_line, np.ones(window_size)/window_size, mode='valid')

# Ajuster les indices pour correspondre à la taille réduite après convolution
i_vals_smooth = i_vals[:len(H_vals_smooth)]

# Tracé avec la courbe moyenne
plt.figure(figsize=(7, 4))
plt.plot(i_vals, H_vals_line, 'o-', label=f'Hamiltonian (j={j_fixed}, k={k_fixed})', alpha=0.6)
plt.plot(i_vals_smooth, H_vals_smooth, '-', linewidth=2, color='red', label="Moyenne mobile")
plt.xlabel("i")
plt.ylabel("Hamiltonian")
plt.title(f"Coupe 1D de l'Hamiltonian avec moyenne mobile (j={j_fixed}, k={k_fixed})")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

