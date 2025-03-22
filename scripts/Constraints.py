import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def load_constraints(filename):
    df = pd.read_csv(filename)
    H = df["Hamiltonian"].to_numpy()
    return H

# Résolutions utilisées
resolutions = [32, 64, 128, 256]
dx_values = [18.0 / N for N in resolutions]

# Charger les erreurs
errors = []
for res in resolutions[:-1]:  # On s'arrête avant la plus fine
    H_coarse = load_constraints(f"Output/constraints_N{res}.csv")
    H_fine = load_constraints(f"Output/constraints_N{res*2}.csv")
    error_L2 = np.sqrt(np.mean((H_coarse - H_fine[::2])**2))  # Comparaison à la grille 2x plus fine
    errors.append(error_L2)

# Calculer l'ordre de convergence
orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
print("Ordres de convergence estimés:", orders)

# Tracé de la convergence
plt.loglog(dx_values[:-1], errors, 'o-', label="Erreur $L_2$")
plt.loglog(dx_values[:-1], errors[0] * (np.array(dx_values[:-1])/dx_values[0])**2, '--', label="Ordre 2")
plt.xlabel(r"$\Delta x$")
plt.ylabel(r"Erreur $L_2$")
plt.legend()
plt.grid()
plt.title("Test de convergence des contraintes BSSN")
plt.show()
