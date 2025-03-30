import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import re
import os

def extract_time_from_filename(fname):
    """
    Extrait la valeur de temps à partir du nom de fichier.
    On suppose un format: tilde_gamma_xx_t=XXXX.csv
    """
    # Regex pour trouver tout ce qui est entre 't=' et '.csv'
    
    match = re.search(r"_t=([0-9]+[,.]?[0-9]*)\.csv", fname)


    if match:
        return float(match.group(1))
    else:
        return 0.0

def load_all_data():
    """
    Cherche tous les fichiers CSV correspondants,
    les charge et les trie par ordre croissant de temps.
    Retourne une liste de (time, x-array, gxx-array).
    """
    file_list = glob.glob("tilde_gamma_xx_t=*.csv")
    if not file_list:
        print("Aucun fichier CSV trouvé dans ./Output/ !")
        return []

    # Trie les fichiers par le temps extrait du nom
    file_list = sorted(file_list, key=lambda f: extract_time_from_filename(f))

    data_times = []
    for fname in file_list:
        t_val = extract_time_from_filename(fname)
        # Charge le CSV : on ignore les lignes de commentaire (#)
        x_vals, gxx_vals = np.loadtxt(fname, comments='#', delimiter=',', unpack=True)
        data_times.append((t_val, x_vals, gxx_vals))

    return data_times

def main():
    # Charge toutes les données
    data_times = load_all_data()
    if not data_times:
        return  # Rien à animer

    # Prépare la figure Matplotlib
    fig, ax = plt.subplots()
    ax.set_xlabel('x')
    ax.set_ylabel(r'$\tilde{\gamma}_{xx}$')

    # Premier jeu de données (pour initialiser)
    t0, x0, gxx0 = data_times[0]
    line, = ax.plot(x0, gxx0, label=f"t = {t0:.2f}")
    ax.legend()

    # Définition de la fonction d’animation
    def update(frame_index):
        t, x_vals, gxx_vals = data_times[frame_index]
        line.set_xdata(x_vals)
        line.set_ydata(gxx_vals)
        line.set_label(f"t = {t:.2f}")
        ax.legend()
        ax.set_title("Animation de tilde_gamma_xx")
        return line,

    # Adapter les limites du plot (optionnel) 
    x_min = min([np.min(xvals) for (_, xvals, _) in data_times])
    x_max = max([np.max(xvals) for (_, xvals, _) in data_times])
    gxx_min = min([np.min(gxx) for (_, _, gxx) in data_times])
    gxx_max = max([np.max(gxx) for (_, _, gxx) in data_times])
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(gxx_min, gxx_max)

    # Creation de l’animation
    ani = FuncAnimation(fig, update, frames=len(data_times), interval=500, blit=False, repeat=True)
    plt.show()

if __name__ == "__main__":
    main()
