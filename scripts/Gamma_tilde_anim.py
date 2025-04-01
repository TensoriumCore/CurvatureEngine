import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import os

# Configuration
output_dir = "."
pattern = "gamma_slice_t*.csv"
components = ['00', '01', '02', '10', '11', '12', '20', '21', '22']
n_components = len(components)

def load_all_slices():
    files = sorted(glob.glob(os.path.join(output_dir, pattern)))
    all_data = []
    times = []
    
    for f in files:
        t = float(f.split('_t')[1].replace('.csv', ''))
        df = pd.read_csv(f)
        all_data.append(df)
        times.append(t)
    
    return all_data, np.array(times)

def setup_plot_grid():
    fig, axs = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle('Composantes de la métrique γ_ij', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig, axs

# Fonction de visualisation statique
def plot_single_slice(df, time):
    """Version robuste de la visualisation"""
    try:
        fig, axs = plt.subplots(3, 3, figsize=(15, 12))
        fig.suptitle(f'Composantes de γ_ij à t = {time:.3f}', fontsize=16)
        
        x = df['x'].unique()
        z = df['z'].unique()
        X, Z = np.meshgrid(x, z)
        
        for idx, comp in enumerate(components):
            ax = axs[idx//3, idx%3]
            gamma = df[f'gamma_{comp}'].values.reshape(len(z), len(x))
            
            im = ax.pcolormesh(X, Z, gamma, shading='auto', cmap='viridis')
            ax.set_title(f'γ_{comp}', fontsize=10)
            ax.set_aspect('equal')
            
            # Barre de couleur alternative (sans make_axes_locatable)
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        
        plt.tight_layout()
        output_file = f'gamma_slice_t{time:.3f}.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Visualisation sauvegardée dans {output_file}")
        
    except Exception as e:
        print(f"Erreur lors de la visualisation : {str(e)}")
        if 'fig' in locals():
            plt.close(fig)

def animate_slices():
    all_data, times = load_all_slices()
    fig, axs = setup_plot_grid()
    
    # Initialisation
    x = all_data[0]['x'].unique()
    z = all_data[0]['z'].unique()
    X, Z = np.meshgrid(x, z)
    images = []
    
    for idx, comp in enumerate(components):
        ax = axs[idx//3, idx%3]
        gamma = all_data[0][f'gamma_{comp}'].values.reshape(len(z), len(x))
        
        im = ax.pcolormesh(X, Z, gamma, shading='auto', cmap='viridis', vmin=-2, vmax=2)
        ax.set_title(f'γ_{comp}')
        ax.set_aspect('equal')
        images.append(im)
    
    # Fonction d'animation
    def update(frame):
        for idx, comp in enumerate(components):
            gamma = all_data[frame][f'gamma_{comp}'].values.reshape(len(z), len(x))
            images[idx].set_array(gamma.ravel())
        fig.suptitle(f't = {times[frame]:.2f}', fontsize=16)
        return images
    
    ani = animation.FuncAnimation(fig, update, frames=len(all_data), interval=100, blit=False)
    ani.save('gamma_evolution.mp4', writer='ffmpeg', fps=5, dpi=150)
    plt.show()
    plt.close()
    print("Animation saved as gamma_evolution.mp4")

def plot_metric_norm():
    all_data, times = load_all_slices()
    x = all_data[0]['x'].unique()
    z = all_data[0]['z'].unique()
    X, Z = np.meshgrid(x, z)
    
    for time_idx, df in enumerate(all_data):
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        
        # Calcul de la norme de Frobenius
        norm = np.zeros(len(df))
        for comp in components:
            norm += df[f'gamma_{comp}'].values**2
        norm = np.sqrt(norm).reshape(len(z), len(x))
        
        im = ax.pcolormesh(X, Z, norm, shading='auto', cmap='viridis', norm='log')
        ax.set_title(f'Norme de γ_ij à t = {times[time_idx]:.2f}')
        ax.set_aspect('equal')
        fig.colorbar(im, ax=ax, label='||γ|| (log scale)')
        plt.savefig(f'gamma_norm_t{times[time_idx]:.3f}.png', dpi=150, bbox_inches='tight')
        plt.close()

# Menu principal
def main():
    print("\nVisualisation des slices de γ_ij")
    print("1. Afficher une slice spécifique")
    print("2. Générer une animation temporelle")
    print("3. Visualiser la norme des métriques")
    
    choice = input("Choix (1-3): ")
    
    if choice == '1':
        all_data, times = load_all_slices()
        time_idx = int(input(f"Index du temps (0-{len(times)-1}): "))
        plot_single_slice(all_data[time_idx], times[time_idx])
    elif choice == '2':
        animate_slices()
    elif choice == '3':
        plot_metric_norm()
    else:
        print("Choix invalide")

if __name__ == "__main__":
    main()
