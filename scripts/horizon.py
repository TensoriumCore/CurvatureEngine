
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_horizon(csv_file):
    data = np.loadtxt(csv_file, delimiter=",", skiprows=1)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, s=5, c='red', alpha=0.5)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Points proches de l'horizon")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_horizon("horizon.csv")
