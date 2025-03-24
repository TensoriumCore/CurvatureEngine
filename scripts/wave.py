import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import re

files = sorted(glob.glob("tilde_gamma_xx_t=*.csv"),
               key=lambda f: float(re.findall(r"t=([\d.]+)", f)[0].rstrip('.')))

if not files:
    print("❌ Aucun fichier trouvé.")
    exit(1)

data_list = []
time_list = []

for file in files:
    data = np.loadtxt(file, delimiter=",")
    data_list.append(data)
    t = float(re.findall(r"t=([\d.]+)", file)[0].rstrip('.'))
    time_list.append(t)

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(data_list[0][:, 0].min(), data_list[0][:, 0].max())
ax.set_ylim(np.min(data_list) * 1.1, np.max(data_list) * 1.1)
ax.set_xlabel("x")
ax.set_ylabel("tilde_gamma_xx")
title = ax.set_title("")

def init():
    line.set_data([], [])
    return line,

def animate(i):
    x = data_list[i][:, 0]
    y = data_list[i][:, 1]
    line.set_data(x, y)
    title.set_text(f"t = {time_list[i]:.3f}")
    return line, title

ani = animation.FuncAnimation(fig, animate, init_func=init,
                              frames=len(data_list), interval=80, blit=False)

plt.tight_layout()
plt.show()
