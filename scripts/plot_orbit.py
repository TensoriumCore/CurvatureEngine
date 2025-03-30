import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

data = np.loadtxt("BH_positions.csv", delimiter=",")
time = data[:, 0]
x1, y1 = data[:, 1], data[:, 2]
x2, y2 = data[:, 4], data[:, 5]

fig, ax = plt.subplots()
line1, = ax.plot([], [], 'ro-', label="BH 1")
line2, = ax.plot([], [], 'bo-', label="BH 2")
ax.set_xlim(np.min([x1, x2]) - 1, np.max([x1, x2]) + 1)
ax.set_ylim(np.min([y1, y2]) - 1, np.max([y1, y2]) + 1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Animation des orbites des trous noirs")
ax.grid(True)
ax.legend()
ax.set_aspect('equal')

def update(frame):
    line1.set_data(x1[:frame], y1[:frame])
    line2.set_data(x2[:frame], y2[:frame])
    return line1, line2

ani = animation.FuncAnimation(fig, update, frames=len(time), interval=100, blit=True)
plt.show()
