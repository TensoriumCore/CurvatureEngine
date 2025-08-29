#!/usr/bin/env python3
import re, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

files = sorted(glob.glob("gauge_xy_t*.csv"),
               key=lambda s: float(re.search(r"t([0-9]+\.[0-9]+)", s).group(1)) if re.search(r"t([0-9]+\.[0-9]+)", s) else np.inf)
if not files:
    raise SystemExit("Aucun fichier gauge_xy_t*.csv trouvé.")

times, frames = [], []
for fp in files:
    t = float(re.search(r"t([0-9]+\.[0-9]+)", fp).group(1))
    df = pd.read_csv(fp)
    xs = np.sort(df['x'].unique()); ys = np.sort(df['y'].unique())
    A  = df.pivot(index='y', columns='x', values='alpha').values
    B0 = df.pivot(index='y', columns='x', values='beta0').values
    B1 = df.pivot(index='y', columns='x', values='beta1').values
    X, Y = np.meshgrid(xs, ys)
    frames.append({'X':X,'Y':Y,'A':A,'B0':B0,'B1':B1})
    times.append(t)
times = np.array(times)

# Échelles globales
A_all = np.stack([fr['A'] for fr in frames], axis=0)
A_vmin, A_vmax = np.nanmin(A_all), np.nanmax(A_all)

# ------ Contrôle de la longueur des flèches ------
# Cap sur la norme de β (percentile pour ignorer les outliers)
Bmag_all = np.stack([np.hypot(fr['B0'], fr['B1']) for fr in frames], axis=0)
Bcap = np.percentile(Bmag_all, 95)  # 95e percentile
eps = 1e-12

def clip_beta(B0, B1):
    mag = np.hypot(B0, B1)
    factor = np.minimum(1.0, Bcap / (mag + eps))
    return B0 * factor, B1 * factor

# Sous-échantillonnage quiver
NY, NX = frames[0]['A'].shape
step_x = max(1, NX // 25)
step_y = max(1, NY // 25)

# Figure
fig, ax = plt.subplots(figsize=(7, 6))
xmin, xmax = frames[0]['X'][0, 0], frames[0]['X'][0, -1]
ymin, ymax = frames[0]['Y'][0, 0], frames[0]['Y'][-1, 0]
extent = (xmin, xmax, ymin, ymax)

im = ax.imshow(frames[0]['A'], origin='lower', extent=extent,
               vmin=A_vmin, vmax=A_vmax, interpolation='nearest', aspect='equal')
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label(r'alpha (α)')

# Quiver initial avec clipping + scale fixe en unités XY
B0c, B1c = clip_beta(frames[0]['B0'], frames[0]['B1'])
Q = ax.quiver(frames[0]['X'][::step_y, ::step_x],
              frames[0]['Y'][::step_y, ::step_x],
              B0c[::step_y, ::step_x], B1c[::step_y, ::step_x],
              pivot='mid', scale_units='xy', scale=Bcap if Bcap > 0 else 1.0)

# Clé de quiver (montre la longueur correspondant à Bcap/2)
if Bcap > 0:
    ax.quiverkey(Q, X=0.88, Y=1.02, U=Bcap/2, label=f"|β| = {Bcap/2:.3g}", labelpos='E')

title = ax.set_title(f"t = {times[0]:.3f}")
ax.set_xlabel("x"); ax.set_ylabel("y")

paused = {'state': False}
def on_click(event): paused['state'] = not paused['state']
def on_key(event): 
    if event.key == ' ': paused['state'] = not paused['state']
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('key_press_event', on_key)

def update(i):
    if paused['state']:
        return im, Q, title
    fr = frames[i]
    im.set_data(fr['A'])
    B0c, B1c = clip_beta(fr['B0'], fr['B1'])
    Q.set_UVC(B0c[::step_y, ::step_x], B1c[::step_y, ::step_x])
    title.set_text(f"t = {times[i]:.3f}")
    return im, Q, title

anim = FuncAnimation(fig, update, frames=len(frames), interval=150, blit=False, repeat=True)
plt.tight_layout()
plt.show()
