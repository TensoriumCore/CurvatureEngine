from paraview.simple import *

reader = LegacyVTKReader(FileNames=["K_full.vtk"])
Show(reader)
rv = GetActiveViewOrCreate("RenderView")
ResetCamera()

# -------------------------------------------------------------
# K_0 : slice + color map + opacité
# -------------------------------------------------------------
calcK0 = Calculator(Input=reader)
calcK0.ResultArrayName = "K_0_diagonal"
calcK0.Function = "K_0"

sliceK0 = Slice(Input=calcK0)
sliceK0.SliceType = "Plane"
sliceK0.SliceType.Origin = [0, 0, 0]
sliceK0.SliceType.Normal = [0, 1, 0]
Show(sliceK0, rv)
sliceK0Rep = GetDisplayProperties(sliceK0, view=rv)
sliceK0Rep.Representation = "Surface"
ColorBy(sliceK0Rep, ('POINTS', 'K_0_diagonal'))
sliceK0Rep.SetScalarBarVisibility(rv, True)
sliceK0Rep.Opacity = 0.3

k0_ctf = GetColorTransferFunction('K_0_diagonal')
k0_ctf.ApplyPreset('Cool to Warm (Extended)', True)
k0_ctf.RescaleTransferFunction(-0.2, 0.2)

# -------------------------------------------------------------
# K_4 : volume + opacité "peak" + autre palette
# -------------------------------------------------------------
# calcK4 = Calculator(Input=reader)
# calcK4.ResultArrayName = "K_4_diagonal"
# calcK4.Function = "K_4"
# Show(calcK4, rv)
# calcK4Rep = GetDisplayProperties(calcK4, view=rv)
# calcK4Rep.Representation = "Volume"
# ColorBy(calcK4Rep, ('POINTS', 'K_4_diagonal'))
# calcK4Rep.SetScalarBarVisibility(rv, True)
#
# k4_ctf = GetColorTransferFunction('K_4_diagonal')
# k4_ctf.ApplyPreset('RdBu', True)
# k4_ctf.RescaleTransferFunction(-0.2, 0.2)
#
# k4_otf = GetOpacityTransferFunction('K_4_diagonal')
# k4_otf.Points = [
#     -0.2, 0.0, 0.5, 0.0,
#     -0.05, 0.7, 0.5, 0.0,
#     0.0, 1.0, 0.5, 0.0,
#     0.05, 0.7, 0.5, 0.0,
#     0.2, 0.0, 0.5, 0.0
# ]
#
# -------------------------------------------------------------
# K_8 : isosurfaces + color map
# -------------------------------------------------------------
calcK8 = Calculator(Input=reader)
calcK8.ResultArrayName = "K_8_diagonal"
calcK8.Function = "K_8"
contK8 = Contour(Input=calcK8)
contK8.ContourBy = ['POINTS', 'K_8_diagonal']
contK8.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]
Show(contK8, rv)
contK8Rep = GetDisplayProperties(contK8, view=rv)
ColorBy(contK8Rep, ('POINTS', 'K_8_diagonal'))
contK8Rep.SetScalarBarVisibility(rv, True)
contK8Rep.Opacity = 0.3

k8_ctf = GetColorTransferFunction('K_8_diagonal')
k8_ctf.ApplyPreset('Rainbow Desaturated', True)
k8_ctf.RescaleTransferFunction(-0.2, 0.2)
print("Nombre de cellules de l'isosurface : ", contK8.GetDataInformation().GetNumberOfCells())
#
# -------------------------------------------------------------
# Trace : K_0 + K_4 + K_8 en volume
# -------------------------------------------------------------
# calcTrace = Calculator(Input=reader)
# calcTrace.ResultArrayName = "K_trace"
# calcTrace.Function = "K_0 + K_4 + K_8"
# Show(calcTrace, rv)
# traceRep = GetDisplayProperties(calcTrace, view=rv)
# traceRep.Representation = "Volume"
# ColorBy(traceRep, ('POINTS', 'K_trace'))
# traceRep.SetScalarBarVisibility(rv, True)
# traceRep.Opacity = 0.3
#
# trace_ctf = GetColorTransferFunction('K_trace')
# trace_ctf.ApplyPreset('Viridis (matplotlib)', True)
# trace_ctf.RescaleTransferFunction(-0.3, 0.3)
#
# trace_otf = GetOpacityTransferFunction('K_trace')
# trace_otf.Points = [
#     -0.3, 0.0, 0.5, 0.0,
#     -0.1, 0.6, 0.5, 0.0,
#     0.0, 1.0, 0.5, 0.0,
#     0.1, 0.6, 0.5, 0.0,
#     0.3, 0.0, 0.5, 0.0
# ]

import numpy as np
# Paramètres Kerr
M = 1.0
a = 0.9
L = 128.0  # Domaine spatial
theta_res = 100  # Résolution angulaire
phi_res = 100    # Résolution polaire

r_H = M + np.sqrt(M**2 - a**2)  # Horizon des événements
r_e_plus = lambda theta: M + np.sqrt(M**2 - a**2 * np.cos(theta)**2)  # Ergosphère externe
r_e_minus = lambda theta: M - np.sqrt(M**2 - a**2 * np.cos(theta)**2)  # Ergosphère interne

# Fonction pour générer des points de surface paramétrique
def generate_surface_points(radius_func, color, opacity=0.3):
    theta_vals = np.linspace(0, np.pi, theta_res)
    phi_vals = np.linspace(0, 2*np.pi, phi_res)
    points = []

    for theta in theta_vals:
        for phi in phi_vals:
            r = radius_func(theta)
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            points.append([x, y, z])

    # Création de la source programmée pour ParaView
    surface_source = ProgrammableSource()
    surface_source.OutputDataSetType = 'vtkPolyData'
    surface_source.Script = f"""
import vtk
import numpy as np

points = {points}
num_points = len(points)

polydata = self.GetOutput()
pts = vtk.vtkPoints()
cells = vtk.vtkCellArray()

for p in points:
    pts.InsertNextPoint(p)

for i in range(num_points - 1):
    line = vtk.vtkLine()
    line.GetPointIds().SetId(0, i)
    line.GetPointIds().SetId(1, i + 1)
    cells.InsertNextCell(line)

polydata.SetPoints(pts)
polydata.SetLines(cells)
"""

    Show(surface_source, rv)
    surfaceRep = GetDisplayProperties(surface_source, view=rv)
    surfaceRep.Representation = "Surface"
    surfaceRep.AmbientColor = color  # Couleur de la surface
    surfaceRep.Opacity = opacity  # Transparence
    return surface_source

# ----------------------------------------
# Affichage des différentes surfaces de Kerr
# ----------------------------------------

# Horizon des événements (rouge)
horizon_source = generate_surface_points(lambda theta: r_H, [1, 0, 0], 0.5)

# Ergosphère externe (bleu)
ergosphere_plus_source = generate_surface_points(r_e_plus, [0, 0, 1], 0.3)

# Ergosphère interne (vert)
ergosphere_minus_source = generate_surface_points(r_e_minus, [0, 1, 0], 0.3)
ResetCamera()
Render()
ResetCamera()
Render()
Interact()
