from paraview.simple import *

reader = LegacyVTKReader(FileNames=["tilde_gamma_full.vtk"])
Show(reader)
rv = GetActiveViewOrCreate("RenderView")
ResetCamera()

print("Available arrays:", reader.PointData.keys())

# Afficher une isosurface pour tilde_gamma_00
calcGamma00 = Calculator(Input=reader)
calcGamma00.ResultArrayName = "Gamma_00"
calcGamma00.Function = "tilde_gamma_00"
contGamma00 = Contour(Input=calcGamma00)
contGamma00.ContourBy = ['POINTS', 'Gamma_00']
contGamma00.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]
Show(contGamma00, rv)
contGamma00Rep = GetDisplayProperties(contGamma00, view=rv)
contGamma00Rep.Representation = "Wireframe"
ColorBy(contGamma00Rep, ('POINTS', 'Gamma_00'))
contGamma00Rep.SetScalarBarVisibility(rv, True)
contGamma00Rep.Opacity = 0.1
gamma00_ctf = GetColorTransferFunction('Gamma_00')
gamma00_ctf.ApplyPreset('Cool to Warm', True)
gamma00_ctf.RescaleTransferFunction(-0.2, 0.2)

# Afficher une isosurface pour tilde_gamma_11
calcGamma11 = Calculator(Input=reader)
calcGamma11.ResultArrayName = "Gamma_11"
calcGamma11.Function = "tilde_gamma_11"
contGamma11 = Contour(Input=calcGamma11)
contGamma11.ContourBy = ['POINTS', 'Gamma_11']
contGamma11.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]
Show(contGamma11, rv)
contGamma11Rep = GetDisplayProperties(contGamma11, view=rv)
contGamma11Rep.Representation = "Wireframe"
ColorBy(contGamma11Rep, ('POINTS', 'Gamma_11'))
contGamma11Rep.SetScalarBarVisibility(rv, True)
contGamma11Rep.Opacity = 0.1
gamma11_ctf = GetColorTransferFunction('Gamma_11')
gamma11_ctf.ApplyPreset('Cool to Warm', True)
gamma11_ctf.RescaleTransferFunction(-0.2, 0.2)

# Afficher une isosurface pour tilde_gamma_22
calcGamma22 = Calculator(Input=reader)
calcGamma22.ResultArrayName = "Gamma_22"
calcGamma22.Function = "tilde_gamma_22"
contGamma22 = Contour(Input=calcGamma22)
contGamma22.ContourBy = ['POINTS', 'Gamma_22']
contGamma22.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]
Show(contGamma22, rv)
contGamma22Rep = GetDisplayProperties(contGamma22, view=rv)
contGamma22Rep.Representation = "Wireframe"
ColorBy(contGamma22Rep, ('POINTS', 'Gamma_22'))
contGamma22Rep.SetScalarBarVisibility(rv, True)
contGamma22Rep.Opacity = 0.1
gamma22_ctf = GetColorTransferFunction('Gamma_22')
gamma22_ctf.ApplyPreset('Cool to Warm', True)
gamma22_ctf.RescaleTransferFunction(-0.2, 0.2)

# Finalisation de la visualisation
ResetCamera()
Render()
Interact()
