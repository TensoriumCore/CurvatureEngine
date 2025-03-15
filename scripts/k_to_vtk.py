from paraview.simple import *

# Charger le fichier VTK
reader = LegacyVTKReader(FileNames=["K_full.vtk"])
Show(reader)
rv = GetActiveViewOrCreate("RenderView")
ResetCamera()

# -------------------------------------------------------------
# K_8 : Isosurfaces + Color map + Transparence
# -------------------------------------------------------------
calcK8 = Calculator(Input=reader)
calcK8.ResultArrayName = "K_8_diagonal"
calcK8.Function = "K_8"

contK8 = Contour(Input=calcK8)
contK8.ContourBy = ['POINTS', 'K_8_diagonal']
contK8.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]

Show(contK8, rv)
contK8Rep = GetDisplayProperties(contK8, view=rv)
contK8Rep.Representation = "Wireframe"  # Mode grille
ColorBy(contK8Rep, ('POINTS', 'K_8_diagonal'))
contK8Rep.SetScalarBarVisibility(rv, True)
contK8Rep.Opacity = 0.3

k8_ctf = GetColorTransferFunction('K_8_diagonal')
k8_ctf.ApplyPreset('Cool to Warm (Extended)', True)  # Palette plus moderne
k8_ctf.RescaleTransferFunction(-0.2, 0.2)

# -------------------------------------------------------------
# Affichage de l'Horizon des Événements
# -------------------------------------------------------------
contHorizon = Contour(Input=reader)
contHorizon.ContourBy = ['POINTS', 'Horizon']
contHorizon.Isosurfaces = [1.0]  # Wireframe de l'horizon

Show(contHorizon, rv)
contHorizonRep = GetDisplayProperties(contHorizon, view=rv)
contHorizonRep.Representation = "Wireframe"
contHorizonRep.Opacity = 0.5  # Transparence partielle
contHorizonRep.DiffuseColor = [0.0, 1.0, 0.0]  # Rouge pour l'horizon
ColorBy(contHorizonRep, ('POINTS', 'Horizon'))
contHorizonRep.SetScalarBarVisibility(rv, False)

# -------------------------------------------------------------
# Affichage de l'Ergosphère Externe
# -------------------------------------------------------------
contErgosphereExt = Contour(Input=reader)
contErgosphereExt.ContourBy = ['POINTS', 'Ergosphere_plus']
contErgosphereExt.Isosurfaces = [1.0]  # Wireframe de l'ergosphère externe

Show(contErgosphereExt, rv)
contErgosphereExtRep = GetDisplayProperties(contErgosphereExt, view=rv)
contErgosphereExtRep.Representation = "Wireframe"
contErgosphereExtRep.Opacity = 0.3  # Transparence élevée
contErgosphereExtRep.DiffuseColor = [0.0, 0.0, 1.0]  # Bleu pour l'ergosphère externe
ColorBy(contErgosphereExtRep, ('POINTS', 'Ergosphere_plus'))
contErgosphereExtRep.SetScalarBarVisibility(rv, False)

# -------------------------------------------------------------
# Affichage de l'Ergosphère Interne
# -------------------------------------------------------------
contErgosphereInt = Contour(Input=reader)
contErgosphereInt.ContourBy = ['POINTS', 'Ergosphere_minus']
contErgosphereInt.Isosurfaces = [1.0]  # Wireframe de l'ergosphère interne

Show(contErgosphereInt, rv)
contErgosphereIntRep = GetDisplayProperties(contErgosphereInt, view=rv)
contErgosphereIntRep.Representation = "Wireframe"
contErgosphereIntRep.Opacity = 0.4  # Transparence intermédiaire
contErgosphereIntRep.DiffuseColor = [0.0, 1.0, 0.0]  # Vert pour l'ergosphère interne
ColorBy(contErgosphereIntRep, ('POINTS', 'Ergosphere_minus'))
contErgosphereIntRep.SetScalarBarVisibility(rv, False)

# -------------------------------------------------------------
# Rendu final et Interactivité
# -------------------------------------------------------------
ResetCamera()
Render()
Interact()
