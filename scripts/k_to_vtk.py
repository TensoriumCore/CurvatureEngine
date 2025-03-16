from paraview.simple import *

reader = LegacyVTKReader(FileNames=["K_full.vtk"])
Show(reader)
rv = GetActiveViewOrCreate("RenderView")
ResetCamera()

calcK8 = Calculator(Input=reader)
calcK8.ResultArrayName = "K_8_diagonal"
calcK8.Function = "K_8"

contK8 = Contour(Input=calcK8)
contK8.ContourBy = ['POINTS', 'K_8_diagonal']
contK8.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]

Show(contK8, rv)
contK8Rep = GetDisplayProperties(contK8, view=rv)
contK8Rep.Representation = "Surface" 
ColorBy(contK8Rep, ('POINTS', 'K_8_diagonal'))
contK8Rep.SetScalarBarVisibility(rv, True)
contK8Rep.Opacity = 0.1

k8_ctf = GetColorTransferFunction('K_8_diagonal')
k8_ctf.ApplyPreset('Hot Desaturated', True) 
k8_ctf.RescaleTransferFunction(-0.2, 0.2)

contHorizon = Contour(Input=reader)
contHorizon.ContourBy = ['POINTS', 'Horizon']
contHorizon.Isosurfaces = [1.0] 

Show(contHorizon, rv)
contHorizonRep = GetDisplayProperties(contHorizon, view=rv)
contHorizonRep.Representation = "Surface"
contHorizonRep.Opacity = 0.5  
contHorizonRep.DiffuseColor = [0.0, 1.0, 0.0]  
ColorBy(contHorizonRep, ('POINTS', 'Horizon'))
contHorizonRep.SetScalarBarVisibility(rv, False)

contErgosphereExt = Contour(Input=reader)
contErgosphereExt.ContourBy = ['POINTS', 'Ergosphere_plus']
contErgosphereExt.Isosurfaces = [1.0] 

Show(contErgosphereExt, rv)
contErgosphereExtRep = GetDisplayProperties(contErgosphereExt, view=rv)
contErgosphereExtRep.Representation = "Wireframe"
contErgosphereExtRep.Opacity = 0.3  
contErgosphereExtRep.DiffuseColor = [1.0, 0.0, 1.0]
ColorBy(contErgosphereExtRep, ('POINTS', 'Ergosphere_plus'))
contErgosphereExtRep.SetScalarBarVisibility(rv, False)

contErgosphereInt = Contour(Input=reader)
contErgosphereInt.ContourBy = ['POINTS', 'Ergosphere_minus']
contErgosphereInt.Isosurfaces = [1.0] 

Show(contErgosphereInt, rv)
contErgosphereIntRep = GetDisplayProperties(contErgosphereInt, view=rv)
contErgosphereIntRep.Representation = "Surface"
contErgosphereIntRep.Opacity = 0.8  
contErgosphereIntRep.DiffuseColor = [1.0, 1.0, 0.0]
ColorBy(contErgosphereIntRep, ('POINTS', 'Ergosphere_minus'))
contErgosphereIntRep.SetScalarBarVisibility(rv, False)


ResetCamera()
Render()
Interact()
