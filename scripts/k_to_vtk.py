from paraview.simple import *

reader = LegacyVTKReader(FileNames=["K_full.vtk"])
Show(reader)
rv = GetActiveViewOrCreate("RenderView")
ResetCamera()
#
# calcTraceK = Calculator(Input=reader)
# calcTraceK.ResultArrayName = "TraceK"
# calcTraceK.Function = "K_0 + K_4 + K_8"
# Show(calcTraceK, rv)
# traceKRep = GetDisplayProperties(calcTraceK, view=rv)
# ColorBy(traceKRep, ('POINTS', 'TraceK'))
# traceKRep.SetScalarBarVisibility(rv, True)
# traceKRep.Opacity = 0.1
#
calcK8 = Calculator(Input=reader)
calcK8.ResultArrayName = "K_8_diagonal"
calcK8.Function = "K_8"
contK8 = Contour(Input=calcK8)
contK8.ContourBy = ['POINTS', 'K_8_diagonal']
contK8.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]
Show(contK8, rv)
contK8Rep = GetDisplayProperties(contK8, view=rv)
contK8Rep.Representation = "Wireframe"
ColorBy(contK8Rep, ('POINTS', 'K_8_diagonal'))
contK8Rep.SetScalarBarVisibility(rv, True)
contK8Rep.Opacity = 0.1
k8_ctf = GetColorTransferFunction('K_8_diagonal')
k8_ctf.ApplyPreset('Hot Desaturated', True)
k8_ctf.RescaleTransferFunction(-0.2, 0.2)


calc_dKt8 = Calculator(Input=reader)
calc_dKt8.ResultArrayName = "dKt_8_diagonal"
calc_dKt8.Function = "dKt_8"
cont_dKt8 = Contour(Input=calc_dKt8)
cont_dKt8.ContourBy = ['POINTS', 'dKt_8_diagonal']
cont_dKt8.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]
Show(cont_dKt8, rv)
cont_dKt8Rep = GetDisplayProperties(cont_dKt8, view=rv)
cont_dKt8Rep.Representation = "Wireframe"
ColorBy(cont_dKt8Rep, ('POINTS', 'dKt_8_diagonal'))
cont_dKt8Rep.SetScalarBarVisibility(rv, True)
cont_dKt8Rep.Opacity = 0.1
dKt8_ctf = GetColorTransferFunction('dKt_8_diagonal')
dKt8_ctf.ApplyPreset('Cool to Warm', True)
dKt8_ctf.RescaleTransferFunction(-0.2, 0.2)
#
# calc_dkt4 = Calculator(Input=reader)
# calc_dkt4.ResultArrayName = "dKt_4_diagonal"
# calc_dkt4.Function = "dKt_4"
# cont_dkt4 = Contour(Input=calc_dkt4)
# cont_dkt4.ContourBy = ['POINTS', 'dKt_4_diagonal']
# cont_dkt4.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]
# Show(cont_dkt4, rv)
# cont_dkt4Rep = GetDisplayProperties(cont_dkt4, view=rv)
# cont_dkt4Rep.Representation = "Wireframe"
# ColorBy(cont_dkt4Rep, ('POINTS', 'dKt_4_diagonal'))
# cont_dkt4Rep.SetScalarBarVisibility(rv, True)
# cont_dkt4Rep.Opacity = 0.1
# dkt4_ctf = GetColorTransferFunction('dKt_4_diagonal')
# dkt4_ctf.ApplyPreset('Cool to Warm', True)
# dkt4_ctf.RescaleTransferFunction(-0.2, 0.2)
#
# calc_dkt0 = Calculator(Input=reader)
# calc_dkt0.ResultArrayName = "dKt_0_diagonal"
# calc_dkt0.Function = "dKt_0"
# cont_dkt0 = Contour(Input=calc_dkt0)
# cont_dkt0.ContourBy = ['POINTS', 'dKt_0_diagonal']
# cont_dkt0.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]
# Show(cont_dkt0, rv)
# cont_dkt0Rep = GetDisplayProperties(cont_dkt0, view=rv)
# cont_dkt0Rep.Representation = "Wireframe"
# ColorBy(cont_dkt0Rep, ('POINTS', 'dKt_0_diagonal'))
# cont_dkt0Rep.SetScalarBarVisibility(rv, True)
# cont_dkt0Rep.Opacity = 0.1
# dkt0_ctf = GetColorTransferFunction('dKt_0_diagonal')
# dkt0_ctf.ApplyPreset('Cool to Warm', True)
# dkt0_ctf.RescaleTransferFunction(-0.2, 0.2)
#
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
