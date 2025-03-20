
from paraview.simple import *

reader = LegacyVTKReader(FileNames=["geodesic.vtk"])
Show(reader)
rv = GetActiveViewOrCreate("RenderView")
ResetCamera()

rep = GetDisplayProperties(reader, view=rv)
rep.Representation = "Wireframe"
ColorBy(rep, ('POINTS', 'lambda'))
rep.SetScalarBarVisibility(rv, True)

Render()
Interact()
