from paraview.simple import *

# Chargement du fichier VTK
reader = LegacyVTKReader(FileNames=["K_full.vtk"])
Show(reader)
rv = GetActiveViewOrCreate("RenderView")
ResetCamera()

contHorizon = Contour(Input=reader)
contHorizon.ContourBy = ['POINTS', 'Horizon']
contHorizon.Isosurfaces = [1.0] 

Show(contHorizon, rv)
contHorizonRep = GetDisplayProperties(contHorizon, view=rv)
contHorizonRep.Representation = "Surface"
contHorizonRep.Opacity = 0.5  
contHorizonRep.DiffuseColor = [0.0, 1.0, 0.0]  # Vert
ColorBy(contHorizonRep, ('POINTS', 'Horizon'))
contHorizonRep.SetScalarBarVisibility(rv, False)

# ========================== Ergosphère ========================== #
contErgosphere = Contour(Input=reader)
contErgosphere.ContourBy = ['POINTS', 'Ergosphere_plus']
contErgosphere.Isosurfaces = [1.0] 

Show(contErgosphere, rv)
contErgosphereRep = GetDisplayProperties(contErgosphere, view=rv)
contErgosphereRep.Representation = "Wireframe"
contErgosphereRep.Opacity = 0.3  
contErgosphereRep.DiffuseColor = [1.0, 0.0, 1.0]  # Magenta
ColorBy(contErgosphereRep, ('POINTS', 'Ergosphere_plus'))
contErgosphereRep.SetScalarBarVisibility(rv, False)

# ========================== Filtrage du disque d'accrétion (fluide) ========================== #
thresholdFluid = Threshold(Input=reader)
thresholdFluid.Scalars = ['POINTS', 'fluid']
thresholdFluid.ThresholdMethod = "Between"  # Sélectionne un intervalle
thresholdFluid.LowerThreshold = 0.05        # Valeur minimale pour le seuil
thresholdFluid.UpperThreshold = 1.0         # Valeur maximale pour le seuil

contFluid = Contour(Input=thresholdFluid)
contFluid.ContourBy = ['POINTS', 'fluid']
contFluid.Isosurfaces = [0.1, 0.5, 1.0]  # Sélection des régions denses

Show(contFluid, rv)
contFluidRep = GetDisplayProperties(contFluid, view=rv)
contFluidRep.Representation = "Surface"
contFluidRep.Opacity = 0.8  
contFluidRep.DiffuseColor = [0.0, 0.0, 1.0]  # Bleu
ColorBy(contFluidRep, ('POINTS', 'fluid'))
contFluidRep.SetScalarBarVisibility(rv, False)

# ========================== Coupe équatoriale du disque ========================== #
sliceDisk = Slice(Input=thresholdFluid)
sliceDisk.SliceType = "Plane"
sliceDisk.SliceType.Origin = [0.0, 0.0, 0.0]  # Centre de la coupe
sliceDisk.SliceType.Normal = [0.0, 0.0, 1.0]  # Coupe perpendiculaire à Z

Show(sliceDisk, rv)
sliceDiskRep = GetDisplayProperties(sliceDisk, view=rv)
sliceDiskRep.Representation = "Surface"
sliceDiskRep.Opacity = 0.8  
sliceDiskRep.DiffuseColor = [0.5, 0.5, 1.0]  # Bleu clair
ColorBy(sliceDiskRep, ('POINTS', 'fluid'))
sliceDiskRep.SetScalarBarVisibility(rv, True)

# ========================== Vitesse du fluide ========================== #
contFluidVelocity = Contour(Input=reader)
contFluidVelocity.ContourBy = ['POINTS', 'fluid_velocity']
contFluidVelocity.Isosurfaces = [0.1, 0.3, 0.5]

Show(contFluidVelocity, rv)
contFluidVelocityRep = GetDisplayProperties(contFluidVelocity, view=rv)
contFluidVelocityRep.Representation = "Wireframe"
contFluidVelocityRep.Opacity = 0.5  
contFluidVelocityRep.DiffuseColor = [1.0, 0.5, 0.0]  # Orange
ColorBy(contFluidVelocityRep, ('POINTS', 'fluid_velocity'))
contFluidVelocityRep.SetScalarBarVisibility(rv, True)

# ========================== Affichage final ========================== #
ResetCamera()
Render()
Interact()
