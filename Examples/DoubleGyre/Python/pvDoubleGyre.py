
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1280, 768]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [1.0, 0.5, 0.0]
renderView1.CameraPosition = [1.0, 0.5, 6.7]
renderView1.CameraFocalPoint = [1.0, 0.5, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.8546757251683884
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

reader = TrivialProducer(registrationName='grid')
readerDisplay = Show(reader, renderView1, 'GeometryRepresentation')
readerDisplay.Representation = 'Outline'

glyph1 = Glyph(registrationName='Glyph1', Input=reader, GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'Velocity']
glyph1.ScaleArray = ['POINTS', 'Velocity']
#glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 0.6
glyph1.GlyphMode = 'Uniform Spatial Distribution (Surface Sampling)'
glyph1.MaximumNumberOfSamplePoints = 500

glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
ColorBy(glyph1Display, ['POINTS', 'Velocity'])

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=reader)

annotateTimeFilter1Display = Show(annotateTimeFilter1)
annotateTimeFilter1Display.Bold = 1
annotateTimeFilter1Display.FontSize = 24

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
pNG1.Trigger = 'TimeStep'
pNG1.Trigger.Frequency = 2
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1280, 768]
pNG1.Writer.Format = 'PNG'

# create extractor
#vTR1 = CreateExtractor('VTI', reader, registrationName='VTI1')
vTR1 = CreateExtractor('VTPD', reader, registrationName='VTPD1')
vTR1.Trigger = 'TimeStep'
vTR1.Trigger.Frequency = 10
vTR1.Writer.FileName = 'doublegyre_{timestep:06d}.vtpd'

# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
