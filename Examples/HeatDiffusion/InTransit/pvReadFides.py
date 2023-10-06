# state file generated using paraview version 5.11.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

from paraview import print_info

from paraview.catalyst import Options
options = Options()

# directory under which to save all extracts
options.ExtractsOutputDirectory = '/tmp'
SaveExtractsUsingCatalystOptions(options)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
view = GetRenderView()
view.InteractionMode = '2D'
view.AxesGrid = 'GridAxes3DActor'
view.CenterOfRotation = [0.5, 0.5, 0.0]
view.StereoType = 'Crystal Eyes'
view.CameraPosition = [0.5, 0.5, 3.35]
view.CameraFocalPoint = [0.5, 0.5, 0.0]
view.CameraFocalDisk = 1.0
view.CameraParallelScale = 0.7071067811865476
view.BackEnd = 'OSPRay raycaster'
view.OSPRayMaterialLibrary = materialLibrary1

if __name__ == "__main__":
  producer = FidesReader(registrationName='diffusion.bp', FileName='diffusion.bp')
else:
  producer = TrivialProducer(registrationName='fides')

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=producer)
contour1.ContourBy = ['POINTS', 'temperature']
contour1.ComputeNormals = 0
contour1.Isosurfaces = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
contour1.PointMergeMethod = 'Uniform Binning'
  
producerDisplay = Show(producer)
producerDisplay.Representation = 'Surface'
ColorBy(producerDisplay, ('POINTS', 'temperature'))
# get color transfer function/color map for 'temperature'
temperatureLUT = GetColorTransferFunction('temperature')
temperatureLUT.RescaleTransferFunction(0.1, 0.9)

contour1Display = Show(contour1, view, 'GeometryRepresentation')
ColorBy(contour1Display, ('POINTS', 'temperature'))
contour1Display.LineWidth = 2.0

temperatureLUTColorBar = GetScalarBar(temperatureLUT, view)
temperatureLUTColorBar.Title = 'temperature'
temperatureLUTColorBar.ComponentTitle = ''

# set color bar visibility
temperatureLUTColorBar.Visibility = 1

contour1Display.SetScalarBarVisibility(view, True)

SetActiveSource(contour1)

# the extractor will save the view on each time step
extractor = CreateExtractor('PNG', view, registrationName='PNG1')
extractor.Writer.FileName = 'output-{timestep}.png'
extractor.Writer.ImageResolution = [800, 800]

if __name__ == "__main__":
  GetAnimationScene().UpdateAnimationUsingDataTimeSteps()

def catalyst_execute(info):
    global producer
    producer.UpdatePipeline()
