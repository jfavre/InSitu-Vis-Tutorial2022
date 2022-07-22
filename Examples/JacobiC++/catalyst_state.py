#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = GetRenderView()
renderView1.ViewSize = [800,800]
renderView1.CenterOfRotation = [0.5, 0.5, 0.0]
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 16.252051019206938

LoadPalette(paletteName='WhiteBackground')

reader = TrivialProducer(registrationName='grid')
rep = Show(reader, renderView1)
rep.Representation = 'Outline'
ColorBy(rep, ['POINTS', ''])

contour1 = Contour(registrationName='Contour1', Input=reader)
contour1.ContourBy = ['POINTS', 'temperature']
contour1.ComputeNormals = 0
contour1.ComputeScalars = 1
contour1.Isosurfaces = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
contour1.PointMergeMethod = 'Uniform Binning'

pid = ProcessIdScalars(Input=contour1)
pid.UpdatePipeline()

contour1Display = Show(pid, renderView1)
contour1Display.LineWidth = 2
ColorBy(contour1Display, ['POINTS', 'ProcessId'])
processIdLUT = GetColorTransferFunction('ProcessId')
processIdLUT.RescaleTransferFunction(0.0, 3.0)

ResetCamera()

pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
pNG1.Trigger = 'TimeStep'
pNG1.Trigger.Frequency = 200
pNG1.Writer.FileName = '/tmp/view-{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [800,800]
pNG1.Writer.Format = 'PNG'

# The vtkPartitionedDataSet type works for all our supported types

vTP1 = CreateExtractor('VTPD', reader, registrationName='VTPD1')
vTP1.Trigger = 'TimeStep'
vTP1.Trigger.Frequency = 500
vTP1.Writer.FileName = 'dataset_{timestep:06d}.vtpd'


SetActiveSource(reader)

from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.EnableCatalystLive = 0
options.CatalystLiveURL = ':22222'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
