# script-version: 2.0
# Catalyst state generated using paraview version 5.10.0-RC1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = GetRenderView()
renderView1.ViewSize = [512,512]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [10.5, 10.5, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [10.5, 10.5, 10000.0]
renderView1.CameraFocalPoint = [10.5, 10.5, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 16.252051019206938

reader = TrivialProducer(registrationName='grid')

rep = Show(reader, renderView1)
rep.Representation = 'Surface With Edges'
ColorBy(rep, ['POINTS', 'temperature'])

ResetCamera()

pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
pNG1.Trigger = 'TimeStep'
pNG1.Trigger.Frequency = 50
pNG1.Writer.FileName = '/tmp/view-{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [512, 512]
pNG1.Writer.Format = 'PNG'

# create extractor
vTP1 = CreateExtractor('VTS', reader, registrationName='VTS1')
vTP1.Trigger = 'TimeStep'
vTP1.Trigger.Frequency = 500
vTP1.Writer.FileName = 'grid{timestep:06d}.pvts'

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
