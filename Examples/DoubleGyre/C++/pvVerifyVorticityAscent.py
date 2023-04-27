# state file generated using paraview version 5.11.0-RC1-243-g7e6c5885b0
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

renderView1 = GetRenderView()
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [1.0, 0.49803924560546875, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [1.0000000591389835, 0.49803924513980746, 6.700000396231189]
renderView1.CameraFocalPoint = [1.0000000591389835, 0.49803924513980746, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.9652828870844443
renderView1.BackEnd = 'OSPRay raycaster'

# create a new 'XML Rectilinear Grid Reader'
visit_ex_dbvtr = XMLRectilinearGridReader(registrationName='visit_ex_db.vtr', FileName=['/dev/shm/visit_ex_db.vtr'])
visit_ex_dbvtr.PointArrayStatus = ['mesh_mesh/velocity_mag2d', 'mesh_mesh/vel_x', 'mesh_mesh/vel_y']
visit_ex_dbvtr.TimeArray = 'None'

# create a new 'RenameArrays'
renameArrays1 = RenameArrays(registrationName='RenameArrays1', Input=visit_ex_dbvtr)
renameArrays1.PointArrays = ['mesh_mesh/vel_x', 'vel_x', 'mesh_mesh/vel_y', 'vel_y', 'mesh_mesh/velocity_mag2d', 'velocity_mag2d']
renameArrays1.FieldArrays = ['CYCLE', 'CYCLE', 'MeshName', 'MeshName', 'TIME', 'TIME', 'avtOriginalBounds', 'avtOriginalBounds']

# create a new 'Python Calculator'
pythonCalculator1 = PythonCalculator(registrationName='PythonCalculator1', Input=renameArrays1)
pythonCalculator1.Expression = 'make_vector(vel_x, vel_y)'

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=pythonCalculator1,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'result']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 0.050000002956949174
glyph1.GlyphTransform = 'Transform2'

# create a new 'Gradient'
gradient1 = Gradient(registrationName='Gradient1', Input=pythonCalculator1)
gradient1.ScalarArray = ['POINTS', 'result']
gradient1.ComputeGradient = 0
gradient1.ComputeVorticity = 1

gradient1Display = Show(gradient1, renderView1, 'UniformGridRepresentation')

# get 2D transfer function for 'Vorticity'
vorticityTF2D = GetTransferFunction2D('Vorticity')
vorticityTF2D.ScalarRangeInitialized = 1
vorticityTF2D.Range = [0.0, 2.009558925790224, 0.0, 1.0]

# get color transfer function/color map for 'Vorticity'
vorticityLUT = GetColorTransferFunction('Vorticity')
vorticityLUT.TransferFunction2D = vorticityTF2D
vorticityLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.004779462895112, 0.865003, 0.865003, 0.865003, 2.009558925790224, 0.705882, 0.0156863, 0.14902]
vorticityLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Vorticity'
vorticityPWF = GetOpacityTransferFunction('Vorticity')
vorticityPWF.Points = [0.0, 0.0, 0.5, 0.0, 2.009558925790224, 1.0, 0.5, 0.0]
vorticityPWF.ScalarRangeInitialized = 1

gradient1Display.Representation = 'Surface'
gradient1Display.ColorArrayName = ['POINTS', 'Vorticity']
gradient1Display.LookupTable = vorticityLUT

vorticityLUTColorBar = GetScalarBar(vorticityLUT, renderView1)
vorticityLUTColorBar.WindowLocation = 'Upper Right Corner'
vorticityLUTColorBar.Title = 'Vorticity'
vorticityLUTColorBar.ComponentTitle = 'Magnitude'

# set color bar visibility
vorticityLUTColorBar.Visibility = 1

# show color legend
gradient1Display.SetScalarBarVisibility(renderView1, True)

SetActiveSource(gradient1)
