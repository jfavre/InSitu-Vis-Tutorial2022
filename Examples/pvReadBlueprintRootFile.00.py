# a prototype of an HDF5 blueprint root file reader for ParaView
#
# https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#mesh-index-protocol
#
# tested with Paraview v5.11.2 Wed 25 Oct 16:58:38 CEST 2023
#
# Written by: Jean M, Favre, Swiss National Supercomputing Center
##############################################################################
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = GetRenderView()
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.5, 0.5, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.5, 0.5, 2.7320508075688776]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.7071067811865476
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

programmableSource1 = ProgrammableSource()
programmableSource1.OutputDataSetType = 'vtkPartitionedDataSet'
programmableSource1.Script = """
import conduit
import conduit.relay.io
import numpy as np
from vtk import vtkPoints, vtkImageData, vtkRectilinearGrid, vtkStructuredGrid
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs

def AddArrays(grid, mesh, nnodes):
  it = conduit.NodeIterator()
  it = mesh["mesh/fields"].children()
  while(it.has_next()):
    f = it.next()
    node = it.node()
    nparr = node["values"]
    nparr2 = np.copy(nparr) # forces the read to take place
    dArray = dsa.numpyTovtkDataArray(nparr2, it.name())
    if node["association"] == "vertex":
      assert nparr.shape[0] == nnodes
      grid.GetPointData().AddArray(dArray)
    else:
      grid.GetCellData().AddArray(dArray)
        
def ConduitNode_to_StructuredGrid(mesh):
  assert mesh["mesh/topologies/mesh/type"] == "structured"
  grid = vtkStructuredGrid()
  xcoords = mesh["mesh/coordsets/coords/values/x"];
  ycoords = mesh["mesh/coordsets/coords/values/y"];
  dims    = [mesh["mesh/topologies/mesh/elements/dims/i"]+1,
            mesh["mesh/topologies/mesh/elements/dims/j"]+1,
            1]
  if mesh.has_path("mesh/topologies/mesh/elements/dims/z"): # this is a 3D grid
    dims[2] = mesh["mesh/topologies/mesh/elements/dims/z"]+1
    zcoords = mesh["mesh/coordsets/coords/values/z"]
    coords = algs.make_vector(np.copy(xcoords), np.copy(ycoords), np.copy(zcoords))
  else:
    coords = algs.make_vector(np.copy(xcoords), np.copy(ycoords), np.zeros_like(xcoords))

  grid.SetDimensions(dims)
  points = vtkPoints()
  grid.SetPoints(points)
  points.SetData(dsa.numpyTovtkDataArray(coords, "coords"))

  nnodes = np.prod(dims)
  if mesh.has_path("mesh/fields"):
    AddArrays(grid, mesh, nnodes)

  return grid

def Make_structured_grid(fname):
  mesh = conduit.Node()
  conduit.relay.io.load(mesh, fname, "hdf5")
  grid = ConduitNode_to_StructuredGrid(mesh)
  return grid
  
def ConduitNode_to_RectilinearGrid(mesh):
  assert mesh["mesh/topologies/mesh/type"] == "rectilinear"
  grid = vtkRectilinearGrid()
  xcoords = mesh["mesh/coordsets/coords/values/x"]
  ycoords = mesh["mesh/coordsets/coords/values/y"]
  dims    = [xcoords.shape[0], ycoords.shape[0], 1]
  if mesh.has_path("mesh/coordsets/coords/values/z"): # this is a 3D grid
    zcoords = mesh["mesh/coordsets/coords/values/z"]
    dims[2] = zcoords.shape[0]
  else:
    zcoords = np.array([0.0])

  grid.SetDimensions(dims)
  grid.SetXCoordinates(dsa.numpyTovtkDataArray(np.copy(xcoords), "xcoords"))
  grid.SetYCoordinates(dsa.numpyTovtkDataArray(np.copy(ycoords), "ycoords"))
  grid.SetZCoordinates(dsa.numpyTovtkDataArray(np.copy(zcoords), "zcoords"))
  
  nnodes = np.prod(dims)              
  if mesh.has_path("mesh/fields"):
    AddArrays(grid, mesh, nnodes)

  return grid
  
def Make_rectilinear_grid(fname):
  mesh = conduit.Node()
  conduit.relay.io.load(mesh, fname, "hdf5")
  grid = ConduitNode_to_RectilinearGrid(mesh)
  return grid

def ConduitNode_to_UniformGrid(mesh):
  assert mesh["mesh/topologies/mesh/type"] == "uniform"
  grid = vtkImageData()
  # initialize for a 2D grid and check later if a 3D grid is present
  origin  = [mesh["mesh/coordsets/coords/origin/x"],
             mesh["mesh/coordsets/coords/origin/y"],
             0.]
  spacing = [mesh["mesh/coordsets/coords/spacing/dx"],
             mesh["mesh/coordsets/coords/spacing/dy"],
             1.0]
  dims    = [mesh["mesh/coordsets/coords/dims/i"],
             mesh["mesh/coordsets/coords/dims/j"],
             1]
  if mesh.has_path("mesh/coordsets/coords/dims/k"): # this is a 3D grid
    origin[2]  = mesh["mesh/coordsets/coords/origin/z"]
    spacing[2] = mesh["mesh/coordsets/coords/spacing/dz"]
    dims[2]    = mesh["mesh/coordsets/coords/dims/k"]

  grid.SetOrigin(origin)
  grid.SetSpacing(spacing)
  grid.SetDimensions(dims)

  nnodes = np.prod(dims)              
  if mesh.has_path("mesh/fields"):
    AddArrays(grid, mesh, nnodes)
  return grid
  
def Make_uniform_grid(fname):
  mesh = conduit.Node()
  conduit.relay.io.load(mesh, fname, "hdf5")
  grid = ConduitNode_to_UniformGrid(mesh)
  return grid
  
root = conduit.Node()
basename = "/dev/shm/"
conduit.relay.io.load(root, basename + "mesh.cycle_001000.root", "hdf5")
number_of_domains = root["blueprint_index/mesh/state/number_of_domains"]
topotype = root["blueprint_index/mesh/topologies/mesh/type"]
if topotype == "uniform":
  if number_of_domains == 1:
    grid = ConduitNode_to_UniformGrid(root)
    output.SetPartition(0, grid)
  else:
    for domain in range(number_of_domains):
      print("domain ", domain, " in ", basename + format(root["file_pattern"] % domain))
      grid = Make_uniform_grid(format(basename + root["file_pattern"] % domain))
      output.SetPartition(domain, grid)
elif topotype == "rectilinear":
  if number_of_domains == 1:
    grid = ConduitNode_to_RectilinearGrid(root)
    output.SetPartition(0, grid)
  else:
    for domain in range(number_of_domains):
      print("domain ", domain, " in ", basename + format(root["file_pattern"] % domain))
      grid = Make_rectilinear_grid(format(basename + root["file_pattern"] % domain))
      output.SetPartition(domain, grid)
elif topotype == "structured":
  if number_of_domains == 1:
    grid = ConduitNode_to_StructuredGrid(root)
    output.SetPartition(0, grid)
  else:
    for domain in range(number_of_domains):
      print("domain ", domain, " in ", basename + format(root["file_pattern"] % domain))
      grid = Make_structured_grid(format(basename + root["file_pattern"] % domain))
      output.SetPartition(domain, grid)
else:
  print("Mesh type not implemented yet")
"""

programmableSource1Display = Show(programmableSource1)
programmableSource1Display.Representation = 'Surface'
ColorBy(programmableSource1Display, ['POINTS', 'temperature'])

