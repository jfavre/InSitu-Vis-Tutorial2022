##############################################################################
# A simple simulator for the heat equation in 2D, with an in-situ coupling
# using Catalyst https://catalyst-in-situ.readthedocs.io/en/latest/index.html
#
# Author: Jean M. Favre, Swiss National Supercomputing Center
#
# Can run in parallel (if in-situ is turned on), splitting the domain in
# the vertical direction.
# There is no error checking on grid resolution and MPI domain splitting.
# it is strongly advised to use grid resolutions like 2^N and
# an even number of MPI partitions, e.g.
#
# Run: mpiexec -n 4 python3 heat_diffusion_insitu_parallel_Catalyst.py \
#                           --res=64 -t 1000 --script ../C++/catalyst_state.py
#
# Tested with Python 3.10.12, Mon  2 Oct 08:16:07 CEST 2023
#
##############################################################################
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import catalyst
import catalyst_conduit as conduit
import catalyst_conduit.blueprint

from mpi4py import MPI


class Simulation:
    """
    A simple 4-point stencil simulation for the heat equation.
    The domain (X, Y) is [0.0, 1.0] x [0.0, 1.0]

    Attributes
    ----------
    resolution : int
        the number of grid points on the I and J axis (default 64)
    iterations : int
        the maximum number of iterations (default 100)
    """
    def __init__(self, resolution=64, iterations=100):
        self.par_size = 1
        self.par_rank = 0
        self.iteration = 0  # current iteration
        self.Max_iterations = iterations
        self.xres = resolution
        self.yres = resolution  # is redefined when splitting the parallel domain
        self.dx = 1.0 / (self.xres + 1)

    def Initialize(self):
        """ 2 additional boundary points are added. Iterations will only touch
        the internal grid points.
        """
        self.rmesh_dims = [self.yres + 2, self.xres + 2]
        self.v = np.zeros(self.rmesh_dims)
        self.ghosts = np.zeros(self.rmesh_dims, dtype=np.ubyte)
        self.vnew = np.zeros([self.yres, self.xres])
        self.set_initial_bc()

    def set_initial_bc(self):
        if self.par_size > 1:
            if self.par_rank == 0:
                self.v[0, :] = [math.sin(math.pi * j * self.dx)
                               for j in range(self.rmesh_dims[1])]
                self.ghosts[-1, :] = 1
            elif self.par_rank == (self.par_size - 1):
                self.v[-1, :] = self.v[0, :] * math.exp(-math.pi)
                self.ghosts[0, :] = 1
            else:
                self.ghosts[0, :] = 1
                self.ghosts[-1, :] = 1
        else:
            # first (bottom) row
            self.v[0, :] = [math.sin(math.pi * j * self.dx)
                           for j in range(self.rmesh_dims[1])]
            # last (top) row
            self.v[-1, :] = self.v[0, :] * math.exp(-math.pi)

    def Finalize(self):
        """plot the scalar field iso-contour lines"""
        fig, ax = plt.subplots()
        CS = ax.contour(self.v, levels=10)
        ax.clabel(CS, inline=True, fontsize=10)
        ax.set_title('Temperature iso-contours')
        fname = f'Temperature-iso-contours.{self.iteration:04d}.png'
        plt.savefig(fname)
        print("Final image \"", fname, "\" written to disk", sep="")

    def SimulateOneTimestep(self):
        # there is no ghost-data exchange. Run in serial-mode only
        self.iteration += 1

        self.vnew = 0.25 * ( self.v[2:, 1:-1]  +  # north neighbor
                             self.v[0:-2, 1:-1] +  # south neighbor
                             self.v[1:-1, 2:] +  # east neighbor
                             self.v[1:-1, :-2])  # west neighbor
        # copy vnew to the interior region of v, leaving the boundary walls untouched.
        self.v[1:-1, 1:-1] = self.vnew.copy()

    def MainLoop(self):
        while self.iteration < self.Max_iterations:
            self.SimulateOneTimestep()

# define a sub-class of Simulation to add a Catalyst in-situ coupling


class ParallelSimulation_With_Catalyst(Simulation):
    """
    Attributes
    ----------
    resolution : int
        the number of grid points on the I and J axis (default 64)
    iterations : int
        the maximum number of iterations (default 100)
    meshtype : string
        can be one of "uniform", "rectilinear", "structured", "unstructured"
        this is for demonstration purposes only. The computation itself is
        independent of the underlying grid since it uses a simple 4-point stencil.
        Our aim is to demonstrate the use of Conduit and verify that Catalyst can process
        each MPI partition correcly with all 4 grid types.
    pv_script : string
        a ParaView Catalyst script file to generate images and other visualization outputs
    verbose : boolean
        prints the Conduit node(s) describing the mesh
    """

    def __init__(self, resolution=64, iterations=100, meshtype="uniform", pv_script="catalyst_state.py", verbose=False):
        self.comm = MPI.COMM_WORLD
        Simulation.__init__(self, resolution, iterations)
        self.MeshType = meshtype

        self.insitu = conduit.Node()
        self.pv_script = pv_script
        self.verbose = verbose
    # Add Catalyst mesh definition

    def Initialize(self):
        self.par_size = self.comm.Get_size()
        self.par_rank = self.comm.Get_rank()
        # split the parallel domain along the Y axis. No error check!
        self.yres = self.xres // self.par_size

        Simulation.Initialize(self)
        self.initialize_catalyst()
        
        self.exec_params = conduit.Node()
        channel = self.exec_params["catalyst/channels/grid"]
        channel["type"] = "mesh"
        mesh = channel["data"]

        # create the coordinate set
        if self.MeshType == "rectilinear":
            xc = np.linspace(0, 1, self.xres + 2)
            y_min = (self.par_rank * self.yres) * self.dx
            y_max = (((self.par_rank + 1) * self.yres) + 1) * self.dx
            yc = np.linspace(y_min, y_max, self.yres + 2)
            mesh["coordsets/coords/type"] = self.MeshType
            mesh["coordsets/coords/values/x"].set_external(xc)
            mesh["coordsets/coords/values/y"].set_external(yc)
            
        elif self.MeshType == "uniform":
            mesh["coordsets/coords/type"] = self.MeshType
            mesh["coordsets/coords/dims/i"] = self.xres + 2
            mesh["coordsets/coords/dims/j"] = self.yres + 2
            mesh["coordsets/coords/origin/x"] = 0.0
            mesh["coordsets/coords/origin/y"] = self.par_rank * self.yres * self.dx
            mesh["coordsets/coords/spacing/dx"] = self.dx
            mesh["coordsets/coords/spacing/dy"] = self.dx
            
        else: # self.MeshType in ('structured', 'unstructured'):
            y_min = (self.par_rank * self.yres) * self.dx
            y_max = (((self.par_rank + 1) * self.yres) + 1)  * self.dx
            self.xc, self.yc = np.meshgrid(np.linspace(0, 1, self.xres + 2),
                                           np.linspace(y_min, y_max, self.yres + 2),
                                           indexing='xy')
            mesh["coordsets/coords/type"] = "explicit"
            mesh["coordsets/coords/values/x"].set_external(self.xc.ravel())
            mesh["coordsets/coords/values/y"].set_external(self.yc.ravel())
            
        mesh["topologies/mesh/type"] = self.MeshType
        mesh["topologies/mesh/coordset"] = "coords"

        if self.MeshType == "structured":
            mesh["topologies/mesh/elements/dims/i"] = np.int32(self.xres + 1)
            mesh["topologies/mesh/elements/dims/j"] = np.int32(self.yres + 1)

        if self.MeshType == "unstructured":
            # we describe the grid as a list of VTK 2D Quadrilaterals
            nbOfQuads = (self.xres + 1) * (self.yres + 1) * 4
            self.connectivity = np.zeros(nbOfQuads, dtype=np.int32)
            i = 0
            for iy in range(self.yres + 1):
                for ix in range(self.xres + 1):
                    self.connectivity[4 * i + 0] = ix + iy * (self.xres + 2)
                    self.connectivity[4 * i + 1] = ix + (iy + 1) * (self.xres + 2)
                    self.connectivity[4 * i + 2] = ix + (iy + 1) * (self.xres + 2) + 1
                    self.connectivity[4 * i + 3] = ix + iy * (self.xres + 2) + 1
                    i += 1
            mesh["topologies/mesh/elements/shape"] = "quad"
            mesh["topologies/mesh/elements/connectivity"].set_external(self.connectivity)

        # create a vertex associated field called "temperature"
        mesh["fields/temperature/association"] = "vertex"
        mesh["fields/temperature/topology"] = "mesh"
        # set_external does not handle multidimensional numpy arrays or
        # multidimensional complex strided views into numpy arrays.
        # Views that are effectively 1D-strided are supported.
        mesh["fields/temperature/values"].set_external(self.v.ravel())

        if self.par_size > 1: # create a vertex associated field called "point_ghosts"
            mesh["fields/vtkGhostType/association"] = "vertex"
            mesh["fields/vtkGhostType/topology"] = "mesh"
            mesh["fields/vtkGhostType/values"].set_external(self.ghosts.ravel())

        # make sure the mesh we created conforms to the blueprint
        verify_info = conduit.Node()
        if not conduit.blueprint.mesh.verify(mesh, verify_info):
            print("Heat mesh verify failed!")
        else:
            if self.verbose:
                print(mesh)

    def SimulateOneTimestep(self):
        Simulation.SimulateOneTimestep(self)

        if self.par_size > 1:
            # if in parallel, exchange ghost cells now
            # define who is my neighbor above and below
            below = self.par_rank - 1
            above = self.par_rank + 1
            if self.par_rank == 0:
                below = MPI.PROC_NULL   # tells MPI not to perform send/recv
            if self.par_rank == (self.par_size - 1):
                above = MPI.PROC_NULL   # should only receive/send from/to below
            self.comm.Sendrecv([self.v[-2,], self.xres + 2, MPI.DOUBLE], dest=above,
                               recvbuf=[self.v[-0,], self.xres + 2, MPI.DOUBLE], source=below)
            self.comm.Sendrecv([self.v[1,], self.xres + 2, MPI.DOUBLE], dest=below,
                               recvbuf=[self.v[-1,], self.xres + 2, MPI.DOUBLE], source=above)

    def MainLoop(self):
        while self.iteration < self.Max_iterations:
            self.SimulateOneTimestep()

            state = self.exec_params["catalyst/state"]
            state["timestep"] = self.iteration
            state["time"] = self.iteration * 0.1
            catalyst.execute(self.exec_params)

    def initialize_catalyst(self):
        """Creates a Conduit node """
        self.insitu["catalyst/scripts/script/filename"] = self.pv_script
        self.insitu["catalyst_load/implementation"] = "paraview"

        # open Catalyst
        catalyst.initialize(self.insitu)

    def finalize_catalyst(self):
        """close"""
        catalyst.finalize(self.insitu)


def main(args):
    # run without in-situ Catalyst coupling and without MPI
    if not args.noinsitu:
        sim0 = Simulation(resolution=args.res, iterations=args.timesteps)
        sim0.Initialize()
        sim0.MainLoop()
        sim0.Finalize()
    else:
        # run with in-situ Catalyst coupling and with MPI
        # meshtype can be one of "uniform", "rectilinear", "structured", "unstructured"
        sim = ParallelSimulation_With_Catalyst(resolution=args.res,
                                               meshtype=args.mesh,
                                               iterations=args.timesteps,
                                               pv_script=args.script,
                                               verbose=args.verbose)
        sim.Initialize()
        sim.MainLoop()
        sim.finalize_catalyst()


parser = argparse.ArgumentParser(
    description="heat diffusion miniapp for ParaView Catalyst (v2) testing")
parser.add_argument("-t", "--timesteps", type=int,
                    help="number of timesteps to run the miniapp (default: 1000)",
                    default=1000)
parser.add_argument("--res", type=int,
                    help="resolution in each coordinate direction (default: 64)", default=64)
parser.add_argument("-m", "--mesh", type=str, default="uniform",
                    choices=["uniform", "rectilinear", "structured", "unstructured"],
                    help="mesh type (default: uniform)")
parser.add_argument("-s", "--script", type=str,
                    help="path to the Catalyst script to use for in situ processing.",
                    default="../C++/catalyst_state.py")
parser.add_argument("-n", "--noinsitu",
                    help="toggle the use of the in-situ vis coupling with Catalyst",
                    action='store_false')  # on/off flag)
parser.add_argument("-v", "--verbose",
                    help="toggle printing of the conduit nodes",
                    action='store_true')  # on/off flag

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
