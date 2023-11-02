##############################################################################
# A simple simulator for the heat equation in 2D, with an in-situ coupling
# with the Ascent library https://ascent.readthedocs.io/en/latest/#
#
# Author: Jean M. Favre, Swiss National Supercomputing Center
#
# this version runs in parallel, splitting the domain in the vertical direction
#
# Run: mpiexec -n 2 python3 heat_diffusion_insitu_parallel_Ascent.py \
#                           --res=64 -t 1000 --mesh=uniform
#
# Tested with Python 3.10.12, Tue 12 Sep 16:28:23 CEST 2023
##############################################################################
import sys
import math
import argparse
import numpy as np
import conduit
import conduit.blueprint
import ascent.mpi
import matplotlib.pyplot as plt

from mpi4py import MPI

class Simulation:
    """
    A simple 4-point stencil simulation for the heat equation in serial mode
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
        self.yres = resolution  # redefined when splitting the parallel domain
        self.dx = 1.0 / (self.xres + 1)

    def Initialize(self):
        """ 2 additional boundary points are added. Iterations will only touch
        the internal grid points.
        """
        self.rmesh_dims = [self.yres + 2, self.xres + 2]
        print("Rank ", self.par_rank, ": dimensions = ", self.rmesh_dims)
        self.v = np.zeros(self.rmesh_dims)  # includes 2 ghosts
        self.ghosts = np.zeros(self.rmesh_dims,  dtype=np.ubyte)
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

        self.vnew = 0.25 * (self.v[2:, 1:-1] +  # north neighbor
                            self.v[0:-2, 1:-1] +  # south neighbor
                            self.v[1:-1, 2:] +  # east neighbor
                            self.v[1:-1, :-2])  # west neighbor
        # copy vnew to the interior region of v, leaving the boundary walls untouched.
        self.v[1:-1, 1:-1] = self.vnew.copy()

    def MainLoop(self, _frequency=100):
        while self.iteration < self.Max_iterations:
            self.SimulateOneTimestep()

# define a sub-class of Simulation to add an Ascent in-situ coupling


class ParallelSimulation_With_Ascent(Simulation):
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
        Our aim is to demonstrate the use of Conduit and verify that Ascent can
        process each MPI partition correcly with all 4 grid types.
    verbose : boolean
        prints the Conduit node(s) describing the mesh
    """

    def __init__(self, resolution=64, iterations=100, meshtype="uniform", verbose=False):
        self.comm = MPI.COMM_WORLD
        Simulation.__init__(self, resolution, iterations)
        self.MeshType = meshtype
        self.verbose = verbose

    # Override Initialize for parallel
    def Initialize(self):
        self.par_size = self.comm.Get_size()
        self.par_rank = self.comm.Get_rank()
        # split the parallel domain along the Y axis. No error check!
        self.yres = self.xres // self.par_size
        Simulation.Initialize(self)

        # Add Conduit node and Ascent actions
        # set options to allow errors propagate to python
        ascent_opts = conduit.Node()
        ascent_opts["mpi_comm"] = MPI.COMM_WORLD.py2f()
        ascent_opts["exceptions"] = "forward"
        if self.MeshType in ('uniform', 'rectilinear'):
            ascent_opts["ghost_field_name"] = "point_ghosts"
        else:
            # not implemented yet
            pass
        # open Ascent
        self.a = ascent.mpi.Ascent()
        self.a.open(ascent_opts)

        # setup a mesh
        self.mesh = conduit.Node()

        if self.MeshType == "uniform":
            # create the coordinate set
            self.mesh["coordsets/coords/type"] = self.MeshType
            self.mesh["coordsets/coords/dims/i"] = self.xres + 2
            self.mesh["coordsets/coords/dims/j"] = self.yres + 2
            #self.mesh["coordsets/coords/dims/k"] = 1
            self.mesh["coordsets/coords/origin/x"] = 0.0
            self.mesh["coordsets/coords/origin/y"] = self.par_rank * self.yres * self.dx
            #self.mesh["coordsets/coords/origin/z"] = 0.0
            self.mesh["coordsets/coords/spacing/dx"] = self.dx
            self.mesh["coordsets/coords/spacing/dy"] = self.dx
            #self.mesh["coordsets/coords/spacing/dz"] = self.dx
        elif self.MeshType == "rectilinear":
            xc = np.linspace(0, 1, self.xres + 2)
            y_min = (self.par_rank * self.yres) * self.dx
            y_max = (((self.par_rank + 1) * self.yres) + 1) * self.dx
            yc = np.linspace(y_min, y_max, self.yres + 2)
            self.mesh["coordsets/coords/type"] = self.MeshType
            self.mesh["coordsets/coords/values/x"].set_external(xc)
            self.mesh["coordsets/coords/values/y"].set_external(yc)

        else: # self.MeshType in ('structured', 'unstructured'):
            y_min = (self.par_rank * self.yres) * self.dx
            y_max = (((self.par_rank + 1) * self.yres) + 1) * self.dx
            self.xc, self.yc = np.meshgrid(np.linspace(0, 1, self.xres + 2),
                                           np.linspace(y_min, y_max, self.yres + 2),
                                           indexing='xy')
            self.mesh["coordsets/coords/type"] = "explicit"
            self.mesh["coordsets/coords/values/x"].set_external(self.xc.ravel())
            self.mesh["coordsets/coords/values/y"].set_external(self.yc.ravel())

        if self.MeshType == "structured":
            self.mesh["topologies/mesh/elements/dims/i"] = np.int32(self.xres + 1)
            self.mesh["topologies/mesh/elements/dims/j"] = np.int32(self.yres + 1)

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
            self.mesh["topologies/mesh/elements/shape"] = "quad"
            self.mesh["topologies/mesh/elements/connectivity"].set_external(self.connectivity)

        # add the topology, implicitly derived from the coordinate set
        self.mesh["topologies/mesh/type"] = self.MeshType
        self.mesh["topologies/mesh/coordset"] = "coords"

        # create a vertex associated field called "temperature"
        self.mesh["fields/temperature/association"] = "vertex"
        self.mesh["fields/temperature/topology"] = "mesh"
        # set_external does not handle multidimensional numpy arrays or
        # multidimensional complex strided views into numpy arrays.
        # Views that are effectively 1D-strided are supported.
        self.mesh["fields/temperature/values"].set_external(self.v.ravel())

        if self.MeshType in ('uniform', 'rectilinear'):
            # create a vertex associated field called "point_ghosts"
            self.mesh["fields/point_ghosts/association"] = "vertex"
            self.mesh["fields/point_ghosts/topology"] = "mesh"
            self.mesh["fields/point_ghosts/values"].set_external(self.ghosts.ravel())

        # make sure the mesh we created conforms to the blueprint
        verify_info = conduit.Node()
        if not conduit.blueprint.mesh.verify(self.mesh, verify_info):
            print("Mesh Verify failed!")
            print(verify_info.to_yaml())
        else:
            if self.verbose:
                print(self.mesh)

        # setup actions
        self.actions = conduit.Node()
        add_act = self.actions.append()
        add_act["action"] = "add_scenes"

        # declare a scene (s1) to render the dataset
        self.scenes = add_act["scenes"]
        self.scenes["s1/plots/p1/type"] = "pseudocolor"
        self.scenes["s1/plots/p1/field"] = "temperature"
        # add a second plot to draw the grid lines
        self.scenes["s1/plots/p2/type"] = "mesh"

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

    def MainLoop(self, frequency=100):
        while self.iteration < self.Max_iterations:
            self.SimulateOneTimestep()
            if not self.iteration % frequency:
                self.mesh["state/cycle"] = self.iteration
                self.mesh["state/time"] = self.iteration * 0.1
                self.mesh["state/title"] = "2D Heat diffusion simulation"
                self.mesh["state/info"] = "In-situ pseudocolor rendering of temperature"

                self.scenes["s1/renders/r1/image_name"] = "temperature-par.%04d" % self.iteration
                # execute the actions
                self.a.publish(self.mesh)
                self.a.execute(self.actions)

    def Finalize(self, savedir="./"):
        """ After the final timestep, we save the solution array to disk
        and we close Ascent"""
        self.a.publish(self.mesh)
        action = conduit.Node()
        add_extr = action.append()
        add_extr["action"] = "add_extracts"
        extracts = add_extr["extracts"]
        extracts["e1/type"] = "relay"
        extracts["e1/params/path"] = savedir + "mesh"
        extracts["e1/params/protocol"] = "blueprint/mesh/hdf5"
        self.a.execute(action)
        self.a.close()


def main(args):
    # run without in-situ coupling and without MPI
    if not args.noinsitu:
        sim0 = Simulation(resolution=args.res, iterations=args.timesteps)
        sim0.Initialize()
        sim0.MainLoop()
        sim0.Finalize()
    else:
        # run with in-situ Ascent coupling and with MPI
        # meshtype can be one of "uniform", "rectilinear", "structured", "unstructured"
        sim = ParallelSimulation_With_Ascent(resolution=args.res,
                                             meshtype=args.mesh,
                                             iterations=args.timesteps,
                                             verbose=args.verbose)
        sim.Initialize()
        sim.MainLoop(frequency=args.frequency)
        sim.Finalize(savedir=args.dir)


parser = argparse.ArgumentParser(
    description="heat diffusion miniapp for in-situ visualization testing with Ascent")
parser.add_argument("-t", "--timesteps", type=int,
                    help="number of timesteps to run the miniapp (default: 1000)",
                    default=1000)
parser.add_argument("--res", type=int,
                    help="resolution in each coordinate direction (default: 64)",
                    default=64)
parser.add_argument("-m", "--mesh", type=str, default="uniform",
                    choices=["uniform", "rectilinear", "structured", "unstructured"],
                    help="mesh type (default: uniform)")
parser.add_argument("-f", "--frequency", type=int, default=50,
                    help="How often should the Ascent script be executed in situ processing.")
parser.add_argument("-d", "--dir", type=str,
                    help="path to a directory where to dump the Blueprint output",
                    default=".")
parser.add_argument("-n", "--noinsitu",
                    help="toggle the use of the in-situ vis coupling with Ascent",
                    action='store_false')  # on/off flag)
parser.add_argument("-v", "--verbose",
                    help="toggle printing of the conduit nodes",
                    action='store_true')  # on/off flag

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
