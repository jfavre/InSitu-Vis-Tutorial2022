##############################################################################
# A simple simulator for the heat equation in 2D, with an in-situ coupling
# with the Ascent library https://ascent.readthedocs.io/en/latest/#
#
# Author: Jean M. Favre, Swiss National Supercomputing Center
#
# this version runs in parallel, splitting the domain in the vertical direction
#
# Run: mpiexec -n 2 python3 jacobi_insitu_parallel_Ascent.py
#
# Tested with Python 3.9.7, Wed Apr  6 11:16:53 AM CEST 2022
##############################################################################
import sys, math
import numpy as np
import conduit
import conduit.blueprint
import ascent.mpi

from mpi4py import MPI

# print details about ascent
print(ascent.mpi.about())

class Simulation:
    """
    A simple 4-point stencil simulation for the heat equation
    
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
        self.iteration = 0 # current iteration
        self.Max_iterations = iterations
        self.xres = resolution
        # self.yres is redefined when splitting the parallel domain
        self.dx = 1.0 / (self.xres + 1)

    def Initialize(self):
        """ 2 additional boundary points are added. Iterations will only touch
        the internal grid points.
        A mesh description using Conduit is prepared, sharing the scalar field
        array via zero-copy
        """
        self.rmesh_dims = [self.yres + 2, self.xres + 2]
        print("dimensions = ", self.rmesh_dims)
        self.v = np.zeros(self.rmesh_dims) # includes 2 ghosts
        #self.v = np.ones(self.rmesh_dims) * self.par_rank
        self.vnew = np.zeros([self.yres, self.xres])
        self.set_initial_bc()

    def set_initial_bc(self):
        if self.par_size > 1:
          if self.par_rank == 0:
            self.v[0,:] = [math.sin(math.pi * j * self.dx)
                           for j in range(self.rmesh_dims[1])]
          if self.par_rank == (self.par_size-1):
            self.v[-1,:] = self.v[0,:] * math.exp(-math.pi)
        else:
          #first (bottom) row
          self.v[0,:] = [math.sin(math.pi * j * self.dx)
                         for j in range(self.rmesh_dims[1])]
          #last (top) row
          self.v[-1,:] = self.v[0,:] * math.exp(-math.pi)

    def Finalize(self):
        pass

    def SimulateOneTimestep(self):
        self.iteration += 1

        self.vnew = 0.25 * ( self.v[2:, 1:-1]  + # north neighbor
                             self.v[0:-2, 1:-1] + # south neighbor
                             self.v[1:-1, 2:] + # east neighbor
                             self.v[1:-1, :-2]) # west neighbor
        # copy now vnew to the interior region of v, leaving the boundary walls untouched.
        self.v[1:-1,1:-1] = self.vnew.copy()

    def MainLoop(self, _frequency=100):
      while self.iteration < self.Max_iterations:
        self.SimulateOneTimestep()

class ParallelSimulation_With_Ascent(Simulation):
    def __init__(self, resolution, iterations):
        self.comm = MPI.COMM_WORLD
        Simulation.__init__(self, resolution, iterations)

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
        # open Ascent
        self.a = ascent.mpi.Ascent()
        self.a.open(ascent_opts)

        # setup a uniform mesh
        self.mesh = conduit.Node()

        # create the coordinate set
        self.mesh["coordsets/coords/type"] = "uniform";
        self.mesh["coordsets/coords/dims/i"] = self.xres + 2
        self.mesh["coordsets/coords/dims/j"] = self.yres + 2
        
        # add origin and spacing to the coordset (optional)
        self.mesh["coordsets/coords/origin/x"] = 0.0
        self.mesh["coordsets/coords/origin/y"] = self.par_rank * self.yres * self.dx

        self.mesh["coordsets/coords/spacing/dx"] = self.dx
        self.mesh["coordsets/coords/spacing/dy"] = self.dx
        
        # add the topology, implicitly derived from the coordinate set
        self.mesh["topologies/mesh/type"] = "uniform";
        self.mesh["topologies/mesh/coordset"] = "coords";
        
        # create a vertex associated field called "temperature"
        self.mesh["fields/temperature/association"] = "vertex";
        self.mesh["fields/temperature/topology"] = "mesh";
        # set_external does not handle multidimensional numpy arrays or
        # multidimensional complex strided views into numpy arrays.
        # Views that are effectively 1D-strided are supported.
        self.mesh["fields/temperature/values"].set_external(self.v.ravel())

        # make sure the mesh we created conforms to the blueprint
        verify_info = conduit.Node()
        if not conduit.blueprint.mesh.verify(self. mesh,verify_info):
          print("Mesh Verify failed!")
          print(verify_info.to_yaml())
        else:
          print("Mesh verify success!")
    
        # print the mesh we created
        # print(self.mesh.to_yaml())

        # setup actions
        self.actions = conduit.Node()
        add_act = self.actions.append()
        add_act["action"] = "add_scenes"

        # declare a scene (s1) to render the dataset
        self.scenes = add_act["scenes"]
        # our first scene (named 's1') will render the field 'temperature'
        self.scenes["s1/plots/p1/type"] = "pseudocolor";
        self.scenes["s1/plots/p1/field"] = "temperature";
        # add a second plot to draw the grid lines
        self.scenes["s1/plots/p2/type"] = "mesh";

    def SimulateOneTimestep(self):
        Simulation.SimulateOneTimestep(self)

        if self.par_size > 1:
          # if in parallel, exchange ghost cells now
          # define who is my neighbor above and below
          below = self.par_rank - 1
          above = self.par_rank + 1
          if self.par_rank == 0:
            below = MPI.PROC_NULL   # tells MPI not to perform send/recv
          if self.par_rank == (self.par_size-1):
            above = MPI.PROC_NULL   # should only receive/send from/to below
          self.comm.Sendrecv([self.v[-2,], self.xres + 2, MPI.DOUBLE],
                             dest=above, recvbuf=[self.v[-0,], self.xres + 2, MPI.DOUBLE], source=below)
          self.comm.Sendrecv([self.v[1,], self.xres + 2, MPI.DOUBLE],
                             dest=below, recvbuf=[self.v[-1,], self.xres + 2, MPI.DOUBLE], source=above)
            
    def MainLoop(self, frequency=100):
      while self.iteration < self.Max_iterations:
        self.SimulateOneTimestep()
        if not self.iteration % frequency:
          self.mesh["state/cycle"] = self.iteration
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
        extracts["e1/type"]="relay"
        extracts["e1/params/path"] = savedir + "mesh";
        extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";
        self.a.execute(action)
        self.a.close()

# Main program
#
def main():
    sim = ParallelSimulation_With_Ascent(resolution=64, iterations=500)
    sim.Initialize()
    sim.MainLoop(frequency=50)
    sim.Finalize(savedir="/scratch/snx3000/jfavre/TAscent/")

main()
