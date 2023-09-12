##############################################################################
# A simple simulator for the heat equation in 2D
#
# Parallel version
#
# Author: Jean M. Favre, Swiss National Supercomputing Center
#
# third version runs in parallel, splitting the domain in the vertical direction
#
# run: mpiexec -n 2 python3 jacobi_insitu_parallel.py
#
# Tested with Python 3.9.7, Wed Apr 27 12:02:12 PM CEST 2022
##############################################################################
import sys, math
import numpy as np
import adios2
from mpi4py import MPI

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
        """
        self.rmesh_dims = [self.yres + 2, self.xres + 2]
        print("dimensions = ", self.rmesh_dims)
        self.v = np.zeros(self.rmesh_dims) # includes 2 ghosts
        #self.vnew = np.zeros([self.yres, self.xres])
        self.set_initial_bc()
        
    def Initialize_ADIOS(self):
        self.adios = adios2.ADIOS(configFile="adios2.xml", comm=self.comm)
        self.io    = self.adios.DeclareIO("writerIO")
        #self.io    = self.adios.DeclareIO("InTransit-vis")
        self.T_id  = self.io.DefineVariable("temperature", self.v,
                                            [1, self.xres+2, self.xres+2], # Shape of global object
                                            [0, self.par_rank * self.yres, 0], # Where to begin writing
                                            [1, self.yres+2, self.xres+2], # Where to end writing
                                            adios2.ConstantDims)
        self.step_id  = self.io.DefineVariable("step", np.array([1], dtype=np.int32))
        self.io.DefineAttribute("Fides_Data_Model", "uniform")
        self.io.DefineAttribute("Fides_Origin", np.array([0, 0., 0.]))
        self.io.DefineAttribute("Fides_Spacing", np.array([self.dx, self.dx, self.dx]))
        self.io.DefineAttribute("Fides_Dimension_Variable", "temperature")
        self.io.DefineAttribute("Fides_Variable_List", ["temperature"])
        self.io.DefineAttribute("Fides_Variable_Associations", ["points"])

    def set_initial_bc(self):
        if self.par_size > 1:
          if self.par_rank == 0:
            self.v[0,:] = [math.sin(math.pi * j * self.dx)
                           for j in range(self.rmesh_dims[1])]
          if self.par_rank == (self.par_size - 1):
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
        
        if self.par_rank == 0:
          pass # print("Simulating time step: iteration=%d" % self.iteration)

        self.vnew = 0.25 * ( self.v[2:, 1:-1]  + # north neighbor
                             self.v[0:-2, 1:-1] + # south neighbor
                             self.v[1:-1, 2:] + # east neighbor
                             self.v[1:-1, :-2]) # west neighbor
        # copy now vnew to the interior region of v, leaving the boundary walls untouched.
        self.v[1:-1,1:-1] = self.vnew.copy()

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
      engine = self.io.Open("diffusion.bp", adios2.Mode.Write)
      while self.iteration < self.Max_iterations:
        # writing the ADIOS data before the first simulation step enables us to
        # verify the boundary conditions by writing the 0-th step
        if self.iteration % frequency == 0:
          engine.BeginStep()
          engine.Put(self.T_id, self.v)
          engine.Put(self.step_id, np.array([self.iteration], dtype=np.int32))
          engine.EndStep()
        self.SimulateOneTimestep()
      engine.Close()

class ParallelSimulation(Simulation):
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

    def Finalize(self):
        Simulation.Finalize(self)

# Main program
#
def main():
    sim = ParallelSimulation(resolution=64, iterations=10000)
    sim.Initialize()
    sim.Initialize_ADIOS()
    sim.MainLoop(frequency=500)
    sim.Finalize()

main()
