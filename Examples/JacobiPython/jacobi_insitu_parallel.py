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
      while self.iteration < self.Max_iterations:
        self.SimulateOneTimestep()

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
    sim = ParallelSimulation(resolution=64, iterations=500)
    sim.Initialize()
    sim.MainLoop(frequency=50)
    sim.Finalize()

main()
