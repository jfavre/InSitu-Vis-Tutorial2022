##############################################################################
# A simple simulator for the heat equation in 2D
#
# Author: Jean M. Favre, Swiss National Supercomputing Center
#
# First version runs until completion and makes a matplotlib figure of the
# final scalar field
#
# run: python3 jacobi_insitu_serial.py
#
# Tested with Python 3.9.7, Wed Apr 27 12:02:12 PM CEST 2022
##############################################################################

import sys, math
import numpy as np
import matplotlib

import matplotlib.pyplot as plt
import adios2

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
        self.iteration = 0 # current iteration
        self.Max_iterations = iterations
        self.xres = resolution
        self.yres = self.xres
        self.dx = 1.0 / (self.xres + 1)

    def Initialize(self):
        """ 2 additional boundary points are added. Iterations will only touch
        the internal grid points.
        """
        self.rmesh_dims = [self.yres + 2, self.xres + 2]
        print("grid dimensions = ", self.rmesh_dims)
        self.v = np.zeros(self.rmesh_dims)
        self.vnew = np.zeros([self.yres, self.xres])
        self.set_initial_bc()
        
    def Initialize_ADIOS(self):
        self.adios = adios2.ADIOS(configFile="adios2.xml")
        self.io    = self.adios.DeclareIO("writerIO")
        self.T_id  = self.io.DefineVariable("temperature", self.v, [1, self.yres+2, self.xres+2],
                                             [0,0,0], [1, self.yres+2, self.xres+2], adios2.ConstantDims)
        self.io.DefineAttribute("Fides_Data_Model", "uniform")
        self.io.DefineAttribute("Fides_Origin", np.array([0, 0., 0.]))
        self.io.DefineAttribute("Fides_Spacing", np.array([self.dx, self.dx, self.dx]))
        self.io.DefineAttribute("Fides_Dimension_Variable", "temperature")
        self.io.DefineAttribute("Fides_Variable_List", ["temperature"])
        self.io.DefineAttribute("Fides_Variable_Associations", ["points"])

    def set_initial_bc(self):
        """ initial values set to 0 except on bottom and top wall """
        #first (bottom) row
        self.v[0,:] = [math.sin(math.pi * j * self.dx)
                       for j in range(self.rmesh_dims[1])]
        #last (top) row
        self.v[-1,:] = self.v[0,:]* math.exp(-math.pi)

    def Finalize(self):
        """plot the scalar field iso-contour lines"""
        fig, ax = plt.subplots()
        CS = ax.contour(self.v, levels=10)
        ax.clabel(CS, inline=True, fontsize=10)
        ax.set_title('Temperature iso-contours')
        plt.savefig('Temperature-iso-contours.00.png')
        #plt.show()
        
    def SimulateOneTimestep(self):
        self.iteration += 1
        # print("Simulating time step: iteration=%d" % self.iteration)

        self.vnew = 0.25 * ( self.v[2:, 1:-1]  + # north neighbor
                             self.v[0:-2, 1:-1] + # south neighbor
                             self.v[1:-1, 2:] + # east neighbor
                             self.v[1:-1, :-2]) # west neighbor
# copy vnew to the interior region of v, leaving the boundary walls untouched.
        self.v[1:-1,1:-1] = self.vnew.copy()

    def MainLoop(self):
      engine = self.io.Open("diffusion.bp", adios2.Mode.Write)
      while self.iteration < self.Max_iterations:
        self.SimulateOneTimestep()
        engine.BeginStep()
        engine.Put(self.T_id, self.v)
        engine.EndStep()
      engine.Close()

def main():
    sim = Simulation(resolution=64, iterations=500)
    sim.Initialize()
    sim.Initialize_ADIOS()
    sim.MainLoop()
    sim.Finalize()

main()
