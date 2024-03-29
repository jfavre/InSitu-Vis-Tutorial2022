"""Module demonstrating the use of Ascent for in-situ visualization"""
##############################################################################
# A simple simulator for a 2D problem, with an in-situ coupling
# The data generation parameters for the vector field
# come from https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html
# with the Ascent library https://ascent.readthedocs.io/en/latest/#
#
# Author: Jean M. Favre, Swiss National Supercomputing Center
#
# this serial version runs until completion and saves images of the scalar field
# at regular intervals. It saves the final solution to a blueprint HDF5 file
#
# run: python3 double_gyre_catalyst.py
#
# Tested with Python 3.10.12, Mon 11 Sep 13:42:19 CEST 2023
#
##############################################################################
import math
import numpy as np
import catalyst
import catalyst_conduit as conduit
import catalyst_conduit.blueprint
import matplotlib.pyplot as plt

class Simulation:
    """
    An animated vector field generated by a math expression

    Attributes
    ----------
    resolution : int, int
        the number of grid points on the I and J axis
    iterations : int
        the maximum number of iterations (default 100)
    """
    def __init__(self, resolution=(32,16), iterations=10):
        self.iteration = 0
        self.timestep = 0.1
        self.max_iterations = iterations

        self.xres = resolution[0] # X horizontal resolution
        self.yres = resolution[1] # Y vertical   resolution
        xaxis = np.linspace(0., 2., self.xres)
        yaxis = np.linspace(0., 1., self.yres)
        self.x_coord, self.y_coord = np.meshgrid(xaxis, yaxis, indexing="xy")

        self.vel_x = np.zeros(self.x_coord.shape)
        self.vel_y = np.zeros(self.x_coord.shape)
        self.vel_z = np.zeros(self.x_coord.shape)
        self.A = 0.1 * np.pi
        self.w = 2.0 * np.pi/10.
        self.E = 0.25

    def compute_loop(self):
        """Computes and updates velocity fields"""
        while self.iteration < self.max_iterations:
            At = self.E * math.sin(self.w * self.iteration * self.timestep)
            Bt = 1.0 - 2.0 * At
            Ft = (At * self.x_coord*self.x_coord + Bt * self.x_coord) * np.pi
            fft = 2.0 * At * self.x_coord + Bt
            self.vel_x = -self.A * np.sin(Ft) * np.cos(np.pi*self.y_coord)
            self.vel_y =  self.A * np.cos(Ft) * np.sin(np.pi*self.y_coord)*fft
            self.iteration += 1

    def draw_matplotlib(self):
        """Draw with mathplotlib"""
        #plot the 'vel_x' field iso-contour lines
        fig, ax = plt.subplots()
        CS = ax.contour(self.vel_x, levels=10)
        ax.clabel(CS, inline=True, fontsize=10)
        ax.set_title('Vx iso-contours')
        plt.savefig(f'Vx-iso-contours.{self.iteration:03d}.png')
        #plot the velocity vectors sub-sampled
        fig1, ax1 = plt.subplots()
        stride = 10
        ax1.quiver(self.x_coord[::stride, ::stride], self.y_coord[::stride, ::stride],
                   self.vel_x[::stride, ::stride], self.vel_y[::stride, ::stride])
        ax1.set_title('Velocity vectors')
        plt.savefig(f'Velocity.{self.iteration:03d}.png')

class SimulationWithCatalyst(Simulation):
    """
    A specialized class to add in-situ visualization support

    Attributes
    ----------

    """
    def __init__(self, resolution=(256,128), iterations=10, pv_script="pvDoubleGyre.py"):
        Simulation.__init__(self, resolution, iterations)
        self.delta_x = 2.0 / (self.xres - 1)
        self.insitu = conduit.Node()
        self.pv_script = pv_script

    def compute_loop(self):
        while self.iteration < self.max_iterations:
            At = self.E * math.sin(self.w * self.iteration * self.timestep)
            Bt = 1.0 - 2.0 * At
            Ft = (At * self.x_coord*self.x_coord + Bt * self.x_coord) * np.pi
            fft = 2.0 * At * self.x_coord + Bt
            self.vel_x = -self.A * np.sin(Ft) * np.cos(np.pi*self.y_coord)
            self.vel_y =  self.A * np.cos(Ft) * np.sin(np.pi*self.y_coord)*fft
            self.iteration += 1
            
            exec_params = conduit.Node()
            channel = exec_params["catalyst/channels/grid"]
            channel["type"] = "mesh"
            mesh = channel["data"]
        
            mesh["coordsets/coords/type"] = "uniform"
            mesh["coordsets/coords/dims/i"] = self.xres
            mesh["coordsets/coords/dims/j"] = self.yres
            mesh["coordsets/coords/dims/k"] = 1

            mesh["topologies/mesh/type"] = "uniform"
            mesh["topologies/mesh/coordset"] = "coords"

            mesh["coordsets/coords/origin/x"] = 0.0
            mesh["coordsets/coords/origin/y"] = 0.0
            mesh["coordsets/coords/origin/z"] = 0.0
            mesh["coordsets/coords/spacing/dx"] = self.delta_x
            mesh["coordsets/coords/spacing/dy"] = self.delta_x
            mesh["coordsets/coords/spacing/dz"] = self.delta_x

            mesh["fields/vel_x/association"] = "vertex"
            mesh["fields/vel_x/topology"] = "mesh"
            mesh["fields/vel_x/values"].set_external(self.vel_x.ravel())

            mesh["fields/vel_y/association"] = "vertex"
            mesh["fields/vel_y/topology"] = "mesh"
            mesh["fields/vel_y/values"].set_external(self.vel_y.ravel())

            mesh["fields/Velocity/association"] = "vertex"
            mesh["fields/Velocity/topology"] = "mesh"
            mesh["fields/Velocity/values/u"].set_external(self.vel_x.ravel())
            mesh["fields/Velocity/values/v"].set_external(self.vel_y.ravel())
            mesh["fields/Velocity/values/w"].set_external(self.vel_z.ravel())

            # verify the mesh we created conforms to the blueprint
            verify_info = conduit.Node()
            if not conduit.blueprint.mesh.verify(mesh, verify_info):
              print("DoubleGyre Mesh Verify failed!")
            else:
              pass
              #print("DoubleGyre Mesh verify success!")
              #print(self.mesh.to_yaml())
            """Computes and updates velocity fields and process in-situ requests"""
            state = exec_params["catalyst/state"]
            state["timestep"] = self.iteration
            state["time"] = self.iteration *0.1
            catalyst.execute(exec_params)

    def initialize_catalyst(self):
        """Creates a Conduit node """
        self.insitu["catalyst/scripts/script/filename"] = self.pv_script
        self.insitu["catalyst_load/implementation"] = "paraview"

        # open Catalyst
        catalyst.initialize(self.insitu)

    def finalize_catalyst(self):
        """close"""
        catalyst.finalize(self.insitu)

#sim = Simulation()
sim = SimulationWithCatalyst(iterations=100, pv_script="pvDoubleGyre.py")
sim.initialize_catalyst()
sim.compute_loop()
sim.draw_matplotlib()
sim.finalize_catalyst()
