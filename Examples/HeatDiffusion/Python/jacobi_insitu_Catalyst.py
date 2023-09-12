##############################################################################
# A simple simulator for the heat equation in 2D, with an in-situ coupling
# with the Catalyst library https://catalyst-in-situ.readthedocs.io/en/latest/index.html
#
# Author: Jean M. Favre, Swiss National Supercomputing Center
#
# this serial version runs until completion and saves images of the scalar field
# at regular intervals. It saves the final solution to a blueprint HDF5 file
#
# run: python3 jacobi_insitu_Catalyst.py
#
# Tested with Python 3.10.12, Mon 11 Sep 13:42:19 CEST 2023
#
#
##############################################################################
import math
import glob
import numpy as np
import matplotlib.pyplot as plt
import catalyst
import catalyst_conduit as conduit
import catalyst_conduit.blueprint


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

    def initialize(self):
        """ 2 additional boundary points are added. Iterations will only touch
        the internal grid points.
        """
        self.rmesh_dims = [self.yres + 2, self.xres + 2]
        #print("grid dimensions = ", self.rmesh_dims)
        self.v = np.zeros(self.rmesh_dims)
        self.vnew = np.zeros([self.yres, self.xres])
        self.set_initial_bc()

    def set_initial_bc(self):
        """ initial values set to 0 except on bottom and top wall """
        #first (bottom) row
        self.v[0,:] = [math.sin(math.pi * j * self.dx)
                       for j in range(self.rmesh_dims[1])]
        #last (top) row
        self.v[-1,:] = self.v[0,:]* math.exp(-math.pi)

    def finalize(self):
        """plot the scalar field iso-contour lines"""
        fig, ax = plt.subplots()
        CS = ax.contour(self.v, levels=10)
        ax.clabel(CS, inline=True, fontsize=10)
        ax.set_title('Temperature iso-contours')
        plt.savefig(f'Temperature-iso-contours.{self.iteration:04d}.png')
        #plt.show()

    def simulate_one_timestep(self):
        self.iteration += 1
        # print("Simulating time step: iteration=%d" % self.iteration)

        self.vnew = 0.25 * ( self.v[2:, 1:-1]  + # north neighbor
                             self.v[0:-2, 1:-1] + # south neighbor
                             self.v[1:-1, 2:] + # east neighbor
                             self.v[1:-1, :-2]) # west neighbor
        # copy vnew to the interior region of v, leaving the boundary walls untouched.
        self.v[1:-1,1:-1] = self.vnew.copy()

    def main_loop(self):
        while self.iteration < self.Max_iterations:
            self.simulate_one_timestep()

# we now define a sub-class of Simulation to add a Catalyst in-situ coupling

class Simulation_With_Catalyst(Simulation):
    def __init__(self, resolution=64, iterations=100, meshtype="uniform", pv_script="pvDoubleGyre.py"):
        Simulation.__init__(self, resolution, iterations)
        self.MeshType = meshtype
        if meshtype == "rectilinear":
            self.xc = np.linspace(0, 1, self.xres + 2)
            self.yc = np.linspace(0, 1, self.yres + 2)
        else:
            if meshtype in ('structured', 'unstructured'):
                self.xc, self.yc = np.meshgrid(np.linspace(0, 1, self.xres + 2),
                                           np.linspace(0, 1, self.yres + 2),
                                           indexing='xy')
        if meshtype == "unstructured":
            self.conn = np.zeros(((self.xres + 1) * (self.yres + 1) * 4), dtype=np.int32)
            i=0
            for iy in range(self.yres+1):
                for ix in range(self.xres+1):
                    self.conn[4*i+0] = ix + iy*(self.xres + 2)
                    self.conn[4*i+1] = ix + (iy+1)*(self.xres + 2)
                    self.conn[4*i+2] = ix + (iy+1)*(self.xres + 2)+ 1
                    self.conn[4*i+3] = ix + iy*(self.xres + 2) + 1
                    i += 1
        self.insitu = conduit.Node()
        self.pv_script = pv_script
        
    # Add Catalyst mesh definition
    def initialize(self):
        Simulation.initialize(self)
        self.initialize_catalyst()

    def main_loop(self, frequency=100):
      while self.iteration < self.Max_iterations:
        self.simulate_one_timestep()

        exec_params = conduit.Node()
        channel = exec_params["catalyst/channels/grid"]
        channel["type"] = "mesh"
        mesh = channel["data"]

        # create the coordinate set
        if self.MeshType in ('structured', 'unstructured'):
            mesh["coordsets/coords/type"] = "explicit"
        else:
            mesh["coordsets/coords/type"] = self.MeshType
            mesh["coordsets/coords/dims/i"] = self.xres + 2
            mesh["coordsets/coords/dims/j"] = self.yres + 2

        mesh["topologies/mesh/type"] = self.MeshType
        mesh["topologies/mesh/coordset"] = "coords"

        if self.MeshType == "uniform":
          # add origin and spacing to the coordset (optional)
            mesh["coordsets/coords/origin/x"] = 0.0
            mesh["coordsets/coords/origin/y"] = 0.0
            mesh["coordsets/coords/spacing/dx"] = self.dx
            mesh["coordsets/coords/spacing/dy"] = self.dx
        else:
            mesh["coordsets/coords/values/x"].set_external(self.xc.ravel())
            mesh["coordsets/coords/values/y"].set_external(self.yc.ravel())

        if self.MeshType == "structured":
            mesh["topologies/mesh/elements/dims/i"] = np.int32(self.xres + 1)
            mesh["topologies/mesh/elements/dims/j"] = np.int32(self.yres + 1)

        if self.MeshType == "unstructured":
            mesh["topologies/mesh/elements/shape"] = "quad"
            mesh["topologies/mesh/elements/connectivity"].set_external(self.conn)
        # create a vertex associated field called "temperature"
        mesh["fields/temperature/association"] = "vertex"
        mesh["fields/temperature/topology"] = "mesh"
        # set_external does not handle multidimensional numpy arrays or
        # multidimensional complex strided views into numpy arrays.
        # Views that are effectively 1D-strided are supported.
        mesh["fields/temperature/values"].set_external(self.v.ravel())

        # make sure the mesh we created conforms to the blueprint
        verify_info = conduit.Node()
        if not conduit.blueprint.mesh.verify(mesh, verify_info):
            print("Jacobi Mesh Verify failed!")
          #print(verify_info.to_yaml())
        else:
            pass

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
        
def main():
    #sim = Simulation(resolution=64, iterations=500)
    # choices are meshtype="uniform", "rectilinear", "structured", "unstructured"
    sim = Simulation_With_Catalyst(meshtype="uniform", iterations=5000, pv_script="../JacobiC++/catalyst_state.py")
    sim.initialize()
    sim.main_loop()
    sim.finalize_catalyst()
 
main()

# list all images which have been rendered to disk
image_files = glob.glob("datasets/*png")
image_files.sort()
print(image_files)


