This directory contains 6 Python examples of a simple 2D heat equation solver

### First version runs without in-situ coupling, and writes the final image with matplotlib
heat_diffusion_insitu_serial.py

### Second version runs with an in-situ coupling with Ascent and writes the final image with matplotlib, the final mesh in Blueprint HDF5 format and temperature pseudo-coloring at regular intervals
heat_diffusion_insitu_Ascent.py

### Third version runs in parallel, without in-situ coupling, without I/O
heat_diffusion_insitu_parallel.py

### Fourth version runs in parallel, with an in-situ coupling with Ascent and with the same output as the 2nd version
heat_diffusion_insitu_parallel_Ascent.py

### Fifth version runs in serial, with an in-situ coupling with Catalyst
heat_diffusion_insitu_Catalyst.py

### Sixth version runs in parallel, with an in-situ coupling with Catalyst
heat_diffusion_insitu_parallel_Catalyst.py



