
. ~/adios2.sh

ccmake -B buildSSTConsumer -S . \
 -DCMAKE_BUILD_TYPE=Release \
 -DFides_DIR=/local/apps/fides-adios2.9-install/lib/cmake/fides \
 -DVTKm_DIR=/local/apps/Ascent/install/vtk-m-v2.0.0/lib/cmake/vtkm-2.0 \
 -DADIOS2_DIR=/local/apps/ADIOS2-v2.9.1-install/lib/cmake/adios2
 
 cmake --build buildSSTConsumer
 
# First run. Classic output to an ADIOS BP file
# make sure adios2.xml contains the string <engine type="BP4">
 
python3 heat_diffusion_insitu_serial.py 

bpls -l diffusion.bp

# Second run. In transit transfer to a data consumer

# make sure adios2.xml contains the string <engine type="SST">

rm -rf diffusion.bp
. ~/vtkm.sh
python3 heat_diffusion_insitu_serial.py &

./buildSSTConsumer/fides-sst-reader
# view images
eog diffusion*png

