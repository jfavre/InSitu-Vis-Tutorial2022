
export PYTHONPATH=/local/apps/Ascent/install/conduit-v0.8.7/python-modules:/local/apps/Ascent/install/ascent-v0.9.1/python-modules:$PYTHONPATH
export LD_LIBRARY_PATH=/local/apps/Ascent/install/ascent-v0.9.1/lib:/local/apps/Ascent/install/vtk-m-v1.9.0/lib

rm temperature-ser.0*png Temperature-iso-contours.01.png 
python3 jacobi_insitu_Ascent.py

