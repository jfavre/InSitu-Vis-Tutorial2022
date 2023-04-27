
export PYTHONPATH=/local/apps/Ascent/install/conduit-v0.8.6/python-modules:/local/apps/Ascent/install/ascent-v0.9.0/python-modules
export LD_LIBRARY_PATH=/local/apps/Ascent/install/ascent-v0.9.0/lib:/local/apps/Ascent/install/vtk-m-v1.9.0/lib

rm *png
python3 double_gyre_ascent.py

