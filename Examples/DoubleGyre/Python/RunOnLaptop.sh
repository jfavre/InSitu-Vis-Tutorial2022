
export PYTHONPATH=/local/apps/Ascent/install/conduit-v0.8.8/python-modules:/local/apps/Ascent/install/ascent-develop/python-modules
export LD_LIBRARY_PATH=/local/apps/Ascent/install/ascent-develop/lib:/local/apps/Ascent/install/vtk-m-v2.0.0/lib

rm *png
python3 double_gyre_ascent.py

