
pushd /local/apps/catalyst-dev-install/lib/catalyst/
ln -s /local/apps/ParaView/dev/lib/catalyst/libcatalyst-paraview.so .

export PYTHONPATH=/local/apps/catalyst-dev-install/lib/python3.10/site-packages:$PYTHONPATH
export CATALYST_IMPLEMENTATION_NAME=paraview
export CATALYST_IMPLEMENTATION_PATHS=""
export LD_LIBRARY_PATH=/local/apps/catalyst-dev-install/lib/catalyst:/local/apps/ParaView/dev/lib:$LD_LIBRARY_PATH

