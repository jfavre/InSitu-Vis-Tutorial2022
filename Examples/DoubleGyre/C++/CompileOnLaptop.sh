
cmake -B buildAscent -S . \
      -DINSITU=Ascent \
      -DAscent_DIR=/local/apps/Ascent/install/ascent-develop/lib/cmake/ascent

cmake --build buildAscent

./buildAscent/bin/double_gyre_ascent 128 64 100

ll vort_mag.0* vel_mag.0*

=============================================================================
cmake -B buildCatalyst -S . \
      -DINSITU=Catalyst \
      -Dcatalyst_DIR=/local/apps/catalyst-dev-install/lib/cmake/catalyst-2.0

cmake --build buildCatalyst

export CATALYST_IMPLEMENTATION_PATHS=/local/apps/ParaView/dev/lib/catalyst
export VTK_SILENCE_GET_VOID_POINTER_WARNINGS=1

./buildCatalyst/bin/double_gyre_catalyst 128 64 10 ../Python/pvDoubleGyre.py
