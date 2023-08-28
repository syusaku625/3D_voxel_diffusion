#!/bin/sh
mkdir build
cd build
cmake -D TP_DIR=/mnt/d/lib/TextParser \
-D CMAKE_INSTALL_PREFIX=/mnt/d/work/O17/3D_voxel_diffusion/bin \
-D CMAKE_EXPORT_COMPILE_COMMANDS=ON \
..

make -j8 && make install