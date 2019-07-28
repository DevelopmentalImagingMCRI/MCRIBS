#!/bin/bash

# get and build VTK
mkdir -p VTK
rm -fr VTK/VTK-build
mkdir -p VTK/VTK-build
rm -fr VTK/VTK-install
cd VTK
#git clone https://github.com/Kitware/VTK.git
cd VTK
git checkout tags/v8.2.0
cd ../VTK-build

# ubuntu 
#apt install -y libglu1-mesa-dev libxt-dev python-dev gcc g++
cmake -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=`pwd`/VTK/VTK-install -DCMAKE_BUILD_TYPE=Release -DVTK_LEGACY_SILENT=ON -DVTK_WRAP_PYTHON=ON ../VTK
