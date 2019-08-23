#!/bin/bash

MCRIBSDIR=`pwd`

if [ ! -f "$MCRIBSDIR/VTK/VTK-install/lib/cmake/vtk-8.2/VTKConfig.cmake" ]
then
	# get and build VTK
	mkdir -p VTK
	rm -fr VTK/VTK-build
	mkdir -p VTK/VTK-build
	rm -fr VTK/VTK-install
	cd VTK
	git clone https://github.com/Kitware/VTK.git
	cd VTK
	git checkout tags/v8.2.0
	cd ../VTK-build
	
	# ubuntu 
	#apt install -y libglu1-mesa-dev libxt-dev python-dev gcc g++
	cmake -DCMAKE_CXX_FLAGS="-fPIC" -DCMAKE_C_FLAGS="-fPIC" -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/VTK/VTK-install -DCMAKE_BUILD_TYPE=Release -DVTK_LEGACY_SILENT=ON -DVTK_WRAP_PYTHON=ON ../VTK
	make -j`nproc`
	make install
	cd $MCRIBSDIR
fi
#apt install -y zlib1g-dev libboost-dev libeigen3-dev libflann-dev

rm -fr MIRTK/MIRTK-build
mkdir -p MIRTK/MIRTK-build
cd MIRTK/MIRTK-build
cmake -DMODULE_Deformable=ON -DMODULE_Mapping=ON -DMODULE_PointSet=ON -DMODULE_Scripting=ON -DWITH_TBB=ON -DMODULE_DrawEM=ON -DWITH_VTK=ON -DDEPENDS_VTK_DIR=$MCRIBSDIR/VTK/VTK-install/lib/cmake/vtk-8.2 -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/MIRTK/MIRTK-install -DCMAKE_BUILD_TYPE=Release -DWITH_FLANN=ON ../MIRTK
make -j`nproc`
make install
