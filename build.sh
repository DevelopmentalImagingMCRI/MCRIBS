#!/bin/bash

MCRIBSDIR=`pwd`

. /etc/os-release

BUILDTYPE=Release

case $ID in
	ubuntu)
		CMAKE=cmake
		MIRTKVTKDEPENDS=$MCRIBSDIR/VTK/VTK-install/lib/cmake/vtk-8.2
		USINGFLANN=ON
	;;
	centos)
		CMAKE=cmake3
		MIRTKVTKDEPENDS=$MCRIBSDIR/VTK/VTK-install/lib64/cmake/vtk-8.2
		USINGFLANN=OFF
		# boost-devel tbb-devel cmake3
	;;
esac

if [ ! -f "$MIRTKVTKDEPENDS/VTKConfig.cmake" ]
then
	# get and build VTK
	mkdir -p VTK
	#rm -fr VTK/VTK-install
	cd VTK
	git clone https://github.com/Kitware/VTK.git
	cd VTK
	git checkout tags/v8.2.0
	cd ..

	
	#LIGHTWEIGHTPYTHON="-DVTK_WRAP_PYTHON:BOOL=ON -DVTK_Group_StandAlone:BOOL=OFF -DVTK_Group_Rendering:BOOL=OFF -DModule_vtkCommonColor:BOOL=ON -DModule_vtkCommonComputationalGeometry:BOOL=ON -DModule_vtkCommonDataModel:BOOL=ON -DModule_vtkCommonExecutionModel:BOOL=ON -DModule_vtkCommonMath:BOOL=ON -DModule_vtkCommonMisc:BOOL=ON -DModule_vtkCommonSystem:BOOL=ON -DModule_vtkCommonTransforms:BOOL=ON -DModule_vtkFiltersCore:BOOL=ON -DModule_vtkFiltersExtraction:BOOL=ON -DModule_vtkFiltersGeneral:BOOL=ON -DModule_vtkFiltersGeneric:BOOL=ON -DModule_vtkFiltersGeometry:BOOL=ON -DModule_vtkFiltersPython:BOOL=ON -DModule_vtkIOCore:BOOL=ON -DModule_vtkIOGeometry:BOOL=ON -DModule_vtkIOLegacy:BOOL=ON -DModule_vtkWrappingPythonCore:BOOL=ON -DModule_vtkIOXML:BOOL=ON -DModule_vtkIOXMLPoly:BOOL=ON"

	#rm -fr VTK-python2-build
	#mkdir -p VTK-python2-build
	#cd VTK-python2-build
	# ubuntu
	#apt install -y libglu1-mesa-dev libxt-dev python-dev gcc g++
#$CMAKE -DCMAKE_CXX_FLAGS="-fPIC" -DCMAKE_C_FLAGS="-fPIC" -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/VTK/VTK-python2-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DVTK_LEGACY_SILENT=ON -DVTK_WRAP_PYTHON=ON -DVTK_PYTHON_VERSION=2 -DVTK_Group_Rendering=OFF ../VTK
	#$CMAKE -DCMAKE_CXX_FLAGS="-fPIC" -DCMAKE_C_FLAGS="-fPIC" -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/VTK/VTK-python2-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DVTK_LEGACY_SILENT=ON -DVTK_WRAP_PYTHON:BOOL=OFF -DVTK_Group_Rendering:BOOL=OFF ../VTK
	#make -j`nproc`
	#make install
	#cd ..
	
	#rm -fr VTK-build VTK-install
	#mkdir -p VTK-build
	cd VTK-build
	$CMAKE -DCMAKE_CXX_FLAGS="-fPIC" -DCMAKE_C_FLAGS="-fPIC" -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/VTK/VTK-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DVTK_LEGACY_SILENT=ON -DVTK_WRAP_PYTHON=ON -DVTK_PYTHON_VERSION=3 -DVTK_Group_Rendering:BOOL=OFF ../VTK
	make -j`nproc`
	make install
	cd $MCRIBSDIR

fi
#exit
#apt install -y zlib1g-dev libboost-dev libeigen3-dev libflann-dev
rm -fr MIRTK/MIRTK-build MIRTK/MIRTK-install
mkdir -p MIRTK/MIRTK-build
cd MIRTK/MIRTK-build
$CMAKE -DMODULE_Deformable=ON -DMODULE_Mapping=ON -DMODULE_PointSet=ON -DMODULE_Scripting=ON -DWITH_TBB=ON -DMODULE_DrawEM=ON -DWITH_VTK=ON -DDEPENDS_VTK_DIR=$MIRTKVTKDEPENDS -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/MIRTK/MIRTK-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DWITH_FLANN=$USINGFLANN ../MIRTK
make -j`nproc`
make install

# the python 2 directories can be deleted since 
rm -fr VTK/VTK-build VTK/VTK-install/include
