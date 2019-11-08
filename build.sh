#!/bin/bash

MCRIBSDIR=`pwd`

. /etc/os-release

BUILDTYPE=Release

case $ID in
	ubuntu)
		CMAKE=cmake
		MIRTKVTKDEPENDS=$MCRIBSDIR/VTK/VTK-install/lib/cmake/vtk-8.2
	;;
	centos)
		CMAKE=cmake3
		MIRTKVTKDEPENDS=$MCRIBSDIR/VTK/VTK-install/lib64/cmake/vtk-8.2
		# boost-devel tbb-devel cmake3
	;;
esac

ARCHFLAGS=""

while [ ! -z "$1" ]
do
	if [ "$1" == "-h" -o "$1" == '--help' ]
	then
		echo "<ENVIRONMENT VARIABLES> $0 [options]"
		echo
		echo -e "ENVIRONMENT VARIABLES"
		echo -e "\tCC\tSets the C compiler."
		echo -e "\tCXX\tSets the C++ compiler."
		echo
		echo "options"
		echo -e "\t-h or --help\tThis help."
		echo -e "\t-debug\t\tBuild VTK and MIRTK with Debug switches, Release if not given."
		echo -e "\t-archnative\tUse -march=native and -mtune=native when building VTK and MIRTK."
		exit
	fi
	if [ "$1" == "-archnative" ]
	then
		ARCHFLAGS="-DCMAKE_C_FLAGS=\"-march=native -mtune=native\" -DCMAKE_CXX_FLAGS=\"-march=native -mtune=native\""
	fi
	if [ "$1" == "-debug" ]
	then
		BUILDTYPE=Debug
	fi
	shift;
done

#if [ ! -f "$MIRTKVTKDEPENDS/VTKConfig.cmake" ]
#then
	# get and build VTK
	mkdir -p VTK
	#rm -fr VTK/VTK-install
	cd VTK
	git clone https://github.com/Kitware/VTK.git
	cd VTK
	git checkout tags/v8.2.0
	cd ../..
	PATCHFILE=lib/vtk_intersection_test.patch
	patch -p0 -N --dry-run --silent < $PATCHFILE 2>/dev/null
	if [ $? -eq 0 ];
	then
	    #apply the patch
	    patch -p0 -N < $PATCHFILE
	fi
	cd VTK

	LIGHTWEIGHTPYTHON="-DVTK_WRAP_PYTHON:BOOL=ON -DVTK_Group_StandAlone:BOOL=OFF -DVTK_Group_Rendering:BOOL=OFF -DModule_vtkCommonColor:BOOL=ON -DModule_vtkCommonComputationalGeometry:BOOL=ON -DModule_vtkCommonDataModel:BOOL=ON -DModule_vtkCommonExecutionModel:BOOL=ON -DModule_vtkCommonMath:BOOL=ON -DModule_vtkCommonMisc:BOOL=ON -DModule_vtkCommonSystem:BOOL=ON -DModule_vtkCommonTransforms:BOOL=ON -DModule_vtkFiltersCore:BOOL=ON -DModule_vtkFiltersExtraction:BOOL=ON -DModule_vtkFiltersGeneral:BOOL=ON -DModule_vtkFiltersGeneric:BOOL=ON -DModule_vtkFiltersGeometry:BOOL=ON -DModule_vtkFiltersPython:BOOL=ON -DModule_vtkIOCore:BOOL=ON -DModule_vtkIOGeometry:BOOL=ON -DModule_vtkIOLegacy:BOOL=ON -DModule_vtkWrappingPythonCore:BOOL=ON -DModule_vtkIOXML:BOOL=ON -DModule_vtkFiltersHybrid:BOOL=ON -DModule_vtkFiltersModeling:BOOL=ON -DModule_vtkImagingStencil:BOOL=ON -DModule_vtkIOPLY:BOOL=ON -DModule_vtkFiltersFlowPaths:BOOL=ON -DModule_vtkFiltersParallel:BOOL=ON"

# MIRTK vtkCommonCore vtkCommonDataModel vtkCommonExecutionModel vtkFiltersCore vtkFiltersExtraction vtkFiltersFlowPaths vtkFiltersGeneral vtkFiltersGeometry vtkFiltersHybrid vtkFiltersModeling vtkImagingCore vtkImagingStencil vtkIOGeometry vtkIOLegacy vtkIOPLY vtkIOXML

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

	rm -fr VTK-build VTK-install
#rm -fr VTK-build VTK-install
	mkdir -p VTK-build
	cd VTK-build
	#$CMAKE -DCMAKE_CXX_FLAGS="-fPIC" -DCMAKE_C_FLAGS="-fPIC" -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/VTK/VTK-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DVTK_LEGACY_SILENT=ON -DVTK_WRAP_PYTHON=ON -DVTK_PYTHON_VERSION=3 -DVTK_Group_Rendering:BOOL=OFF ../VTK
#eval $CMAKE $ARCHFLAGS -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/VTK/VTK-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DVTK_LEGACY_SILENT=ON -DVTK_WRAP_PYTHON=ON -DVTK_PYTHON_VERSION=3 -DVTK_Group_Rendering:BOOL=OFF ../VTK
	eval $CMAKE $ARCHFLAGS -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/VTK/VTK-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DVTK_LEGACY_SILENT=ON -DVTK_WRAP_PYTHON=ON -DVTK_PYTHON_VERSION=3 $LIGHTWEIGHTPYTHON ../VTK

	make -j`nproc`
	make install
	cd $MCRIBSDIR

#fi
#exit
#apt install -y zlib1g-dev libboost-dev libeigen3-dev libflann-dev

if [ ! -z "`which ccache`" ]
then
    WITHCCACHE=YES
else
    WITHCCACHE=NO
fi

rm -fr MIRTK/MIRTK-build MIRTK/MIRTK-install
mkdir -p MIRTK/MIRTK-build
cd MIRTK/MIRTK-build
eval $CMAKE $ARCHFLAGS -DMODULE_Deformable=ON -DMODULE_Mapping=ON -DMODULE_PointSet=ON -DMODULE_Scripting=ON -DWITH_TBB=ON -DMODULE_DrawEM=ON -DWITH_VTK=ON -DDEPENDS_VTK_DIR=$MIRTKVTKDEPENDS -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/MIRTK/MIRTK-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DWITH_FLANN=OFF -DWITH_CCACHE=$WITHCCACHE ../MIRTK
make -j`nproc`
make install

# the python 2 directories can be deleted since
#rm -fr VTK/VTK-build VTK/VTK-install/include MIRTK/MIRTK-build MIRTK/MIRTK-install/include
