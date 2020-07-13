#!/bin/bash

MCRIBSDIR=`pwd`

. /etc/os-release

BUILDTYPE=Release

case $ID in
	ubuntu)
		CMAKE=cmake
		MIRTKVTKDEPENDS=$MCRIBSDIR/VTK/VTK-install/lib/cmake/vtk-8.1
	;;
	centos)
		CMAKE=cmake3
		MIRTKVTKDEPENDS=$MCRIBSDIR/VTK/VTK-install/lib64/cmake/vtk-8.1
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


mkdir -p ITK
cd ITK
git clone https://github.com/InsightSoftwareConsortium/ITK.git
cd ITK
git checkout tags/v5.1.0
cd ..
mkdir -p ITK-build
cd ITK-build
pwd
eval $CMAKE $ARCHFLAGS -DBUILD_EXAMPLES:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/ITK/ITK-install ../ITK
make -j`nproc`
make install
cd $MCRIBSDIR

 # get and build VTK
mkdir -p VTK
#rm -fr VTK/VTK-install VTK/VTK-build
cd VTK
git clone https://github.com/Kitware/VTK.git
cd VTK
git checkout tags/v8.1.2
cd ../..
PATCHFILE=lib/vtk-8.1_intersection_test.patch
patch -p0 -N --dry-run --silent < $PATCHFILE 2>/dev/null
if [ $? -eq 0 ];
then
    #apply the patch
    patch -p0 -N < $PATCHFILE
fi

PYTHON3VERSION=`python3 --version | awk '{ print $2; }'`
PYTHON3MAJOR=`echo $PYTHON3VERSION | cut -f1 -d.`
PYTHON3MINOR=`echo $PYTHON3VERSION | cut -f2 -d.`

echo $PYTHON3VERSION $PYTHON3MAJOR $PYTHON3MINOR

if [ "$PYTHON3MINOR" -ge "8" ]
then
	PATCHFILE=lib/0001-Compatibility-for-Python-3.8.patch
	patch -p0 -N --dry-run --silent < $PATCHFILE 2>/dev/null
	if [ $? -eq 0 ];
	then
	    #apply the patch
	    patch -p0 -N < $PATCHFILE
	fi
fi

cd VTK

LIGHTWEIGHTPYTHON="-DVTK_WRAP_PYTHON:BOOL=ON -DVTK_Group_StandAlone:BOOL=OFF -DVTK_Group_Rendering:BOOL=OFF -DModule_vtkCommonColor:BOOL=ON -DModule_vtkCommonComputationalGeometry:BOOL=ON -DModule_vtkCommonDataModel:BOOL=ON -DModule_vtkCommonExecutionModel:BOOL=ON -DModule_vtkCommonMath:BOOL=ON -DModule_vtkCommonMisc:BOOL=ON -DModule_vtkCommonSystem:BOOL=ON -DModule_vtkCommonTransforms:BOOL=ON -DModule_vtkFiltersCore:BOOL=ON -DModule_vtkFiltersExtraction:BOOL=ON -DModule_vtkFiltersGeneral:BOOL=ON -DModule_vtkFiltersGeneric:BOOL=ON -DModule_vtkFiltersGeometry:BOOL=ON -DModule_vtkFiltersPython:BOOL=ON -DModule_vtkIOCore:BOOL=ON -DModule_vtkIOGeometry:BOOL=ON -DModule_vtkIOLegacy:BOOL=ON -DModule_vtkWrappingPythonCore:BOOL=ON -DModule_vtkIOXML:BOOL=ON -DModule_vtkFiltersHybrid:BOOL=ON -DModule_vtkFiltersModeling:BOOL=ON -DModule_vtkImagingStencil:BOOL=ON -DModule_vtkIOPLY:BOOL=ON -DModule_vtkFiltersFlowPaths:BOOL=ON -DModule_vtkFiltersParallel:BOOL=ON"

rm -fr VTK-build VTK-install
mkdir -p VTK-build
cd VTK-build
eval $CMAKE $ARCHFLAGS -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/VTK/VTK-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DVTK_LEGACY_SILENT=ON -DVTK_WRAP_PYTHON=ON -DVTK_PYTHON_VERSION=3 $LIGHTWEIGHTPYTHON ../VTK

make -j`nproc`
make install
cd $MCRIBSDIR

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
eval $CMAKE $ARCHFLAGS -DMODULE_Deformable=ON -DMODULE_Mapping=ON -DMODULE_PointSet=ON -DMODULE_Scripting=ON -DWITH_TBB=ON -DMODULE_DrawEM=ON -DWITH_VTK=ON -DITK_DIR=$MCRIBSDIR/ITK/ITK-install/lib/cmake/ITK-5.1 -DDEPENDS_VTK_DIR=$MIRTKVTKDEPENDS -DCMAKE_INSTALL_PREFIX=$MCRIBSDIR/MIRTK/MIRTK-install -DCMAKE_BUILD_TYPE=$BUILDTYPE -DWITH_FLANN=OFF -DWITH_CCACHE=$WITHCCACHE ../MIRTK
make -j`nproc`
make install

# the python 2 directories can be deleted since
#rm -fr VTK/VTK-build VTK/VTK-install/include MIRTK/MIRTK-build MIRTK/MIRTK-install/include
