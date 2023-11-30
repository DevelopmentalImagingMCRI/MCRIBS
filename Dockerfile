FROM ubuntu:22.04

RUN sed --in-place --regexp-extended "s/(\/\/)(archive\.ubuntu)/\1au.\2/" /etc/apt/sources.list && \
	apt-get update && apt-get upgrade --yes

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get install -y bash \
	zlib1g-dev \
	libboost-dev \
	libglu1-mesa-dev \
	libxt-dev \
	python3-dev \
	libtbb2-dev \
	libflann-dev \
	libeigen3-dev \
	python3-contextlib2 \
	python3-imageio \
	python3-numpy \
	python3-scipy \
	python3-pandas \
	python3-numexpr \
	python3-skimage \
	python3-vtk7 \
	python3-h5py \
	cmake \
	libglvnd-dev \
	build-essential \
	gcc gcc-12 \
	g++ g++-12 \
	gfortran \
	git \
	parallel \
	libxcursor-dev \
	curl \
	file \
	xz-utils \
	dc \
	bc \
	tcsh

RUN mkdir -p /opt/MCRIBS
RUN mkdir -p /opt/MCRIBS/MIRTK

COPY MIRTK/MIRTK /opt/MCRIBS/MIRTK/MIRTK
COPY build.sh /opt/MCRIBS

COPY .git /opt/MCRIBS/.git

RUN ls -la /opt/MCRIBS
WORKDIR /opt/MCRIBS

RUN /opt/MCRIBS/build.sh

RUN mkdir -p /opt/ANTs

# build ANTs
WORKDIR /opt/ANTs
RUN git clone https://github.com/ANTsX/ANTs.git
WORKDIR /opt/ANTs/ANTs
RUN git checkout fe3a0e3

RUN mkdir -p /opt/ANTs/ANTs-build
WORKDIR /opt/ANTs/ANTs-build
RUN cmake -DCMAKE_C_COMPILER=gcc-12 -DCMAKE_CXX_COMPILER=g++-12 -DBUILD_TESTING=OFF ../ANTs
RUN make -j`nproc`

RUN mkdir -p /opt/ANTs/bin

RUN cp -r /opt/ANTs/ANTs-build/ANTS-build/Examples/* /opt/ANTs/bin
RUN cp /opt/ANTs/ANTs/Scripts/* /opt/ANTs/bin
RUN rm -fr /opt/ANTs/bin/*.cxx /opt/ANTs/bin/*.cmake /opt/ANTs/bin/*.a /opt/ANTs/bin/CMakeFiles

WORKDIR /opt

# copy in FSL
ENV FSLDIR="/opt/fsl"
RUN curl -s https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py -o fslinstaller.py
RUN python3 fslinstaller.py -d /opt/fsl -V 6.0.6.5
RUN rm fslinstaller.py

# copy in Freesurfer
COPY freesurfer-7.4.1.tar.xz /opt
RUN tar xpf freesurfer-7.4.1.tar.xz 

RUN touch /opt/entrypoint.sh &&\
echo "#!/bin/bash" >> /opt/entrypoint.sh &&\
echo "export FSLDIR=/opt/fsl" >> /opt/entrypoint.sh &&\
echo ". \$FSLDIR/etc/fslconf/fsl.sh" >> /opt/entrypoint.sh &&\
echo "export FREESURFER_HOME=/opt/freesurfer-7.4.1" >> /opt/entrypoint.sh &&\
echo ". \$FREESURFER_HOME/SetUpFreeSurfer.sh" >> /opt/entrypoint.sh &&\
echo "export MCRIBS_HOME=/opt/MCRIBS" >> /opt/entrypoint.sh &&\
echo ". \$MCRIBS_HOME/SetUpMCRIBS.sh" >> /opt/entrypoint.sh &&\
echo "export ANTSPATH=/opt/ANTs/bin" >> /opt/entrypoint.sh &&\
echo "export PATH=\$PATH:\$ANTSPATH" >> /opt/entrypoint.sh &&\
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$FREESURFER_HOME/lib/vtk" >> /opt/entrypoint.sh &&\
echo "export QT_QPA_PLATFORM_PLUGIN_PATH=/usr/lib/x86_64-linux-gnu/qt5/plugins/platforms" >> /opt/entrypoint.sh &&\
echo "export QT_PLUGIN_PATH=/usr/lib/x86_64-linux-gnu/qt5/plugins" >> /opt/entrypoint.sh &&\
echo "export FS_LICENSE=/opt/freesurfer-license.txt" >> /opt/entrypoint.sh &&\
echo "CMD=\$1" >> /opt/entrypoint.sh &&\
echo "shift;" >> /opt/entrypoint.sh &&\
echo "/opt/MCRIBS/bin/\$CMD \$@" >> /opt/entrypoint.sh &&\
chmod 777 /opt/entrypoint.sh

RUN rm -fr /opt/MCRIBS/ITK/ITK-build /opt/MCRIBS/VTK/VTK-build /opt/MCRIBS/MIRTK/MIRTK-build
RUN rm -f freesurfer-7.4.1.tar.xz
RUN rm -fr /opt/ANTs/ANTs /opt/ANTs/ANTs-build

COPY bin /opt/MCRIBS/bin
COPY lib /opt/MCRIBS/lib

COPY SetUpMCRIBS.sh /opt/MCRIBS
RUN mkdir -p /work
WORKDIR /work

ENTRYPOINT ["/opt/entrypoint.sh"]
