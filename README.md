# M-CRIB-S
Surface version of M-CRIB.

This software uses DrawEM and Deformable within MIRTK to perform cortical surface extraction.
A customised Freesurfer-like pipeline is provided to perform cortical parcellation on the surfaces with M-CRIB compatible labelling schemes.

The current paper is at https://www.biorxiv.org/content/10.1101/544304v1.

MIRTK was downloaded from https://github.com/BioMedIA/MIRTK. A customised version is included in this repo.

# Installation

The script build.sh will build VTK and MIRTK. Install the following dependencies prior to running:

- For VTK
  - zlib1g-dev libboost-dev libglu1-mesa-dev libxt-dev python3-dev libeigen3-dev
- For MIRTK
  - libtbb-dev libflann-dev python-contextlib2
- For Python
  - python-contextlib2

This will install all of the dependencies for ubuntu:

`apt install zlib1g-dev libboost-dev libglu1-mesa-dev libxt-dev python-dev python3-dev libtbb-dev libflann-dev libeigen3-dev python-contextlib2`

In centos, do this:

`yum install zlib1g-devel boost-devel libglu1-mesa-devel libxt-devel python-devel python3-devel tbb-devel eigen3-devel python-contextlib2`

The script build.sh will checkout VTK 8.2.0 and build it. TODO: get MIRTK and python to use bindings in the VTK 8.2.0 build directory.

# Setting up

Set up the environment variable `MCRIBS_HOME` to the directory you checked out the git repo. Then source the file SetUpMCRIBS.sh in the root directory as follows:

`. $MCRIBS_HOME/SetUpMCRIBS.sh`

The main script to run is MCRIBReconAll. To get usage run MCRIBReconAll --help
