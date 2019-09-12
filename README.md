# M-CRIB-S
Surface version of M-CRIB.

This software uses DrawEM and Deformable within MIRTK to perform cortical surface extraction.
A customised Freesurfer-like pipeline is provided to perform cortical parcellation on the surfaces with M-CRIB compatible labelling schemes.

The current paper is at https://www.biorxiv.org/content/10.1101/544304v1.

MIRTK was downloaded from https://github.com/BioMedIA/MIRTK. A customised version is included in this repo.

# Installation

The script build.sh will build VTK and MIRTK. Install the following dependencies prior to running:

- For VTK
  - zlib1g-dev libboost-dev libglu1-mesa-dev libxt-dev python-dev libeigen3-dev
- For MIRTK
  - libtbb-dev libflann-dev
- For python
  - python-vtk python3-vtk
  
This will install all of the dependencies

`apt install zlib1g-dev libboost-dev libglu1-mesa-dev libxt-dev python-dev libtbb-dev libflann-dev libeigen3-dev python-vtk python3-vtk`

The script build.sh will checkout VTK 8.2.0 and build it.

# Setting up

Set up the environment variable `MCRIBS_HOME` to the directory you checked out the git repo. Then source the file SetUpMCRIBS.sh in the root directory as follows:

`. $MCRIBS_HOME/SetUpMCRIBS.sh`

The main script to run is MCRIBReconAll. To get usage run MCRIBReconAll --help
