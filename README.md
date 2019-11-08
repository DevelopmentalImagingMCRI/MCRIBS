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

`yum install zlib-devel boost-devel mesa-libGLU-devel libXt-devel python-devel python3-devel tbb-devel eigen3-devel python-contextlib2`

The script build.sh will checkout VTK 8.2.0 and build it. TODO: get MIRTK and python to use bindings in the VTK 8.2.0 build directory.

# Setting up

Set up the environment variable `MCRIBS_HOME` to the directory you checked out the git repo. Then source the file SetUpMCRIBS.sh in the root directory as follows:

`. $MCRIBS_HOME/SetUpMCRIBS.sh`

The main script to run is MCRIBReconAll. To get usage run MCRIBReconAll --help

## TODO

- [x] Labels unmyelinated WM are myelinated PLIC and ALIC???. The regions were correctly labelled as myelinated WM. The main WM labels are kept as "cerebral WM" despite them being unmyelinated. It makes it easier for display purposes.
- [ ] Make a MCRIB version of the ALBERTs? Working on this.
- [ ] Folding index
- [ ] Try MANTIS as a stitching thing with DrawEM.
- [ ] Make a CSH version of the setup script.
- [ ] Screenshot generators.
- [x] Investigate warning message in VTK 8.2.0 during merge-surfaces. vtkCleanPolyData causes an OOB lookup on one of the arrays. I have checked this. It seems to be a warning message that is handled correctly. I have compared vtk-8.1 vs. vtk-8.2 outputs of merge-surface in 4 subjects and they are identical. So this error can be ignored.
- [x] Fix memory leak in gcc 8 and gcc 9. evaluate-mesh, VTK issue [17722](https://gitlab.kitware.com/vtk/vtk/issues/17722), fixed in commit 95a252c04b52deb4448c02c4bd26d39d3e3688c7
