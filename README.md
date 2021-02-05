# M-CRIB-S
Surface version of M-CRIB.

This software uses DrawEM and Deformable within MIRTK to perform cortical surface extraction.
A customised Freesurfer-like pipeline is provided to perform cortical parcellation on the surfaces with M-CRIB compatible labelling schemes.

The current paper is at https://www.biorxiv.org/content/10.1101/544304v1.

MIRTK was downloaded from https://github.com/BioMedIA/MIRTK. A customised version is included in this repo.

# Installation

The script build.sh will build ITK, VTK and MIRTK. Install the following dependencies prior to running:

- For VTK
  - zlib1g-dev libboost-dev libglu1-mesa-dev libxt-dev python3-dev libeigen3-dev
- For MIRTK
  - libtbb-dev libflann-dev python-contextlib2
- For Python
  - python-contextlib2 python3-contextlib2 python3-imageio python3-numpy python3-scipy python3-pandas python3-numexpr

This will install all of the dependencies for ubuntu:

`apt install zlib1g-dev libboost-dev libglu1-mesa-dev libxt-dev python-dev python3-dev libtbb-dev libflann-dev libeigen3-dev python-contextlib2 python3-contextlib2 python3-imageio python3-numpy python3-scipy python3-pandas python3-numexpr`

In centos, do this:

`yum install zlib-devel boost-devel mesa-libGLU-devel libXt-devel python-devel python3-devel tbb-devel eigen3-devel python-contextlib2 python3-contextlib2 python3-imageio python3-numpy python3-scipy python3-pandas python3-numexpr`

The script build.sh will checkout VTK 8.1.0, ITK 5.0.1 and build it.

# Setting up

Assuming you have just executed the `git checkout` command:

`git checkout https://github.com/DevelopmentalImagingMCRI/MCRIBS.git`

Then depending on the shell you are using source the configuration files as follows:

- Bash: `. MCRIBS/SetUpMCRIBS.sh`
- TCSH shell: `source MCRIBS/SetUpMCRIBS.csh`

The main script to run is MCRIBReconAll. To get usage run MCRIBReconAll --help

## Set up

Place your raw structural images, for subject `subjid` into subdirectories as follows:  

* T2 structural: `RawT2/subjid.nii.gz`
* T1 structural: `RawT1/subjid.nii.gz` (optional)

## Running

There is a wrapper script called `MCRIBReconAll` that runs the pipeline of segmentation, surface extraction. The synopsis of the command is:

`MCRIBReconAll [processing directive] [options] <subject id>`

### Processing directives and options

- --conform: Reorients T2 image to radiological orientation, axial slices. Resamples to isotropic voxels. Input <RawT2/subjid.nii.gz> Output <RawT2RadiologicalIsotropic/subjid.nii.gz>
  - Options:
  - --voxelsize VOXELSIZE: Voxel size to use for isotropic resampling. Use "volumepreserve" to preserve original voxel volume
- --tissueseg: Tissue type segmentation, depends on --conform. Input <RawT2RadiologicalIsotropic/subjid.nii.gz> Outputs <TissueSeg>
  - Options:
  - --tissuesegmethod {DrawEM}: Specify tissue segmentation method, only DrawEM is supported
  - --subjectage {28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44}: Subject age in weeks
- --surfrecon: Cortical surface extraction, depends --tissueseg. Inputs <TissueSeg>, Outputs <SurfReconDeformable/subjid>
  - Options:
  - --surfreconmethod {Deformable}: Specify cortical surface extraction method
  - --deformablejointhresh DEFORMABLEJOINTHRESH: Join threshold parameter for Deformable, default=1.
  - --deformablefastcollision: Use Deformable fast collision test
- --inflatesphere: Perform inflation, spherical mapping, curv, area, depends --surfrecon. Inputs <SurfReconDeformable/subjid>, Outputs <freesurfer/subjid>
- --surfreg: Perform surface registration to the spherical template, depends --surfrecon
- --surfvol: Perform surface volume calculations, depends --surfrecon
- --cortribbon: Perform cortical ribbon volume generation, depends --surfrecon
- --cortparc: Perform cortical parcellation, depends on --surfreg
  - Options:
  - --cortparcatlases {aparc,aparc+DKTatlas} [{aparc,aparc+DKTatlas} ...]: Parcellation scheme(s) to use
- --aparc2aseg: Transfer cortical parcellations to volume, depends on --surfreg
  - Options:
  - --aparc2asegatlases {aparc,aparc+DKTatlas} [{aparc,aparc+DKTatlas} ...] Parcellation scheme(s) to use
- --apas2aseg: Refine aseg.presurf using cortical parcellations, depends on --aparc2aseg with aparc
- --cortthickness: Compute cortical thickness
- --segstats: Perform segstats on the aseg image, depends on --apas2aseg
- --parcstats: Perform stats on the cortical surfaces, depends on --apas2aseg

The following shorthand options may be used:

- --all: executes --conform, --tissueseg, --surfrecon, --inflatesphere, --surfreg, --surfvol, --cortribbon, --cortparc, --aparc2aseg, --apas2aseg, --cortthickness, --segstats, --parcstats
- -autoreconaftersurf: executes all steps post surface extraction, i.e. --inflatesphere, --surfreg, --surfvol, --cortribbon, --cortparc, --aparc2aseg, --apas2aseg, --cortthickness, --segstats, --parcstats

Other options:
- --use-gpu, -use-gpu: Use GPU for --inflatesphere steps
- -openmp, --openmp nthreads. Use nthreads for multithreading where possible.

## dHCP interoperability

### Import data

This has recently changed, need to update the documentation here.

### Export data

MCRIB data may be exported back to dHCP format using the following command:

`MCRIBDHCPExport <subj id> <dHCP dir>`

The parameter `<subj id>` is the subject id of a subject that has been processed with MCRIBReconAll. The following files are currently converted:

 - ?h.white
 - ?h.inflated
 - ?h.pial
 - ?h.sphere
 - ?h.thickness
 - ?h.sulc
 - ?h.curv
 - ?h.aparc.annot (DK parcellation)
 - ?h.aparc+DKTatlas.annot (DKT parcellation)

