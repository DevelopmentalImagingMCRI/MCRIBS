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
  - python-contextlib2 python3-contextlib2 python3-imageio python3-numpy python3-scipy python3-pandas python3-numexpr python3-h5py

This will install all of the dependencies for ubuntu:

`apt install zlib1g-dev libboost-dev libglu1-mesa-dev libxt-dev python-dev python3-dev libtbb-dev libflann-dev libeigen3-dev python-contextlib2 python3-contextlib2 python3-imageio python3-numpy python3-scipy python3-pandas python3-numexpr  python3-h5py`

In centos, do this:

`yum install zlib-devel boost-devel mesa-libGLU-devel libXt-devel python-devel python3-devel tbb-devel eigen3-devel python-contextlib2 python3-contextlib2 python3-imageio python3-numpy python3-scipy python3-pandas python3-numexpr python3-h5py`

The script build.sh will checkout VTK 8.2.0 and build it.

# Setting up

Assuming you have just executed the `git checkout` command:

`git checkout https://github.com/DevelopmentalImagingMCRI/MCRIBS.git`

Then depending on the shell you are using source the configuration files as follows:

- Bash: `. MCRIBS/SetUpMCRIBS.sh`
- TCSH shell: `source MCRIBS/SetUpMCRIBS.csh`

The main script to run is MCRIBReconAll. To get usage run MCRIBReconAll --help

## dHCP interoperability

### Import data

If a subject has already been put through the dHCP subject reconstruction pipeline, data can be imported with the command

`MCRIBDHCPImport <subj id> <dHCP dir>`

The `<subj id>` is the intended id of the subject after importing. The parameter `<dHCP dir>` is the "anat" directory in the derivatives directory. It should be something like `../derivatives/sub-subject1/ses-session1/anat`. The prefix `sub-subject1_ses-session1` is at the start of each file, the import tool looks for the following files in the `<dHCP dir>` directory:

* `sub-subject1_ses-session1_T2w.nii.gz`
* `sub-subject1_ses-session1_drawem_all_labels.nii.gz`
* `sub-subject1_ses-session1_brainmask_drawem.nii.gz`
* `Native/sub-subject1_ses-session1_right_corr_thickness.shape.gii`
* `Native/sub-subject1_ses-session1_right_curvature.shape.gii`
* `Native/sub-subject1_ses-session1_right_drawem.label.gii`
* `Native/sub-subject1_ses-session1_right_inflated.surf.gii`
* `Native/sub-subject1_ses-session1_right_midthickness.surf.gii`
* `Native/sub-subject1_ses-session1_right_pial.surf.gii`
* `Native/sub-subject1_ses-session1_right_roi.shape.gii`
* `Native/sub-subject1_ses-session1_right_sphere.surf.gii`
* `Native/sub-subject1_ses-session1_right_sulc.shape.gii`
* `Native/sub-subject1_ses-session1_right_thickness.shape.gii`
* `Native/sub-subject1_ses-session1_right_very_inflated.surf.gii`
* `Native/sub-subject1_ses-session1_right_white.surf.gii`
* `Native/sub-subject1_ses-session1_thickness.dscalar.nii`
* `Native/sub-subject1_ses-session1_sulc.dscalar.nii`

After running the import tool with the following indicative command:

`MCRIBDHCPImport subj ../derivatives/sub-subject1/ses-session1/anat`

, use the following command to execute the MCRIB pipeline

`MCRIBReconAll -all subj`

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

## TODO

- [x] Labels unmyelinated WM are myelinated PLIC and ALIC???. The regions were correctly labelled as myelinated WM. The main WM labels are kept as "cerebral WM" despite them being unmyelinated. It makes it easier for display purposes.
- [ ] Make a MCRIB version of the ALBERTs? Working on this.
- [ ] Folding index
- [ ] Try MANTIS as a stitching thing with DrawEM.
- [ ] Make a CSH version of the setup script.
- [ ] Screenshot generators.
- [x] Investigate warning message in VTK 8.2.0 during merge-surfaces. vtkCleanPolyData causes an OOB lookup on one of the arrays. I have checked this. It seems to be a warning message that is handled correctly. I have compared vtk-8.1 vs. vtk-8.2 outputs of merge-surface in 4 subjects and they are identical. So this error can be ignored.
- [x] Fix memory leak in gcc 8 and gcc 9. evaluate-mesh, VTK issue [17722](https://gitlab.kitware.com/vtk/vtk/issues/17722), fixed in commit [95a252c04b52deb4448c02c4bd26d39d3e3688c7](https://github.com/DevelopmentalImagingMCRI/MCRIBS/commit/95a252c04b52deb4448c02c4bd26d39d3e3688c7)
