#!/bin/bash

export SUBJECTS_DIR=freesurfer

if [ -z "$2" ]
then
    echo "Usage: $0 <subject id> <annot>"
    echo
    echo "Transforms an annot from native space to dHCP42 space"
    echo "Example: $0 foo aparc+DKTatlas"
    exit
fi

SUBJ=$1
ATLAS=$2

for HEMI in lh rh
do
    mri_surf2surf --cortex --srcsubject $SUBJ --trgsubject fsaverage --hemi $HEMI --sval-annot $ANNOT.annot --o $SUBJECTS_DIR/fsaverage/label/$HEMI.${SUBJ}.tmp.annot --srcsurfreg sphere.reg2 --trgsurfreg sphere.reg2
    mri_surf2surf --srcsubject fsaverage --trgsubject dHCP42 --hemi $HEMI --sval-annot ${SUBJ}.tmp.annot --o $SUBJECTS_DIR/$1/label/$HEMI.$ATLAS.dHCP42.annot --srcsurfreg sphere.reg.dHCP42 --trgsurfreg sphere.reg 
    rm -f fsaverage/label/$HEMI.${SUBJ}.tmp.annot
done

