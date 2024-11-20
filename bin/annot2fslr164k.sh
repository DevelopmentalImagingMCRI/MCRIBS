#!/bin/bash

export SUBJECTS_DIR=freesurfer

rm -fr freesurfer/fsaverage freesurfer/fs_LR_164k
cp -r `dirname $0`/../lib/fs_LR_164k `dirname $0`/../lib/fsaverage freesurfer

if [ -z "$2" ]
then
    echo "Usage: $0 <subject id> <annot>"
    echo
    echo "Transforms an annot from native space to fs_LR_164k space"
    echo "Example: $0 foo aparc+DKTatlas"
    exit
fi

SUBJ=$1
ATLAS=$2

for HEMI in lh rh
do
    mri_surf2surf --cortex --srcsubject $SUBJ --trgsubject fsaverage --hemi $HEMI --sval-annot $ANNOT.annot --o $SUBJECTS_DIR/fsaverage/label/$HEMI.${SUBJ}.tmp.annot --srcsurfreg sphere.reg2 --trgsurfreg sphere.reg2
    mri_surf2surf --srcsubject fsaverage --trgsubject fs_LR_164k --hemi $HEMI --sval-annot ${SUBJ}.tmp.annot --o $SUBJECTS_DIR/$1/label/$HEMI.$ATLAS.fslr164k.annot --srcsurfreg sphere.reg.fslr164k --trgsurfreg sphere.reg 
    rm -f fsaverage/label/$HEMI.${SUBJ}.tmp.annot
done

