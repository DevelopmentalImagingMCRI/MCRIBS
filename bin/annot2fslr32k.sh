#!/bin/bash

export SUBJECTS_DIR=freesurfer

rm -fr freesurfer/fsaverage freesurfer/fs_LR_32k
cp -r `dirname $0`/../lib/fs_LR_32k `dirname $0`/../lib/fsaverage freesurfer

if [ -z "$2" ]
then
    echo "Usage: $0 <subject id> <annot>"
    echo
    echo "Transforms an annot from native space to fs_LR_32k space"
    echo "Example: $0 foo aparc+DKTatlas"
    exit
fi

SUBJ=$1
ATLAS=$2

for HEMI in lh rh
do
    mri_surf2surf --srcsubject $SUBJ --trgsubject fsaverage --hemi $HEMI --sval-annot $ATLAS.annot --o $SUBJECTS_DIR/fsaverage/label/$HEMI.${SUBJ}.$ATLAS.tmp.annot --srcsurfreg sphere.reg2 --trgsurfreg sphere.reg2
    #mri_surf2surf --cortex --srcsubject fsaverage --trgsubject fs_LR_32k --hemi $HEMI --sval-annot ${SUBJ}.$ATLAS.tmp.annot --o $SUBJECTS_DIR/$1/label/$HEMI.$ATLAS.fslr32k.annot --srcsurfreg sphere.reg.fslr32k --trgsurfreg sphere
    mri_surf2surf --cortex --srcsubject fsaverage --trgsubject fs_LR_32k --hemi $HEMI --sval-annot ${SUBJ}.$ATLAS.tmp.annot --o $SUBJECTS_DIR/$1/label/$HEMI.$ATLAS.fslr32k.annot --srcsurfreg sphere.reg.fslr32k --trgsurfreg sphere
    rm -f fsaverage/label/$HEMI.${SUBJ}.$ATLAS.tmp.annot
done

