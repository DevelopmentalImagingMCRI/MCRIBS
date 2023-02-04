#!/bin/bash

OUTDIR=RichardsScalper
mkdir -p $OUTDIR

for i in `cat subjects.list`
do
	if [ ! -f "$OUTDIR/${i}_brain.nii.gz" ]
	then
		echo $i
		neonateScalper -i T2NeckCroppedIsotropic/${i}.nii.gz -o $OUTDIR/${i}_brain.nii.gz -m ${OUTDIR}/${i}_brain_mask.nii.gz
	fi
done
