#!/bin/bash

if [ -z "$1" ]
then
	echo "Usage: $0 <subject id>"
	exit
fi

SUBJID=$1

#/extract_all_vtp_indir.sh DrawEM/$SUBJID/T2/${SUBJID}.nii.gz SurfReconDeformable/$SUBJID/temp

VTPExtractAll --surf-volgeom=RawT2RadiologicalIsotropic/${SUBJID}.nii.gz SurfReconDeformable/$SUBJID/temp/cerebrum-1.vtp &
VTPExtractAll --surf-volgeom=RawT2RadiologicalIsotropic/${SUBJID}.nii.gz SurfReconDeformable/$SUBJID/temp/cerebrum-lh.vtp &
VTPExtractAll --surf-volgeom=RawT2RadiologicalIsotropic/${SUBJID}.nii.gz SurfReconDeformable/$SUBJID/temp/cerebrum-rh.vtp &
wait;

MakeBadFaces SurfReconDeformable/$SUBJID/temp/cerebrum-1_tkr.surf SurfReconDeformable/$SUBJID/temp/cerebrum-1_tkr_bad.curv

freeview \
	-v RawT2RadiologicalIsotropic/${SUBJID}.nii.gz \
	-v TissueSeg/${SUBJID}_all_labels.nii.gz:colormap=lut:lut=$MCRIBS_HOME/lib/drawem_labels_FSLUT.txt:opacity=0.3 \
	-v SurfReconDeformable/$SUBJID/recon/regions.nii.gz:opacity=0.1 \
	-f SurfReconDeformable/$SUBJID/temp/cerebrum-lh_tkr.surf:edgecolor=red \
	-f SurfReconDeformable/$SUBJID/temp/cerebrum-rh_tkr.surf:edgecolor=green \
	-f SurfReconDeformable/$SUBJID/temp/cerebrum-1_tkr.surf:overlay=SurfReconDeformable/$SUBJID/temp/cerebrum-1_tkr_bad.curv:edgecolor=overlay
