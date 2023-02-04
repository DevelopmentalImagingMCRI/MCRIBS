#!/bin/bash

export SUBJECTS_DIR=`pwd`/freesurfer
export TISSUESEGDIR=TissueSegMCRIBS

if [ -z "$1" ]
then
	echo "Usage $0 <subjid> [options]"
	exit
fi

SUBJID=$1

# make the MADs image, ala DrawEM

nn=5
mkdir -p $TISSUESEGDIR/$SUBJID/MADs
calculate-gradients $TISSUESEGDIR/$SUBJID/${SUBJID}_t2w_restore.nii.gz $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_t2w_restore_grad.nii.gz 0 
calculate-filtering $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_t2w_restore_grad.nii.gz -kernel $nn -median $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz
#cp $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur1.nii.gz
calculate-element-wise $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_t2w_restore_grad.nii.gz -sub $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz -abs -out $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz 
calculate-filtering  $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz -kernel $nn -median $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz  
calculate-element-wise $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_t2w_restore_grad.nii.gz -div-with-zero $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz -div 1.4826 -sq -mul 0.5 -add 1 -log -out $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz 
calculate-element-wise $TISSUESEGDIR/$SUBJID/${SUBJID}_t2w_restore.nii.gz -div-with-zero $TISSUESEGDIR/$SUBJID/${SUBJID}_t2w_restore.nii.gz -mul $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz -add 1 -out $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz 
calculate-element-wise $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz -mul 0 -add 1 -div-with-zero $TISSUESEGDIR/$SUBJID/MADs/${SUBJID}_cur.nii.gz -out $TISSUESEGDIR/$SUBJID/${SUBJID}_t2w_restore_mad.nii.gz
rm -fr $TISSUESEGDIR/$SUBJID/MADs


# make gm posterior image, add the posteriors for all GM and subcortical GM labels
GMPOSTERIORIMAGES=""
NUMGMPOSTERIORIMAGES=0
GMMATCH=
for j in `seq 1000 1035`
do
    if [ -f "$TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_${j}.nii.gz" ]
    then
	    GMPOSTERIORIMAGES="$GMPOSTERIORIMAGES $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_${j}.nii.gz"
        NUMGMPOSTERIORIMAGES=`expr $NUMGMPOSTERIORIMAGES + 1`
    fi
    if [ -f "$TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_`expr ${j} + 1000`.nii.gz" ]
    then
	    GMPOSTERIORIMAGES="$GMPOSTERIORIMAGES $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_`expr ${j} + 1000`.nii.gz"
        NUMGMPOSTERIORIMAGES=`expr $NUMGMPOSTERIORIMAGES + 1`
    fi
done
for j in 9 11 12 13 28 48 50 51 52 60 53 54 17 18
do
    jstr=$(printf "%04d" $j)
    if [ -f "$TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_${jstr}.nii.gz" ]
    then
        GMPOSTERIORIMAGES="$GMPOSTERIORIMAGES $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_${jstr}.nii.gz"
        NUMGMPOSTERIORIMAGES=`expr $NUMGMPOSTERIORIMAGES + 1`
    fi
done

mri_binarize --i $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited.nii.gz --o $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_subcort_gm.nii.gz --match 53 54 17 18 9 11 12 13 28 48 50 51 52 60  --dilate 2

fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_t2w_restore_mad.nii.gz -mul $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_subcort_gm.nii.gz $TISSUESEGDIR/$SUBJID/${SUBJID}_t2w_restore_subspace.nii.gz

#fslmerge -a $TISSUESEGDIR/$SUBJID/all_gm_posteriors.nii.gz $GMPOSTERIORIMAGES
AverageImages 3 $TISSUESEGDIR/$SUBJID/all_gm_posteriors.nii.gz 0 $GMPOSTERIORIMAGES
fslmaths $TISSUESEGDIR/$SUBJID/all_gm_posteriors.nii.gz -mul $NUMGMPOSTERIORIMAGES $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_gm
rm -f $TISSUESEGDIR/$SUBJID/all_gm_posteriors.nii.gz

mri_binarize --i $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited.nii.gz --o $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_cerebellum.nii.gz --match 91 93 170 53 54 17 18 9 11 12 13 28 48 50 51 52 60 --dilate 5 --erode 5

#mri_binarize --i $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited.nii.gz --o $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_cerebellum.nii.gz --match 91 93 --dilate 3 --erode 3

fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_brain_mask.nii.gz -sub $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_cerebellum.nii.gz $TISSUESEGDIR/$SUBJID/${SUBJID}_drawem_mask.nii.gz

# make wm
fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_0002 -add $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_0041 -add $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_0170 $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_wm

# make csf
fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_0024 -add $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_0004 -add $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_0043 -add $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_0014 $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_csf
# make background
fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_0165 -add $TISSUESEGDIR/$SUBJID/${SUBJID}_labelposterior_0258 $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_background
#mri_binarize --i 

mri_binarize --i $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited.nii.gz --o $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_background.nii.gz --match 165 258 --dilate 2

fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_background.nii.gz -binv -mas $TISSUESEGDIR/$SUBJID/${SUBJID}_segmentation_gm.nii.gz -add $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_gm $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_gm
rm -f $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_background.nii.gz

#fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_cerebellum.nii.gz -div 2 -add $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_wm $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_wm
fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_labelfusion_dkt_edited_cerebellum.nii.gz -div 2 -add $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_gm $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_gm

fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_background \
    -add $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_wm \
    -add $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_gm \
    -add $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_csf $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_sum

fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_sum -binv -add $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_sum $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_sum
fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_background -div $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_sum $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_background
fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_wm -div $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_sum $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_wm
fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_gm -div $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_sum $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_gm
fslmaths $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_csf -div $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_sum $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_csf


#rm -f $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_sum.nii.gz
#exit


# mri_binarize --i $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}.nii.gz --o $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_gm.nii.gz $GMMATCH
# mri_binarize --i $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}.nii.gz --o $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_wm.nii.gz --match 2 41
# mri_binarize --i $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}.nii.gz --o $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_background.nii.gz --match 165 258
# mri_binarize --i $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}.nii.gz --o $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_subcort_gm.nii.gz --match 9 11 12 13 28 48 50 51 52 60
# mri_binarize --i $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}.nii.gz --o $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_latvent.nii.gz --match 4 43
# mri_binarize --i $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}.nii.gz --o $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_csf.nii.gz --match 24 14 4 43
# #mri_binarize --i $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}.nii.gz --o $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_cerebellum.nii.gz --match 91 93 --dilate 3 --erode 3

# fslmaths $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_gm -add $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_subcort_gm -Tmean $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_gm_prob
# fslmaths $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_wm -Tmean $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_wm_prob
# fslmaths $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_csf -add $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_latvent -Tmean $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_csf_prob
# fslmaths $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_background -Tmean $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_background_prob
# #f#slmaths $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_latvent -Tmean $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_latvent_prob
# #fslmaths $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_cerebellum -Tmean $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_cerebellum_prob

# fslmaths $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_gm_prob -mul 0 -add 1 $TISSUESEGDIR/$SUBJID/all_dkt_to_${SUBJID}_gmm

draw-em $TISSUESEGDIR/$SUBJID/${SUBJID}_t2w_restore.nii.gz 4  \
    $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_background.nii.gz \
    $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_gm.nii.gz \
    $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_wm.nii.gz \
    $TISSUESEGDIR/$SUBJID/${SUBJID}_posterior_csf.nii.gz \
    draw-em.nii.gz -padding 0 -mrf connectivities.mrf -tissues 1 0 1 1 1 2 1 3 -hui -postpenalty $TISSUESEGDIR/$SUBJID/${SUBJID}_t2w_restore_subspace.nii.gz -mask $TISSUESEGDIR/$SUBJID/${SUBJID}_drawem_mask.nii.gz