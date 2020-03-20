#!/bin/bash
# ============================================================================
# Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Andreas Schuh
# Copyright 2013-2016 Antonios Makropoulos
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================


[ $# -ge 1 -a $# -le 3 ] || { echo "usage: $(basename "$0") <subject> <#jobs>" 1>&2; exit 1; }
subj=$1
njobs=1
if [ $# -gt 1 ];then njobs=$2;fi

run(){
  echo "$@"
  "$@" || exit 1
}

scriptdir=$(dirname "$BASH_SOURCE")
sdir=segmentations-data
echo $DRAWEMDIR
mkdir -p dofs
#  T=`tmpnam`
#  rm -f $T
for i in {01..10};do
  atlas=P${i}
  dof=dofs/$subj-$atlas-n.dof.gz
  
  if [ ! -f $dof ];then
#echo "mirtk register N4/$subj.nii.gz $DRAWEMDIR/atlases/ALBERTs/T2/$atlas.nii.gz $sdir/gm-posteriors/$subj.nii.gz $DRAWEMDIR/atlases/ALBERTs/gm-posteriors-v3/$atlas.nii.gz -parin $DRAWEMDIR/parameters/ireg-multichannel-structural.cfg  -dofout $dof" >> $T
    run mirtk register N4/$subj.nii.gz $DRAWEMDIR/atlases/MCRIBs/T2/$atlas.nii.gz $sdir/gm-posteriors/$subj.nii.gz $DRAWEMDIR/atlases/MCRIBs/gm-posteriors/$atlas.nii.gz -parin $DRAWEMDIR/parameters/ireg-multichannel-structural.cfg  -dofout $dof -threads $njobs
  fi

done
# if [ -f "$T" ]
# then
# 	parallel -j$njobs < $T
#rm -f $T
# fi

