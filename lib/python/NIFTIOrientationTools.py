#!/usr/bin/env python

import numpy
import nibabel
import os
import sys

import getopt

# apply an orientation transformation according to the ornt
# ornt is a [3, 2] array
# the format is
# [newx, flipx] 
# [newy, flipy] 
# [newz, flipz] 
# each row, the index is the new column
def applyOrntToNIIAffine(NII, ornt_transform):
    NIIAffine = NII.get_affine()
    
    # use fsl's method for 
    # make a transformation affine matrix
    transformAffine = numpy.zeros_like(NIIAffine)
    transformAffine[3, 3] = 1

    for curDim in range(3):
        newDim = int(ornt_transform[curDim, 0])
        transformAffine[curDim, newDim] = ornt_transform[curDim, 1]
        #print str(curDim) + " " + str(newDim)
        if ornt_transform[curDim, 1] < 0:
            transformAffine[curDim, 3] = (NII.shape[newDim] - 1) * NII.header.get_zooms()[newDim]

    pixDimsVector = numpy.concatenate((numpy.array(NII.header.get_zooms()), [1]))

    return numpy.matrix(NIIAffine) * numpy.diag(1.0 / pixDimsVector) * numpy.matrix(transformAffine) * numpy.diag(pixDimsVector)
    
