#!/usr/bin/env python

import numpy
import nibabel
import os
import sys

import getopt
import scipy.linalg

# apply an orientation transformation according to the ornt
# ornt is a [3, 2] array
# the format is
# [newx, flipx]
# [newy, flipy]
# [newz, flipz]
# each row, the index is the new column
def applyOrntToNIIAffine(NII, ornt_transform):
    NIIAffine = NII.affine

    # use fsl's method for
    # make a transformation affine matrix
    transformAffine = numpy.zeros_like(NIIAffine)
    transformAffine[3, 3] = 1
    for curDim in range(3):
        newDim = int(ornt_transform[curDim, 0])
        transformAffine[curDim, newDim] = ornt_transform[curDim, 1]
        if ornt_transform[curDim, 1] < 0:
            transformAffine[curDim, 3] = (NII.shape[curDim] - 1) * NII.header.get_zooms()[curDim]
    
    pixDimsVector = numpy.concatenate((numpy.array(NII.header.get_zooms())[:3], [1]))
    pixDimsVectorSwapped = numpy.concatenate((numpy.array(NII.header.get_zooms())[:3][numpy.int32(ornt_transform[:, 0])], [1]))

    return numpy.matrix(NIIAffine) * numpy.diag(1.0 / pixDimsVector) * numpy.matrix(transformAffine) * numpy.diag(pixDimsVectorSwapped)

if __name__ == "__main__":
    import nibabel

    inputNII = nibabel.load('orig_image.nii.gz')
    inputAXCodes = nibabel.aff2axcodes(inputNII.affine)
    inputOrnt = nibabel.orientations.axcodes2ornt(inputAXCodes)

    for refIMG in ['RLPAIS.nii.gz', 'RLSIPA.nii.gz', 'ISRLPA.nii.gz']:
        refNII = nibabel.load(refIMG)
        refAXCodes = nibabel.aff2axcodes(refNII.affine)
        refOrnt = nibabel.orientations.axcodes2ornt(refAXCodes)

        inputToRefTransformOrnt = nibabel.orientations.ornt_transform(inputOrnt, refOrnt)
        inputIMG = nibabel.orientations.apply_orientation(numpy.asanyarray(inputNII.dataobj), inputToRefTransformOrnt)

        outputAffine = applyOrntToNIIAffine(inputNII, inputToRefTransformOrnt)
        print(refIMG)
        print(refNII.affine)
        print(outputAffine)
        print(numpy.max(numpy.abs(refNII.affine - outputAffine)))

