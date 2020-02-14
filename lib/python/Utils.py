import numpy

import imageio

import os

import scipy
import scipy.misc

import errno

import math
import numexpr

# for the label image V, return a dict() whose keys are IDX and for each key, return a tuple of indices for

def labelIDXList(V, IDX):

    In1D = numpy.in1d(V, IDX)

    In1DIDX = numpy.where(In1D)[0]

    # get label voxels that are in1D
    In1DV = V[In1D.reshape(V.shape)]

    SortedIn1DVIDX = numpy.argsort(In1DV)
    SortedIn1DV = In1DV[SortedIn1DVIDX]

    # use the histogram method to get counts
    LabelCounts = numpy.bincount(SortedIn1DV)

    outIDX = dict.fromkeys(IDX.tolist())

    curIDX = 0

    for curLabelIDX in range(LabelCounts.size):
        if LabelCounts[curLabelIDX] > 0:
                curLabel = SortedIn1DV[curIDX]
                outIDX[curLabel] = In1DIDX[SortedIn1DVIDX[curIDX:curIDX + LabelCounts[curLabelIDX]]]
                if V.ndim > 1:
                        outIDX[curLabel] = numpy.unravel_index(outIDX[curLabel], V.shape)
                curIDX += LabelCounts[curLabelIDX]
    return outIDX

V = numpy.tile(numpy.arange(5), [4, 1])

IDX = numpy.array([1, 2, 3, 6])

T = labelIDXList(V, IDX)
# regionProps
# for the label image returns a dictionary with
# 'area' (vector): the number of pixels for each label
# 'mask' (list): a binary mask for each label
# 'pixelList' (list): numpy.where() for each label
# 'centroid' (list): mean of pixel coordinates for each label
def regionProps(labelImage, properties):

    numLabels = numpy.max(labelImage)
    regionProps = dict()

    pixelLists = list()

    if numLabels >= 1:
        for curLabel in range(numLabels):
            pixelLists.append(numpy.where(labelImage == (curLabel + 1)))

    for z in range(len(properties)):
        if properties[z].lower() == 'area':
            if numLabels >= 1:
                #regionProps['area'] = scipy.ndimage.measurements.labeled_comprehension(labelImage, labelImage, numpy.arange(1, numLabels + 1), numpy.size, numpy.uint32, 0)
                regionProps['area'] = numpy.zeros((numLabels))
                #print numpy.size(pixelLists[curLabel][0])
                #print regionProps['area'].shape
                for curLabel in range(numLabels):
                    #print numpy.size(pixelLists[curLabel][0])
                    regionProps['area'][curLabel] = numpy.size(pixelLists[curLabel][0])
            else:
                regionProps['area'] = None
        elif properties[z].lower() == 'pixellist':
            regionProps['pixelList'] = pixelLists[:]
            #if numLabels >= 1:
            #    for curLabel in range(numLabels):
            #        regionProps['pixelList'].append(numpy.where(labelImage == (curLabel + 1)))
        elif properties[z].lower() == 'mask':
            regionProps['mask'] = list()
            if numLabels >= 1:
                for curLabel in range(numLabels):
                    T = numpy.zeros(labelImage.shape, dtype = numpy.bool)
                    T[pixelLists[curLabel]] = 1
                    regionProps['mask'].append(numpy.array(T))
        elif properties[z].lower() == 'centroid':
            if numLabels >= 1:
                regionProps['centroid'] = list()
                for curLabel in range(numLabels):
                    T = numpy.zeros(len(pixelLists[curLabel]))
                    for curDim in range(len(pixelLists[curLabel])):
                        T[curDim] = numpy.mean(pixelLists[curLabel][curDim])
                    regionProps['centroid'].append(T)
    return regionProps


def largestComponent(CCSeg):
    CCSegLabels, numCCSegLabels = scipy.ndimage.measurements.label(CCSeg, structure = numpy.ones([3, 3]))
    if numCCSegLabels > 1:
        CCSegRegionProps = regionProps(CCSegLabels, ['Area'])
        M = numpy.argmax(CCSegRegionProps['area'])
        newCCSeg = (CCSegLabels == (M + 1))
    else:
        newCCSeg = numpy.array(CCSeg)
    return newCCSeg

def mostCentralComponent(CCSeg):
    CCSegLabels, numCCSegLabels = scipy.ndimage.measurements.label(CCSeg, structure = numpy.ones([3, 3]))
    if numCCSegLabels > 1:
        CCSegRegionProps = regionProps(CCSegLabels, ['centroid'])
        allCentroids = numpy.stack(CCSegRegionProps['centroid'])
    #    print allCentroids
        XC = allCentroids - (numpy.atleast_2d(CCSeg.shape) / 2.0)
    #    print numpy.atleast_2d(CCSeg.shape) / 2.0
    #    print XC
        D = numpy.sqrt(numpy.sum(XC * XC, axis = 1))
        I = numpy.argmin(D)
        newCCSeg = (CCSegLabels == (I + 1))
    else:
        newCCSeg = numpy.array(CCSeg)
    return newCCSeg

def ismember(a, b):
    tf = numpy.in1d(a,b) # for newer versions of numpy
    #f = numpy.array([i in b for i in a])
    u = numpy.unique(a[tf])
    index = numpy.array([(numpy.where(b == i))[0][-1] if t else 0 for i,t in zip(a,tf)])
    return (tf, index)

#function [F, HalfMaximum] = gaussian_fwhm2d(SIGMA)

# imitates the imglob utility in fsl
# takes a image filename and checks to see if the corresponding nifti file exists, returns the prefix

supportedExtensions = ['nii', 'nii.gz']

def findNIFTIFromPrefix(filePrefix):
    for curExtension in supportedExtensions:
        if os.path.isfile(filePrefix + '.' + curExtension):
            return filePrefix + '.' + curExtension
    return None

def imglob(fileName):

    if isinstance(fileName, list):
        T = list()
        for curFileName in fileName:
            T.append(imglob(curFileName))
        return T
    else:

        filePrefix = fileName
        # strip the extension
        for curExtension in supportedExtensions:
            if fileName.endswith("." + curExtension):
                filePrefix = filePrefix[0:-(len(curExtension) + 1)]
                break

        for curExtension in supportedExtensions:
            if os.path.isfile(filePrefix + '.' + curExtension):
                return filePrefix

        return None

def gaussianFWHM2D(SIGMA):

    assert(isinstance(SIGMA, numpy.ndarray)),"SIGMA must be an array"

    DetSIGMA = SIGMA[0, 0] * SIGMA[1, 1] - SIGMA[1, 0] * SIGMA[1, 0]

    assert(DetSIGMA > 0),"SIGMA must be positive-definite"

    precisionMatrix = numpy.array([[SIGMA[1, 1], -SIGMA[0, 1]], [-SIGMA[1, 0], SIGMA[0, 0]]]) / DetSIGMA

    XWidth = numpy.ceil(numpy.abs(SIGMA[0, 0]) / 3);
    YWidth = numpy.ceil(numpy.abs(SIGMA[1, 1]) / 3);

    xx = numpy.arange(-XWidth, XWidth + 1)
    yy = numpy.arange(-YWidth, YWidth + 1)

    X, Y = numpy.meshgrid(xx, yy)

    XY = numpy.concatenate((numpy.atleast_2d(numpy.ravel(X)), numpy.atleast_2d(numpy.ravel(Y))), axis = 0)

    maximum = numpy.sqrt(DetSIGMA) / (2.0 * numpy.pi)
    halfMaximum = maximum / 2.0

    quadForm = -0.5 * (numpy.matrix(XY.T) * numpy.matrix(precisionMatrix))

    F = numpy.sum(numpy.array(quadForm) * numpy.array(XY.T), axis = 1)
    F = numpy.reshape(F, (numpy.size(yy), numpy.size(xx)))
    F = numpy.exp(F) * maximum

    return (F, halfMaximum)
#YWidth = ceil(abs(SIGMA(2, 2)) / 3);
#xx = -XWidth:XWidth;
#yy = -YWidth:YWidth;
#
#[X, Y] = meshgrid(xx, yy);
#
#XY = [X(:)'; Y(:)'];
#clear X Y;
#
#Maximum = 1 ./ ((2 * pi) .* sqrt(DetSIGMA));
#HalfMaximum = Maximum / 2;
#QuadForm = -0.5 * (XY' * PrecisionMatrix);
#
#F = sum(QuadForm .* XY', 2);
#F = reshape(F, length(yy), length(xx));
#F = exp(F) .* Maximum;

def parcellationStats(labelImage, pixelArea, arcLengths, arcLengthLabels):

    assert(numpy.size(arcLengthLabels) == numpy.size(arcLengths)),"arcLengths and labels must be the same size"
    numLabels = numpy.max(arcLengthLabels)

    statsToCompute = ['min', 'max', 'mean', 'median', 'std', 'var', 'area']

    STATS = dict()

    for z in range(len(statsToCompute)):
        STATS[statsToCompute[z]] = numpy.zeros(numLabels)

    for z in range(1, numLabels + 1):
        I = numpy.where(arcLengthLabels == z)
        if numpy.size(I) > 0:
            STATS['min'][z - 1] = numpy.min(arcLengths[I])
            STATS['max'][z - 1] = numpy.max(arcLengths[I])
            STATS['mean'][z - 1] = numpy.mean(arcLengths[I])
            STATS['median'][z - 1] = numpy.median(arcLengths[I])
            STATS['std'][z - 1] = numpy.std(arcLengths[I])
            STATS['var'][z - 1] = numpy.var(arcLengths[I])
        STATS['area'][z - 1] = numpy.count_nonzero(labelImage == z) * pixelArea
    return STATS

def normPDF(X, MU, SIGMA):
    XC = X - MU

    return numpy.exp(-XC * XC / SIGMA / SIGMA / 2.0) / numpy.sqrt(2.0 * numpy.pi) / SIGMA

def empty2DList(SZ):
    #assert(isinstance(SZ, tuple)),"SZ must be a tuple"
    assert(len(SZ) == 2),"SZ must have 2 elements"

    I = list()
    for z in range(SZ[0]):
        I.append([None,] * SZ[1])
    return I

def removeWhiteRowsColsPNG(IMGFileName, padding = (10, 10), blackImage = False):
    assert(os.path.isfile(IMGFileName)),"PNG file does not exist"

    IMG = imageio.imread(IMGFileName)

    if blackImage == True:
        IMG = 1 - IMG

    if IMG.ndim < 3:
        T = (IMG < 1)
    elif IMG.shape[2] == 4:
        T = numpy.squeeze(numpy.take(IMG, [0, 1, 2], axis = 2))
        T = numpy.any(T < 1, axis = 2)
    elif IMG.shape[2] == 3:
        T = numpy.any(IMG < 1, axis = 2)

    I = numpy.where(T)

    croppedIMG = numpy.array(IMG)

    for z in range(2):
        H = numpy.bincount(I[z])
        croppedIMG = numpy.take(croppedIMG, numpy.where(H > 0)[0], axis = z)

    if IMG.ndim < 3:
        cornerPadding = numpy.ones((padding[0], padding[1]))
        topBottomPadding = numpy.ones((padding[0], croppedIMG.shape[1]))
        leftRightPadding = numpy.ones((croppedIMG.shape[0], padding[1]))
    else:
        cornerPadding = numpy.ones((padding[0], padding[1], croppedIMG.shape[2]))
        topBottomPadding = numpy.ones((padding[0], croppedIMG.shape[1], croppedIMG.shape[2]))
        leftRightPadding = numpy.ones((croppedIMG.shape[0], padding[1], croppedIMG.shape[2]))

    T = numpy.concatenate((
    numpy.concatenate((cornerPadding, topBottomPadding, cornerPadding), axis = 1),
    numpy.concatenate((leftRightPadding, croppedIMG, leftRightPadding), axis = 1),
    numpy.concatenate((cornerPadding, topBottomPadding, cornerPadding), axis = 1)), axis = 0)

    if blackImage == True:
        T = 1 - T

    imageio.imwrite(IMGFileName, T)

def cropAutoWhitePNG(IMGFileName, padding = (10, 10)):
    assert(os.path.isfile(IMGFileName)),"PNG file does not exist"

    IMG = imageio.imread(IMGFileName)

    if IMG.shape[2] == 4:
        T = numpy.squeeze(numpy.take(IMG, [0, 1, 2], axis = 2))
        T = numpy.any(T < 1, axis = 2)
    elif IMG.shape[2] == 3:
        T = numpy.any(IMG < 1, axis = 2)
    elif IMG.ndim < 3:
        T = (IMG < 1)

    I = numpy.where(T)

    croppedIMG = numpy.array(IMG)

    for z in range(2):
        croppedIMG = numpy.take(croppedIMG, numpy.arange(numpy.min(I[z]), numpy.max(I[z]) + 1), axis = z)

    cornerPadding = numpy.ones((padding[0], padding[1], croppedIMG.shape[2]))
    topBottomPadding = numpy.ones((padding[0], croppedIMG.shape[1], croppedIMG.shape[2]))
    leftRightPadding = numpy.ones((croppedIMG.shape[0], padding[1], croppedIMG.shape[2]))

    T = numpy.concatenate((
    numpy.concatenate((cornerPadding, topBottomPadding, cornerPadding), axis = 1),
    numpy.concatenate((leftRightPadding, croppedIMG, leftRightPadding), axis = 1),
    numpy.concatenate((cornerPadding, topBottomPadding, cornerPadding), axis = 1)), axis = 0)

    imageio.imwrite(IMGFileName, T)

def MNI152FLIRTSymNeckCropTemplate(skullStripped = False):
    scriptPath = os.path.realpath(__file__)
    (head, tail) = os.path.split(scriptPath)

    if skullStripped == False:
        T = 'MNI152_T1_2mm_centered'
    else:
        T = 'MNI152_T1_2mm_centered_brain'

    return os.path.join(head, 'data', T)

# if we are using the MNI template,
# this is done by detecting three asterisks
# if there is no asterisks, just return the file name without modification

def MNI152FLIRTTemplate(skullStripped = False):
    scriptPath = os.path.realpath(__file__)
    (head, tail) = os.path.split(scriptPath)

    if skullStripped == False:
        T = 'MNI152_T1_1mm_centered.nii.gz'
    else:
        T = 'MNI152_T1_1mm_centered_brain.nii.gz'

    return os.path.join(head, 'data', T)
#import CCSegUtilsInterpCython

#@profile
def interp3q(xx, yy, zz, V, xi, yi, zi, interpmethod = 'linear', extrapval = numpy.nan):

    assert(numpy.size(xx) == V.shape[1] and numpy.size(yy) == V.shape[0] and numpy.size(zz) == V.shape[2]),"numpy.size(xx) == V.shape[1] and numpy.size(yy) == V.shape[0] and numpy.size(zz) == V.shape[2]"
    assert(numpy.array_equal(xi.shape, yi.shape)),"(numpy.array_equal(xi.shape, yi.shape)"
    assert(numpy.array_equal(xi.shape, zi.shape)),"(numpy.array_equal(xi.shape, zi.shape)"

    #outV = numpy.zeros_like(xi)
    outV = numpy.empty_like(xi)
    outV.fill(extrapval)

    #outOfMask = numpy.logical_or(numpy.logical_or(numpy.logical_or(numpy.logical_or(numpy.logical_or(xi < xx[0], xi > xx[-1]), yi < yy[0]), yi > yy[-1]), zi < zz[0]), zi > zz[-1])

    #minxx = numpy.min(xx)
    #maxxx = numpy.max(xx)

    #outOfMask = numpy.logical_or(numpy.logical_or(numpy.logical_or(numpy.logical_or(numpy.logical_or(xi < numpy.min(xx), xi > numpy.max(xx)), yi < numpy.min(yy)), yi > numpy.max(yy)), zi < numpy.min(zz)), zi > numpy.max(zz))

    inMaskIDX = numpy.logical_and(xi >= numpy.min(xx), xi <= numpy.max(xx))
    inMaskIDX = numpy.logical_and(inMaskIDX, yi >= numpy.min(yy))
    inMaskIDX = numpy.logical_and(inMaskIDX, yi <= numpy.max(yy))
    inMaskIDX = numpy.logical_and(inMaskIDX, zi >= numpy.min(zz))
    inMaskIDX = numpy.logical_and(inMaskIDX, zi <= numpy.max(zz))

    #outV[numpy.where(outOfMask)] = extrapval
    #if numpy.any(outOfMask):
    #    outV[outOfMask] = extrapval

    #inMaskIDX = numpy.where(numpy.logical_not(outOfMask))
    #inMaskIDX = numpy.logical_not(outOfMask)
    #del outOfMask
    #F = CCSegUtilsInterpCython.interp3q(xx, yy, zz, V, xi[inMaskIDX], yi[inMaskIDX], zi[inMaskIDX], interpmethod)
    #print F
    #print inMaskIDX
    #if numpy.size(inMaskIDX) > 0:
    if numpy.any(inMaskIDX):

        # save memory, do these individually
        #XI = (xi[inMaskIDX] - xx[0]) / (xx[1] - xx[0])
        #YI = (yi[inMaskIDX] - yy[0]) / (yy[1] - yy[0])
        #ZI = (zi[inMaskIDX] - zz[0]) / (zz[1] - zz[0])


        if interpmethod == 'linear':

            numInSlice = V.shape[0] * V.shape[1]

            XI = (xi[inMaskIDX] - xx[0]) / float(xx[1] - xx[0])
            XIDX = numpy.uint32(XI)
            XFrac = XI - XIDX
            del XI
            #I = numpy.where(XI == xx.size - 1)
            I = (XIDX == xx.size - 1)

            #if not numpy.size(I) == 0:
            if numpy.any(I):
                XFrac[I] = 1.0
                XIDX[I] = xx.size - 2
            linearIDX = XIDX * V.shape[0]
            del I
            del XIDX

            YI = (yi[inMaskIDX] - yy[0]) / float(yy[1] - yy[0])
            YIDX = numpy.uint32(YI)
            YFrac = YI - YIDX
            del YI
            #I = numpy.where(YI == yy.size - 1)
            I = (YIDX == yy.size - 1)

            #if not numpy.size(I) == 0:
            if numpy.any(I):
                YFrac[I] = 1.0
                YIDX[I] = yy.size - 2
            linearIDX += YIDX
            del I
            del YIDX

            ZI = (zi[inMaskIDX] - zz[0]) / float(zz[1] - zz[0])
            ZIDX = numpy.uint32(ZI)
            ZFrac = ZI - ZIDX
            del ZI
            #I = numpy.where(YI == yy.size - 1)
            I = (ZIDX == zz.size - 1)

            #if not numpy.size(I) == 0:
            if numpy.any(I):
                ZFrac[I] = 1.0
                ZIDX[I] = zz.size - 2
            linearIDX += ZIDX * numInSlice
            del I
            del ZIDX
            #print str(XI) + " " + str(XFrac)
            #print str(YI) + " " + str(YFrac)
            #print str(ZI) + " " + str(ZFrac)

            #outV[inMaskIDX] = \
            #    (1 - YFrac) * (1 - XFrac) * (1 - ZFrac) * V[(YIDX    , XIDX    , ZIDX    )] + \
            #    (1 - YFrac) * (1 - XFrac) * (    ZFrac) * V[(YIDX    , XIDX    , ZIDX + 1)] + \
            #    (1 - YFrac) * (    XFrac) * (1 - ZFrac) * V[(YIDX    , XIDX + 1, ZIDX    )] + \
            #    (1 - YFrac) * (    XFrac) * (    ZFrac) * V[(YIDX    , XIDX + 1, ZIDX + 1)] + \
            #    (    YFrac) * (1 - XFrac) * (1 - ZFrac) * V[(YIDX + 1, XIDX    , ZIDX    )] + \
            #    (    YFrac) * (1 - XFrac) * (    ZFrac) * V[(YIDX + 1, XIDX    , ZIDX + 1)] + \
            #    (    YFrac) * (    XFrac) * (1 - ZFrac) * V[(YIDX + 1, XIDX + 1, ZIDX    )] + \
            #    (    YFrac) * (    XFrac) * (    ZFrac) * V[(YIDX + 1, XIDX + 1, ZIDX + 1)]
            flatV = numpy.ravel(V, order = 'F')
            #flatinMaskIDX = numpy.ravel(inMaskIDX, order = 'F')
            #print flatV.shape
            #outV[inMaskIDX] = (1 - YFrac) * (1 - XFrac) * (1 - ZFrac) * flatV[linearIDX] + (1 - YFrac) * (1 - XFrac) * (    ZFrac) * flatV[linearIDX + numInSlice] + (1 - YFrac) * (    XFrac) * (1 - ZFrac) * flatV[linearIDX + V.shape[0]] + (1 - YFrac) * (    XFrac) * (    ZFrac) * flatV[linearIDX + V.shape[0] + numInSlice] + (    YFrac) * (1 - XFrac) * (1 - ZFrac) * flatV[linearIDX + 1] + (    YFrac) * (1 - XFrac) * (    ZFrac) * flatV[linearIDX + 1 + numInSlice] + (    YFrac) * (    XFrac) * (1 - ZFrac) * flatV[linearIDX + 1 + V.shape[0]] + (    YFrac) * (    XFrac) * (    ZFrac) * flatV[linearIDX + 1 + V.shape[0] + numInSlice]
#
            outV[inMaskIDX] = numexpr.evaluate('(1 - YFrac) * (1 - XFrac) * (1 - ZFrac) * a + (1 - YFrac) * (1 - XFrac) * (    ZFrac) * b + (1 - YFrac) * (    XFrac) * (1 - ZFrac) * c + (1 - YFrac) * (    XFrac) * (    ZFrac) * d + (    YFrac) * (1 - XFrac) * (1 - ZFrac) * e + (    YFrac) * (1 - XFrac) * (    ZFrac) * f + (    YFrac) * (    XFrac) * (1 - ZFrac) * g + (    YFrac) * (    XFrac) * (    ZFrac) * h', {'a': flatV[linearIDX], 'b': flatV[linearIDX + numInSlice], 'c': flatV[linearIDX + V.shape[0]], 'd': flatV[linearIDX + V.shape[0] + numInSlice], 'e': flatV[linearIDX + 1], 'f': flatV[linearIDX + 1 + numInSlice], 'g': flatV[linearIDX + 1 + V.shape[0]], 'h': flatV[linearIDX + 1 + V.shape[0] + numInSlice], 'XFrac': XFrac, 'YFrac': YFrac, 'ZFrac': ZFrac})

            #I = YIDX + XIDX * V.shape[0] + ZIDX * numInSlice

#            outV[inMaskIDX] =  (1 - YFrac) * (1 - XFrac) * (1 - ZFrac) * flatV[linearIDX]
#            outV[inMaskIDX] += (1 - YFrac) * (1 - XFrac) * (    ZFrac) * flatV[linearIDX + numInSlice] # V[(YIDX    , XIDX    , ZIDX + 1)]
#            outV[inMaskIDX] += (1 - YFrac) * (    XFrac) * (1 - ZFrac) * flatV[linearIDX + V.shape[0]] # V[(YIDX    , XIDX + 1, ZIDX    )]
#            outV[inMaskIDX] += (1 - YFrac) * (    XFrac) * (    ZFrac) * flatV[linearIDX + V.shape[0] + numInSlice] # V[(YIDX    , XIDX + 1, ZIDX + 1)]
#            outV[inMaskIDX] += (    YFrac) * (1 - XFrac) * (1 - ZFrac) * flatV[linearIDX + 1] # V[(YIDX + 1, XIDX    , ZIDX    )]
#            outV[inMaskIDX] += (    YFrac) * (1 - XFrac) * (    ZFrac) * flatV[linearIDX + 1 + numInSlice] # V[(YIDX + 1, XIDX    , ZIDX + 1)]
#            outV[inMaskIDX] += (    YFrac) * (    XFrac) * (1 - ZFrac) * flatV[linearIDX + 1 + V.shape[0]] # V[(YIDX + 1, XIDX + 1, ZIDX    )]
#            outV[inMaskIDX] += (    YFrac) * (    XFrac) * (    ZFrac) * flatV[linearIDX + 1 + V.shape[0] + numInSlice] # V[(YIDX + 1, XIDX + 1, ZIDX + 1)]

#            T =  (1 - YFrac) * (1 - XFrac) * (1 - ZFrac) * flatV[linearIDX]
#            T += (1 - YFrac) * (1 - XFrac) * (    ZFrac) * flatV[linearIDX + numInSlice] # V[(YIDX    , XIDX    , ZIDX + 1)]
#            T += (1 - YFrac) * (    XFrac) * (1 - ZFrac) * flatV[linearIDX + V.shape[0]] # V[(YIDX    , XIDX + 1, ZIDX    )]
#            T += (1 - YFrac) * (    XFrac) * (    ZFrac) * flatV[linearIDX + V.shape[0] + numInSlice] # V[(YIDX    , XIDX + 1, ZIDX + 1)]
#            T += (    YFrac) * (1 - XFrac) * (1 - ZFrac) * flatV[linearIDX + 1] # V[(YIDX + 1, XIDX    , ZIDX    )]
#            T += (    YFrac) * (1 - XFrac) * (    ZFrac) * flatV[linearIDX + 1 + numInSlice] # V[(YIDX + 1, XIDX    , ZIDX + 1)]
#            T += (    YFrac) * (    XFrac) * (1 - ZFrac) * flatV[linearIDX + 1 + V.shape[0]] # V[(YIDX + 1, XIDX + 1, ZIDX    )]
#            T += (    YFrac) * (    XFrac) * (    ZFrac) * flatV[linearIDX + 1 + V.shape[0] + numInSlice] # V[(YIDX + 1, XIDX + 1, ZIDX + 1)]
#            outV[inMaskIDX] = T
#

            #T = outV[inMaskIDX]
            #outV[inMaskIDX] =  (1 - YFrac) * (1 - XFrac) * (1 - ZFrac) * V[(YIDX    , XIDX    , ZIDX    )]
            #outV[inMaskIDX] += (1 - YFrac) * (1 - XFrac) * (    ZFrac) * V[(YIDX    , XIDX    , ZIDX + 1)]
            #outV[inMaskIDX] += (1 - YFrac) * (    XFrac) * (1 - ZFrac) * V[(YIDX    , XIDX + 1, ZIDX    )]
            #outV[inMaskIDX] += (1 - YFrac) * (    XFrac) * (    ZFrac) * V[(YIDX    , XIDX + 1, ZIDX + 1)]
            #outV[inMaskIDX] += (    YFrac) * (1 - XFrac) * (1 - ZFrac) * V[(YIDX + 1, XIDX    , ZIDX    )]
            #outV[inMaskIDX] += (    YFrac) * (1 - XFrac) * (    ZFrac) * V[(YIDX + 1, XIDX    , ZIDX + 1)]
            #outV[inMaskIDX] += (    YFrac) * (    XFrac) * (1 - ZFrac) * V[(YIDX + 1, XIDX + 1, ZIDX    )]
            #outV[inMaskIDX] += (    YFrac) * (    XFrac) * (    ZFrac) * V[(YIDX + 1, XIDX + 1, ZIDX + 1)]
            #print numpy.array_equal(T, outV[inMaskIDX])
        elif interpmethod == 'nearest':
            XI = numpy.int32(numpy.round((xi[inMaskIDX] - xx[0]) / (xx[1] - xx[0])))
            YI = numpy.int32(numpy.round((yi[inMaskIDX] - yy[0]) / (yy[1] - yy[0])))
            ZI = numpy.int32(numpy.round((zi[inMaskIDX] - zz[0]) / (zz[1] - zz[0])))
            outV[inMaskIDX] = V[(YI, XI, ZI)]

    return outV

#@profile
#def interp3q2D(xx, yy, zz, V, xi, yi, zi, interpmethod = 'linear', extrapval = numpy.nan):
#
#    assert(numpy.size(xx) == V.shape[1] and numpy.size(yy) == V.shape[0] and numpy.size(zz) == V.shape[2]),"numpy.size(xx) == V.shape[1] and numpy.size(yy) == V.shape[0] and numpy.size(zz) == V.shape[2]"
#    assert(numpy.array_equal(xi.shape, yi.shape)),"(numpy.array_equal(xi.shape, yi.shape)"
#    assert(numpy.array_equal(xi.shape, zi.shape)),"(numpy.array_equal(xi.shape, zi.shape)"
#
#    outV = numpy.zeros_like(xi)
#
#    #outOfMask = numpy.logical_or(numpy.logical_or(numpy.logical_or(numpy.logical_or(numpy.logical_or(xi < xx[0], xi > xx[-1]), yi < yy[0]), yi > yy[-1]), zi < zz[0]), zi > zz[-1])
#
#    #minxx = numpy.min(xx)
#    #maxxx = numpy.max(xx)
#
#    outOfMask = numpy.logical_or(numpy.logical_or(numpy.logical_or(numpy.logical_or(numpy.logical_or(xi < numpy.min(xx), xi > numpy.max(xx)), yi < numpy.min(yy)), yi > numpy.max(yy)), zi < numpy.min(zz)), zi > numpy.max(zz))
#
#    #outV[numpy.where(outOfMask)] = extrapval
#    if numpy.any(outOfMask):
#        outV[outOfMask] = extrapval
#
#    #inMaskIDX = numpy.where(numpy.logical_not(outOfMask))
#    inMaskIDX = numpy.logical_not(outOfMask)
#    del outOfMask
#    #F = CCSegUtilsInterpCython.interp3q(xx, yy, zz, V, xi[inMaskIDX], yi[inMaskIDX], zi[inMaskIDX], interpmethod)
#    #print F
#    #print inMaskIDX
#    #if numpy.size(inMaskIDX) > 0:
#    if numpy.any(inMaskIDX):
#
#        # save memory, do these individually
#        #XI = (xi[inMaskIDX] - xx[0]) / (xx[1] - xx[0])
#        #YI = (yi[inMaskIDX] - yy[0]) / (yy[1] - yy[0])
#        #ZI = (zi[inMaskIDX] - zz[0]) / (zz[1] - zz[0])
#
#
#        if interpmethod == 'nearest':
#            XI = numpy.int32(numpy.round((xi[inMaskIDX] - xx[0]) / (xx[1] - xx[0])))
#            YI = numpy.int32(numpy.round((yi[inMaskIDX] - yy[0]) / (yy[1] - yy[0])))
#            ZI = numpy.int32(numpy.round((zi[inMaskIDX] - zz[0]) / (zz[1] - zz[0])))
#            outV[inMaskIDX] = V[(YI, XI, ZI)]
#        elif interpmethod == 'linear':
#
#            numInSlice = V.shape[0] * V.shape[1]
#
#            XI = (xi[inMaskIDX] - xx[0]) / (xx[1] - xx[0])
#            XIDX = numpy.uint32(XI)
#            XFrac = XI - XIDX
#            del XI
#            #I = numpy.where(XI == xx.size - 1)
#            I = (XIDX == xx.size - 1)
#
#            #if not numpy.size(I) == 0:
#            if numpy.any(I):
#                XFrac[I] = 1.0
#                XIDX[I] = xx.size - 2
#            linearIDX = numpy.array(XIDX) * V.shape[0]
#            del I
#            del XIDX
#
#            YI = (yi[inMaskIDX] - yy[0]) / (yy[1] - yy[0])
#            YIDX = numpy.uint32(YI)
#            YFrac = YI - YIDX
#            del YI
#            #I = numpy.where(YI == yy.size - 1)
#            I = (YIDX == yy.size - 1)
#
#            #if not numpy.size(I) == 0:
#            if numpy.any(I):
#                YFrac[I] = 1.0
#                YIDX[I] = yy.size - 2
#            linearIDX += numpy.array(YIDX)
#            del I
#            del YIDX
#
#            ZI = (zi[inMaskIDX] - zz[0]) / (zz[1] - zz[0])
#            ZIDX = numpy.uint32(ZI)
#            ZFrac = ZI - ZIDX
#            del ZI
#            #I = numpy.where(YI == yy.size - 1)
#            I = (ZIDX == zz.size - 1)
#
#            #if not numpy.size(I) == 0:
#            if numpy.any(I):
#                ZFrac[I] = 1.0
#                ZIDX[I] = zz.size - 2
#            linearIDX += numpy.array(ZIDX) * V.shape[0]    * V.shape[1]
#            del I
#            del ZIDX
#            #print str(XI) + " " + str(XFrac)
#            #print str(YI) + " " + str(YFrac)
#            #print str(ZI) + " " + str(ZFrac)
#
#            #outV[inMaskIDX] = \
#            #    (1 - YFrac) * (1 - XFrac) * (1 - ZFrac) * V[(YIDX    , XIDX    , ZIDX    )] + \
#            #    (1 - YFrac) * (1 - XFrac) * (    ZFrac) * V[(YIDX    , XIDX    , ZIDX + 1)] + \
#            #    (1 - YFrac) * (    XFrac) * (1 - ZFrac) * V[(YIDX    , XIDX + 1, ZIDX    )] + \
#            #    (1 - YFrac) * (    XFrac) * (    ZFrac) * V[(YIDX    , XIDX + 1, ZIDX + 1)] + \
#            #    (    YFrac) * (1 - XFrac) * (1 - ZFrac) * V[(YIDX + 1, XIDX    , ZIDX    )] + \
#            #    (    YFrac) * (1 - XFrac) * (    ZFrac) * V[(YIDX + 1, XIDX    , ZIDX + 1)] + \
#            #    (    YFrac) * (    XFrac) * (1 - ZFrac) * V[(YIDX + 1, XIDX + 1, ZIDX    )] + \
#            #    (    YFrac) * (    XFrac) * (    ZFrac) * V[(YIDX + 1, XIDX + 1, ZIDX + 1)]
#
#            flatV = numpy.ravel(V, order = 'F')
#
#            #I = YIDX + XIDX * V.shape[0] + ZIDX * numInSlice
#
#            outV[inMaskIDX] =  (1 - YFrac) * (1 - XFrac) * (1 - ZFrac) * flatV[linearIDX]
#            outV[inMaskIDX] += (1 - YFrac) * (1 - XFrac) * (    ZFrac) * flatV[linearIDX + numInSlice] # V[(YIDX    , XIDX    , ZIDX + 1)]
#            outV[inMaskIDX] += (1 - YFrac) * (    XFrac) * (1 - ZFrac) * flatV[linearIDX + V.shape[0]] # V[(YIDX    , XIDX + 1, ZIDX    )]
#            outV[inMaskIDX] += (1 - YFrac) * (    XFrac) * (    ZFrac) * flatV[linearIDX + V.shape[0] + numInSlice] # V[(YIDX    , XIDX + 1, ZIDX + 1)]
#            outV[inMaskIDX] += (    YFrac) * (1 - XFrac) * (1 - ZFrac) * flatV[linearIDX + 1] # V[(YIDX + 1, XIDX    , ZIDX    )]
#            outV[inMaskIDX] += (    YFrac) * (1 - XFrac) * (    ZFrac) * flatV[linearIDX + 1 + numInSlice] # V[(YIDX + 1, XIDX    , ZIDX + 1)]
#            outV[inMaskIDX] += (    YFrac) * (    XFrac) * (1 - ZFrac) * flatV[linearIDX + 1 + V.shape[0]] # V[(YIDX + 1, XIDX + 1, ZIDX    )]
#            outV[inMaskIDX] += (    YFrac) * (    XFrac) * (    ZFrac) * flatV[linearIDX + 1 + V.shape[0] + numInSlice] # V[(YIDX + 1, XIDX + 1, ZIDX + 1)]
#
#            #T = outV[inMaskIDX]
#            #outV[inMaskIDX] =  (1 - YFrac) * (1 - XFrac) * (1 - ZFrac) * V[(YIDX    , XIDX    , ZIDX    )]
#            #outV[inMaskIDX] += (1 - YFrac) * (1 - XFrac) * (    ZFrac) * V[(YIDX    , XIDX    , ZIDX + 1)]
#            #outV[inMaskIDX] += (1 - YFrac) * (    XFrac) * (1 - ZFrac) * V[(YIDX    , XIDX + 1, ZIDX    )]
#            #outV[inMaskIDX] += (1 - YFrac) * (    XFrac) * (    ZFrac) * V[(YIDX    , XIDX + 1, ZIDX + 1)]
#            #outV[inMaskIDX] += (    YFrac) * (1 - XFrac) * (1 - ZFrac) * V[(YIDX + 1, XIDX    , ZIDX    )]
#            #outV[inMaskIDX] += (    YFrac) * (1 - XFrac) * (    ZFrac) * V[(YIDX + 1, XIDX    , ZIDX + 1)]
#            #outV[inMaskIDX] += (    YFrac) * (    XFrac) * (1 - ZFrac) * V[(YIDX + 1, XIDX + 1, ZIDX    )]
#            #outV[inMaskIDX] += (    YFrac) * (    XFrac) * (    ZFrac) * V[(YIDX + 1, XIDX + 1, ZIDX + 1)]
#            #print numpy.array_equal(T, outV[inMaskIDX])
#    return outV
#

#if
#V = numpy.random.uniform(size = (4, 6, 7))

#xx = numpy.arange(V.shape[1])
#yy = numpy.arange(V.shape[0])
#zz = numpy.arange(V.shape[2]) * 5

#xi = numpy.random.uniform(size = (250, 250, 250))
#yi = numpy.random.uniform(size = (250, 250, 250))
#zi = numpy.random.uniform(size = (250, 250, 250))

#xi = numpy.array([-1, 0, 3])
#yi = numpy.array([-1, 0, 3])
#zi = numpy.array([-1, 1.4, 3])
#G = interp3q(xx, yy, zz, V, xi, yi, zi)
#@profile
def interp2q(xx, yy, V, xi, yi, interpmethod = 'linear', extrapval = numpy.nan):

    assert(numpy.size(xx) == V.shape[1] and numpy.size(yy) == V.shape[0]),"numpy.size(xx) == V.shape[1] and numpy.size(yy) == V.shape[0]"
    assert(numpy.array_equal(xi.shape, yi.shape)),"(numpy.array_equal(xi.shape, yi.shape)"
#    print "xx = "
#    print xx
#    print "yy = "
#    print yy
#    print "xi = "
#    print xi
#    print "yi = "
#    print yi
#
    outV = numpy.empty_like(xi)
    outV.fill(extrapval)

    inMaskIDX = numpy.logical_and(xi >= numpy.min(xx), xi <= numpy.max(xx))
    inMaskIDX = numpy.logical_and(inMaskIDX, yi >= numpy.min(yy))
    inMaskIDX = numpy.logical_and(inMaskIDX, yi <= numpy.max(yy))

    #outV[numpy.where(outOfMask)] = extrapval
    #f numpy.any(outOfMask):
    #outV[outOfMask] = extrapval

    #inMaskIDX = numpy.where(numpy.logical_not(outOfMask))

    #if numpy.size(inMaskIDX) > 0:
    if numpy.any(inMaskIDX):

        #xSpacing = xx[1] - xx[0]
        #ySpacing = yy[1] - yy[0]
        XI = (xi[inMaskIDX] - xx[0]) / float(xx[1] - xx[0])
        YI = (yi[inMaskIDX] - yy[0]) / float(yy[1] - yy[0])

        if interpmethod == 'nearest':
            XI = numpy.int32(numpy.round(XI))
            YI = numpy.int32(numpy.round(YI))
            outV[inMaskIDX] = V[(YI, XI)]
        elif interpmethod == 'linear':

            XFrac = XI - numpy.floor(XI)
            #I = numpy.where(XI == xx.size - 1)
            I = (XI == xx.size - 1)

            #if not numpy.size(I) == 0:
            if numpy.any(I):
                XFrac[I] = 1.0
                XI[I] = xx.size - 2
            linearIDX = XI * V.shape[0]

            YFrac = YI - numpy.floor(YI)
            #I = numpy.where(YI == yy.size - 1)
            I = (YI == yy.size - 1)

            #if not numpy.size(I) == 0:
            if numpy.any(I):
                YFrac[I] = 1.0
                YI[I] = yy.size - 2

            XI = numpy.int32(XI)
            YI = numpy.int32(YI)
            #linearIDX = YI * V.shape[0] + XI

            #flatV = numpy.ravel(V, order = 'F')

            outV[inMaskIDX] =    (1 - YFrac) * (1 - XFrac) * V[(YI    , XI    )] + (1 - YFrac) * (    XFrac) * V[(YI    , XI + 1)] + (    YFrac) * (1 - XFrac) * V[(YI + 1, XI    )] + (    YFrac) * (    XFrac) * V[(YI + 1, XI + 1)]
            #outV[inMaskIDX] =    numexpr.evaluate('(1 - YFrac) * (1 - XFrac) * a + (1 - YFrac) * (    XFrac) * b + (    YFrac) * (1 - XFrac) * c + (    YFrac) * (    XFrac) * d', {'a': flatV[linearIDX], 'b': flatV[linearIDX + 1], 'c': flatV[linearIDX + V.shape[1]], 'd': flatV[linearIDX + V.shape[1] + 1], 'XFrac': XFrac, 'YFrac': YFrac})
    return outV

def maxGaussian1D(SIGMA, Derivative = 0):
    GaussianDieOff = 0.0001;

    SIGMASQ = SIGMA * SIGMA;

    W = numpy.arange(1, 501)

    FirstGaussian = numpy.exp(-(W * W) / (2 * SIGMASQ))

    #MaxWidth = find(FirstGaussian > GaussianDieOff, 1, 'last');
    MaxWidth = numpy.where(FirstGaussian > GaussianDieOff);

    if(numpy.size(MaxWidth) == 0):
        MaxWidth = 1;
    else:
        MaxWidth = numpy.size(MaxWidth)

    X = numpy.arange(-MaxWidth, MaxWidth + 1)
    Exponential = numpy.exp(-(X * X) / (2 * SIGMASQ))

    if Derivative == 0:
        return Exponential / (numpy.sqrt(2 * numpy.pi) * SIGMA)
    elif Derivative == 1:
        return Exponential * (-X / (numpy.sqrt(2 * numpy.pi) * SIGMA * SIGMASQ))

def mkdirSafe(D):
    try:
        os.makedirs(D)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(D):
            pass
        else:
            raise Exception

def parasagittalSlicesAndGradients(TIMG, axialNIIPixdims, numSlices = 3):

    midSlice = int(math.floor(TIMG.shape[0] / 2))

    # store the parasagittal slices
    parasagittalIDX = list(range(numpy.maximum(0, midSlice - numSlices), numpy.minimum(TIMG.shape[0], midSlice + numSlices + 1)))
    #print parasagittalIDX
    parasagittalSlices = numpy.take(TIMG, parasagittalIDX, axis = 0)
#
    # get gradients
    SIGMA = 0.5 / numpy.array(axialNIIPixdims)
    STIMG = scipy.ndimage.filters.gaussian_filter(TIMG, SIGMA, mode = 'constant', cval = 0)

    gaussianFilterDeriv = numpy.atleast_3d(maxGaussian1D(SIGMA[0], Derivative = 1))
    T = numpy.reshape(gaussianFilterDeriv, (gaussianFilterDeriv.size, 1, 1))
    STIMGFX = scipy.ndimage.convolve(STIMG, T, mode = 'nearest')

    gaussianFilterDeriv = numpy.atleast_3d(maxGaussian1D(SIGMA[1], Derivative = 1))
    T = numpy.reshape(gaussianFilterDeriv, (1, gaussianFilterDeriv.size, 1))
    STIMGFY = scipy.ndimage.convolve(STIMG, T, mode = 'nearest')

    gaussianFilterDeriv = numpy.atleast_3d(maxGaussian1D(SIGMA[2], Derivative = 1))
    T = numpy.reshape(gaussianFilterDeriv, (1, 1, gaussianFilterDeriv.size))
    STIMGFZ = scipy.ndimage.convolve(STIMG, T, mode = 'nearest')

    parasagittalFX = numpy.take(STIMGFX, parasagittalIDX, axis = 0)
    parasagittalFY = numpy.take(STIMGFY, parasagittalIDX, axis = 0)
    parasagittalFZ = numpy.take(STIMGFZ, parasagittalIDX, axis = 0)

    P = [1, 2, 0]
    parasagittalSlices = numpy.transpose(parasagittalSlices, P)
    parasagittalFX = numpy.transpose(parasagittalFX, P)
    parasagittalFY = numpy.transpose(parasagittalFY, P)
    parasagittalFZ = numpy.transpose(parasagittalFZ, P)

    parasagittalSlices = numpy.rot90(parasagittalSlices, 1)
    parasagittalFX = numpy.rot90(parasagittalFX, 1)
    parasagittalFY = numpy.rot90(parasagittalFY, 1)
    parasagittalFZ = numpy.rot90(parasagittalFZ, 1)

    return (parasagittalSlices, parasagittalFX, parasagittalFY, parasagittalFZ)
