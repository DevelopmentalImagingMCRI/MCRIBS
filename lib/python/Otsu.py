#!/usr/bin/python

import numpy

# sets the threshold THRESH so that class1 <= THRESH < class2
def otsu1(IMG):
    hist, bin_edges = numpy.histogram(IMG, bins=range(257), range=None, normed=False, weights=None, density=True)
    
    OmegaZeros = numpy.cumsum(hist)
    #OmegaK = numpy.array(OmegaZeros)
    # includes the current element in the sum , OmegaOnes[I] = sum(hist[I:end])
    #OmegaOnes = numpy.cumsum(hist[::-1])[::-1]
    # excludes the current element in the sum OmegaOnes[I] = sum(hist[I + 1:end])
    OmegaOnes = 1 - OmegaZeros
    
    #rint OmegaZeros
    #print OmegaOnes
    #print T

    OmegaMask = numpy.logical_and(OmegaZeros > 0, OmegaOnes > 0)
    OmegaMask[-1] = False
    OmegaMask = numpy.logical_not(OmegaMask)

    XTimesHist = numpy.arange(256) * hist
    CumSumXTimesHist = numpy.cumsum(XTimesHist)
    
    OmegaZeros[OmegaMask] = 1
    OmegaOnes[OmegaMask] = 1
    
    MuZeros = CumSumXTimesHist / OmegaZeros
    # this excludes the left element
    MuOnes = (CumSumXTimesHist[-1] - CumSumXTimesHist) / OmegaOnes
    #numpy.cumsum(XTimesHist[::-1])[::-1]
    
    T = (MuOnes - MuZeros)

    SigmaB = OmegaZeros * OmegaOnes * T * T
    SigmaB[OmegaMask] = 0
    return numpy.argmax(SigmaB)
    
# fast reimplemented version
def otsu2(IMG, returnWorkingValues = False):
    
    hist, bin_edges = numpy.histogram(IMG, bins=range(257), range=None, normed=False, weights=None, density=True)
    CumSumHist = numpy.cumsum(hist)
    
    CumSumHistBack = numpy.cumsum(hist[::-1])
    
    ARange = numpy.arange(256)
    XTimesHist = ARange * hist

    #MuZero = (1 / OmegaZero) * numpy.sum(XTimesHist[0:(I + 1)])
    
    # the thresholds are compted so that it includes the element on the right, so
    # X <= THRESH(0) < X <= THRESH(1) < X

    OmegaZeros = numpy.array(CumSumHist)
    # this excludes the left element
    OmegaOnes = numpy.triu(numpy.atleast_2d(CumSumHist) - numpy.atleast_2d(CumSumHist).T)
    
    # put ones on the lower triangular in order to prevent divide by zero errors
    OmegaOnes = OmegaOnes + numpy.tril(numpy.ones(OmegaOnes.shape))
    OmegaTwos = 1 - OmegaZeros
    #print OmegaTwos    == 0
    #OmegaZerosRepeated = numpy.tile(numpy.atleast_2d(OmegaZeros).T, (1, 256))
    #OmegaTwosRepeated = numpy.tile(numpy.atleast_2d(OmegaTwos), (256, 1))
    # the last element of T will be False, due to machine precision, despite OmegaZeros[-1] == 1, so manually set it to True
    T = (OmegaTwos == 0)
    T[-1] = True
    #print T
    OmegaMask = numpy.logical_or(numpy.atleast_2d(OmegaZeros).T == 0, numpy.logical_or(OmegaOnes == 0, numpy.atleast_2d(T)))
    del T
    #print numpy.where(OmegaMask)

    OmegaZeros[OmegaZeros == 0] = 1
    OmegaOnes[OmegaMask] = 1
    OmegaTwos[OmegaTwos == 0] = 1

    CumSumXTimesHist = numpy.cumsum(XTimesHist)
    CumSumXTimesHistGrid = numpy.triu(numpy.atleast_2d(CumSumXTimesHist) - numpy.atleast_2d(CumSumXTimesHist).T, 1) 
    
    CumSumXTimesHistBack = numpy.cumsum(XTimesHist[::-1])[::-1] 
    # this removes the left element of the 
    # CumSumXTimesHistBack[I] = numpy.sum(XTimesHist[I + 1:])
    CumSumXTimesHistBack = CumSumXTimesHistBack - XTimesHist

    MuZeros = CumSumXTimesHist / OmegaZeros
    # this adds in the left most element to the range
    MuOnes = CumSumXTimesHistGrid / OmegaOnes
    #MuOnes = numpy.triu(MuOnes, 1) + numpy.tril(numpy.ones(MuOnes.shape))
    
    #MuOneHistArray[I, J] = numpy.sum(XTimesHist[I:(J + 1)])
    MuTwos = CumSumXTimesHistBack / OmegaTwos
    
    #MuZeros = numpy.zeros((256))
    #MuOnes = numpy.zeros((256, 256))
    #MuTwos = numpy.zeros((256, 256))
    
    #SigmaZeros = numpy.zeros((256))
    #SigmaOnes = numpy.zeros((256, 256))
    #SigmaTwos = numpy.zeros((256, 256))
    
    #SigmaBs = numpy.zeros((256, 256))
#
#    #MaxSigmaB = 0
#    #THRESH = None
#
#    IDXForXC = numpy.uint8(numpy.arange(256))
#
#    colIDX = numpy.tile(numpy.atleast_2d(IDXForXC).T, (1, 256))
#    
#    # for the first threshold, the valid values are less than the current threshold
#    # so we use upper triangular thresholding on colIDX
#    validColIDX = numpy.triu(numpy.ones(colIDX.shape, dtype = numpy.bool)) > 0
#    
#    # each column contains the X index minus MuZero for that threshold
#    XCZeros = (colIDX - numpy.atleast_2d(MuZeros)) * numpy.double(validColIDX)
#    
#    # now we need to multiply and sum with the histogram
#    SigmaZeros = numpy.sum((XCZeros * XCZeros) * numpy.atleast_2d(hist).T, axis = 0) / OmegaZeros
#    
#    # for the second threshold, the valid values are greater than the current threshold
#    # so we use lower triangular thresholding on colIDX
#    validColIDX = numpy.tril(numpy.ones(colIDX.shape, dtype = numpy.bool)) > 0
#    
#    # each column contains the X index minus MuZero for that threshold
#    XCZeros = (colIDX - numpy.atleast_2d(MuTwos)) * numpy.double(validColIDX)
#    
#    # now we need to multiply and sum with the histogram
#    SigmaTwos = numpy.sum((XCZeros * XCZeros) * numpy.atleast_2d(hist).T, axis = 0) / OmegaTwos
#    
#    del colIDX
#    del XCZeros
#    del validColIDX
#
#    # to create the index arrays for the two threshold case, we create row, column and slice index arrays
#    T = numpy.uint8(numpy.arange(256))
#    #import time
#
#    #start_time = time.time()
#    #rowIDX, colIDX, sliceIDX = numpy.meshgrid(T, T, T, indexing = 'ij')
#
#    #validIDX = numpy.logical_and(sliceIDX >= rowIDX, sliceIDX <= colIDX)
#    #elapsed_time = time.time() - start_time
#    #print "Meshgrid version: " + str(elapsed_time)
#    #del rowIDX
#    #del colIDX
#    #del sliceIDX
#    
#    # WE DONT NEED TO COMPUTE THE SIGMAS
#    # the broadcast version is 100 times faster than using meshgrid
#    #start_time = time.time()
#    
#    rowIDX = numpy.reshape(numpy.atleast_3d(T), (256, 1, 1))
#    colIDX = numpy.reshape(numpy.atleast_3d(T), (1, 256, 1))
#    sliceIDX = numpy.reshape(numpy.atleast_3d(T), (1, 1, 256))
#    
#    validIDX = numpy.logical_and(sliceIDX >= rowIDX, sliceIDX <= colIDX)
#    
#    del rowIDX
#    del colIDX
#    
##    elapsed_time = time.time() - start_time
#    
##    print "Broadcast version: " + str(elapsed_time)
#    
##    print numpy.array_equal(validIDX, validIDX2)
#
#    del T
#    
#    #print sliceIDX.shape
#    XC = (numpy.atleast_3d(MuOnes) - sliceIDX) * numpy.double(validIDX)
#    del sliceIDX
#    
#    del validIDX
#    
#    SigmaOnes = numpy.triu(numpy.sum((XC * XC) * numpy.reshape(numpy.atleast_3d(hist), (1, 1, 256)), axis = 2) / OmegaOnes, 1)
#        
#    del XC
#
    # replicate the zeros (as a column, repeat in axis = 1) and twos (as a row, repeat in axis = 0)
    # mask them with the OmegaMasks to mask out cases where the original omega values were zeros

    # then compute SigmaB and get the maximum

    #C = (colIDX - numpy.atleast_2d(MuTwos)) * numpy.double(validColIDX)

    #print T[0:10]
    #olIDX = numpy.tril(colIDX, 1)

    #XC = ARange[0:(I + 1)] - MuZero
    #SigmaZero = (1.0 / OmegaZero) * numpy.sum((XC * XC) * hist[0:(I + 1)])
    
    #MuT = numpy.atleast_2d(MuZeros).T * numpy.atleast_2d(OmegaZeros).T + MuOnes * OmegaOnes + numpy.atleast_2d(MuTwos) * numpy.atleast_2d(OmegaTwos)
    
    #MuT = numpy.triu(MuT, 1) * numpy.double(numpy.logical_not(OmegaMask))
    #print MuT
    MuT = CumSumXTimesHist[-1]
    #print MuT
    TZeros = numpy.atleast_2d(MuZeros).T - MuT
    TOnes = MuOnes - MuT
    TTwos = numpy.atleast_2d(MuTwos) - MuT

    SigmaB = numpy.atleast_2d(OmegaZeros).T * TZeros * TZeros + OmegaOnes * TOnes * TOnes + numpy.atleast_2d(OmegaTwos) * TTwos * TTwos
    SigmaB = numpy.triu(SigmaB, 1) * numpy.double(numpy.logical_not(OmegaMask))
    #SumMu = numpy.atleast_2d(OmegaZeros).T + OmegaOnes + numpy.atleast_2d(OmegaTwos)
    #print numpy.triu(MuT, 1)
    #print numpy.sum(XTimesHist)

    #T = (numpy.atleast_2d(MuZeros).T - MuOnes - numpy.atleast_2d(MuTwos))

    #SigmaB = numpy.atleast_2d(OmegaZeros).T * OmegaOnes * numpy.atleast_2d(OmegaTwos) * T * T
    
    #SigmaB = numpy.triu(SigmaB, 1) * numpy.double(numpy.logical_not(OmegaMask))
    
    # unravel_index is the equivalent of ind2sub
    THRESH = numpy.unravel_index(numpy.argmax(SigmaB), SigmaB.shape)
    #THRESH2 = numpy.unravel_index(numpy.argmax(SigmaB2), SigmaB.shape)
    
    #print THRESH
    #print THRESH2
    #print SigmaB
    #print I

#    print "MuZeros: " + str(MuZeros.shape)
#    print "OmegaZeros: " + str(OmegaZeros.shape)
#    print "SigmaZeros: " + str(SigmaZeros.shape)
#    
#    print "MuOnes: " + str(MuOnes.shape)
#    print "OmegaOnes: " + str(OmegaOnes.shape)
#    print "SigmaOnes: " + str(SigmaOnes.shape)
#    
#    print "MuTwos: " + str(MuTwos.shape)
#    print "OmegaTwos: " + str(OmegaTwos.shape)
#    print "SigmaTwos: " + str(SigmaTwos.shape)
#    
    if returnWorkingValues:
        SigmaZeros = None
        SigmaOnes = None
        SigmaTwos = None
        SigmaB2 = numpy.array(MuT)
        return (THRESH, OmegaZeros, MuZeros, SigmaZeros, OmegaOnes, MuOnes, SigmaOnes, OmegaTwos, MuTwos, SigmaTwos, CumSumXTimesHistGrid, SigmaB, SigmaB2)
    else:
        return THRESH

# slow version
def otsu2Loops(IMG, returnWorkingValues = False):
    
    hist, bin_edges = numpy.histogram(IMG, bins=range(257), range=None, normed=False, weights=None, density=True)
    CumSumHist = numpy.cumsum(hist)
    
    ARange = numpy.arange(256)
    XTimesHist = ARange * hist
    
    OmegaZeros = numpy.zeros((256))
    OmegaOnes = numpy.zeros((256, 256))
    OmegaTwos = numpy.zeros((256, 256))
    
    MuZeros = numpy.zeros((256))
    MuOnes = numpy.zeros((256, 256))
    MuTwos = numpy.zeros((256, 256))

    #SigmaZeros = numpy.zeros((256))
    #SigmaOnes = numpy.zeros((256, 256))
    #SigmaTwos = numpy.zeros((256, 256))
    # WE DONT NEED THE SIGMAS
    SigmaZeros = None
    SigmaOnes = None
    SigmaTwos = None

    MuOneHistArray = numpy.zeros((256, 256))
    SigmaBs = numpy.zeros((256, 256))
    
    Mask = numpy.zeros((256, 256), dtype = numpy.bool)
    MaxSigmaB = 0
    THRESH = None
    for I in range(256):
        OmegaZero = CumSumHist[I]
        OmegaZeros[I] = OmegaZero
        if OmegaZero > 0:
            #T = numpy.arange(0, I + 1, 1)
            MuZero = (1 / OmegaZero) * numpy.sum(XTimesHist[0:(I + 1)])
            MuZeros[I] = MuZero
            #XC = ARange[0:(I + 1)] - MuZero
            #SigmaZero = (1.0 / OmegaZero) * numpy.sum((XC * XC) * hist[0:(I + 1)])
            #SigmaZeros[I] = SigmaZero

            for J in range(I + 1, 256):
                OmegaOne = CumSumHist[J] - CumSumHist[I]
                OmegaTwo = CumSumHist[-1] - CumSumHist[J]
                OmegaOnes[I, J] = OmegaOne
                OmegaTwos[I, J] = OmegaTwo
                
                #if I == 6 and J == 255:
                #    print str(OmegaOne) + " " + str(OmegaOne) + " " + str(OmegaTwo)
                if OmegaOne > 0 and OmegaTwo > 0:
                    Mask[I, J] = True
                    #T = numpy.arange(I, J + 1, 1)
                    MuOneHistArray[I, J] = numpy.sum(XTimesHist[I:(J + 1)])
                    MuOne = MuOneHistArray[I, J] / OmegaOne
                    #XC = ARange[I:(J + 1)] - MuOne
                    #SigmaOne = (1.0 / OmegaOne) * numpy.sum((XC * XC) * hist[I:(J + 1)])
                    
                    #T = numpy.arange(J, 256, 1)
                    MuTwo = (1.0 / OmegaTwo) * numpy.sum(XTimesHist[J:])
                    #XC = ARange[J:] - MuTwo
                    #SigmaTwo = (1.0 / OmegaTwo) * numpy.sum((XC * XC) * hist[J:])
                    
                    MuOnes[I, J] = MuOne
                    MuTwos[I, J] = MuTwo
                
                    #SigmaOnes[I, J] = SigmaOne
                    #SigmaTwos[I, J] = SigmaTwo

                    T = (MuZero - MuOne - MuTwo)
                    SigmaB = OmegaZero * OmegaOne * OmegaTwo * T * T
                    SigmaBs[I, J] = SigmaB
                #    print str(MuZero) + " " + str(SigmaZero) + " " + str(MuOne) + " " + str(SigmaOne) + " " + str(MuTwo) + " " + str(SigmaTwo)
                #    if I > 10:
                #        quit()
                    
                    if SigmaB > MaxSigmaB:
                        MaxSigmaB = SigmaB
                        THRESH = numpy.array([I, J])
    if returnWorkingValues:
        return (THRESH, OmegaZeros, MuZeros, SigmaZeros, OmegaOnes, MuOnes, SigmaOnes, OmegaTwos, MuTwos, SigmaTwos, MuOneHistArray, SigmaBs, Mask) 
    else:
        return THRESH


# this routine simply cuts off the voxels with intensities outside the percentages given in PercentageCutoff, then runs the Otsu2 function and produces a segmentation image
def robustOtsu(IMG, PercentageCutoff, NumberClasses=2, maskOutZeros = False):

    T = numpy.double(IMG)

    if maskOutZeros:
        nonZeros = IMG.nonzero()
        T = T[nonZeros]
        
        if not PercentageCutoff is None:
            T = numpy.sort(numpy.ravel(T))
            
            bottomIDX = numpy.floor(T.size * PercentageCutoff[0]) - 1
            topIDX    = numpy.ceil(T.size * PercentageCutoff[1])
            
            T = T[int(bottomIDX):int(topIDX)]
    originalRange = numpy.max(T) - numpy.min(T)
    originalMin = numpy.min(T)

    T = (T - originalMin) / originalRange
    T = numpy.uint8(numpy.round(T * 255))

    if NumberClasses == 2:
        OtsuThresh = otsu1(T)
    elif NumberClasses == 3:
        
        OtsuThresh = otsu2(T)
        #import Otsu2Cython
        #OtsuThresh = Otsu2Cython.otsu2(T)
    
    OtsuThresh = numpy.array(OtsuThresh)
    OtsuThresh = OtsuThresh / 255.0 * originalRange + originalMin
    
    #outIMG = numpy.zeros(IMG.shape, dtype=numpy.uint8)
    
    if OtsuThresh.size == 1:
        #I = numpy.where(IMG > OtsuThresh)
        #T = numpy.zeros(IMG.shape, dtype=numpy.uint8)
        #T[I] = 1
                T = numpy.uint8(IMG > OtsuThresh)
    else:
        T = numpy.digitize(numpy.ravel(IMG), OtsuThresh)
        T = numpy.reshape(T, IMG.shape)
    
    if maskOutZeros:
        T[nonZeros] = T[nonZeros] + 1
        T[IMG == 0] = 0

    return T
