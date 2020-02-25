#!/usr/bin/python3

import numpy

import os
import sys
import struct

import copy

surfaceMagicNumber = b'\xff\xff\xfe'

APARCLabels = {'ctx-lh-bankssts': 'lBSTS',\
'ctx-lh-caudalanteriorcingulate': 'lCAC',\
'ctx-lh-caudalmiddlefrontal': 'lCMF',\
'ctx-lh-cuneus': 'lCUN',\
'ctx-lh-entorhinal': 'lENT',\
'ctx-lh-fusiform': 'lFUS',\
'ctx-lh-inferiorparietal': 'lINFP',\
'ctx-lh-inferiortemporal': 'lIT',\
'ctx-lh-isthmuscingulate': 'lISTC',\
'ctx-lh-lateraloccipital': 'lLOCC',\
'ctx-lh-lateralorbitofrontal': 'lLORB',\
'ctx-lh-lingual': 'lLIN',\
'ctx-lh-medialorbitofrontal': 'lMORB',\
'ctx-lh-middletemporal': 'lMT',\
'ctx-lh-parahippocampal': 'lPARH',\
'ctx-lh-paracentral': 'lPARC',\
'ctx-lh-parsopercularis': 'lPOPE',\
'ctx-lh-parsorbitalis': 'lPORB',\
'ctx-lh-parstriangularis': 'lPTRI',\
'ctx-lh-pericalcarine': 'lPCAL',\
'ctx-lh-postcentral': 'lPSTC',\
'ctx-lh-posteriorcingulate': 'lPC',\
'ctx-lh-precentral': 'lPREC',\
'ctx-lh-precuneus': 'lPCUN',\
'ctx-lh-rostralanteriorcingulate': 'lRAC',\
'ctx-lh-rostralmiddlefrontal': 'lRMF',\
'ctx-lh-superiorfrontal': 'lSF',\
'ctx-lh-superiorparietal': 'lSP',\
'ctx-lh-superiortemporal': 'lST',\
'ctx-lh-supramarginal': 'lSMAR',\
'ctx-lh-frontalpole': 'lFP',\
'ctx-lh-temporalpole': 'lTP',\
'ctx-lh-transversetemporal': 'lTT',\
'ctx-rh-bankssts': 'rBSTS',\
'ctx-rh-caudalanteriorcingulate': 'rCAC',\
'ctx-rh-caudalmiddlefrontal': 'rCMF',\
'ctx-rh-cuneus': 'rCUN',\
'ctx-rh-entorhinal': 'rENT',\
'ctx-rh-fusiform': 'rFUS',\
'ctx-rh-inferiorparietal': 'rINFP',\
'ctx-rh-inferiortemporal': 'rIT',\
'ctx-rh-isthmuscingulate': 'rISTC',\
'ctx-rh-lateraloccipital': 'rLOCC',\
'ctx-rh-lateralorbitofrontal': 'rLORB',\
'ctx-rh-lingual': 'rLIN',\
'ctx-rh-medialorbitofrontal': 'rMORB',\
'ctx-rh-middletemporal': 'rMT',\
'ctx-rh-parahippocampal': 'rPARH',\
'ctx-rh-paracentral': 'rPARC',\
'ctx-rh-parsopercularis': 'rPOPE',\
'ctx-rh-parsorbitalis': 'rPORB',\
'ctx-rh-parstriangularis': 'rPTRI',\
'ctx-rh-pericalcarine': 'rPCAL',\
'ctx-rh-postcentral': 'rPSTC',\
'ctx-rh-posteriorcingulate': 'rPC',\
'ctx-rh-precentral': 'rPREC',\
'ctx-rh-precuneus': 'rPCUN',\
'ctx-rh-rostralanteriorcingulate': 'rRAC',\
'ctx-rh-rostralmiddlefrontal': 'rRMF',\
'ctx-rh-superiorfrontal': 'rSF',\
'ctx-rh-superiorparietal': 'rSP',\
'ctx-rh-superiortemporal': 'rST',\
'ctx-rh-supramarginal': 'rSMAR',\
'ctx-rh-frontalpole': 'rFP',\
'ctx-rh-temporalpole': 'rTP',\
'ctx-rh-transversetemporal': 'rTT'}

APARCLabelsWithInsula = {'ctx-lh-bankssts': 'lBSTS',\
'ctx-lh-caudalanteriorcingulate': 'lCAC',\
'ctx-lh-caudalmiddlefrontal': 'lCMF',\
'ctx-lh-cuneus': 'lCUN',\
'ctx-lh-entorhinal': 'lENT',\
'ctx-lh-fusiform': 'lFUS',\
'ctx-lh-inferiorparietal': 'lINFP',\
'ctx-lh-insula': 'lINS',\
'ctx-lh-inferiortemporal': 'lIT',\
'ctx-lh-isthmuscingulate': 'lISTC',\
'ctx-lh-lateraloccipital': 'lLOCC',\
'ctx-lh-lateralorbitofrontal': 'lLORB',\
'ctx-lh-lingual': 'lLIN',\
'ctx-lh-medialorbitofrontal': 'lMORB',\
'ctx-lh-middletemporal': 'lMT',\
'ctx-lh-parahippocampal': 'lPARH',\
'ctx-lh-paracentral': 'lPARC',\
'ctx-lh-parsopercularis': 'lPOPE',\
'ctx-lh-parsorbitalis': 'lPORB',\
'ctx-lh-parstriangularis': 'lPTRI',\
'ctx-lh-pericalcarine': 'lPCAL',\
'ctx-lh-postcentral': 'lPSTC',\
'ctx-lh-posteriorcingulate': 'lPC',\
'ctx-lh-precentral': 'lPREC',\
'ctx-lh-precuneus': 'lPCUN',\
'ctx-lh-rostralanteriorcingulate': 'lRAC',\
'ctx-lh-rostralmiddlefrontal': 'lRMF',\
'ctx-lh-superiorfrontal': 'lSF',\
'ctx-lh-superiorparietal': 'lSP',\
'ctx-lh-superiortemporal': 'lST',\
'ctx-lh-supramarginal': 'lSMAR',\
'ctx-lh-frontalpole': 'lFP',\
'ctx-lh-temporalpole': 'lTP',\
'ctx-lh-transversetemporal': 'lTT',\
'ctx-rh-bankssts': 'rBSTS',\
'ctx-rh-caudalanteriorcingulate': 'rCAC',\
'ctx-rh-caudalmiddlefrontal': 'rCMF',\
'ctx-rh-cuneus': 'rCUN',\
'ctx-rh-entorhinal': 'rENT',\
'ctx-rh-fusiform': 'rFUS',\
'ctx-rh-inferiorparietal': 'rINFP',\
'ctx-rh-insula': 'rINS',\
'ctx-rh-inferiortemporal': 'rIT',\
'ctx-rh-isthmuscingulate': 'rISTC',\
'ctx-rh-lateraloccipital': 'rLOCC',\
'ctx-rh-lateralorbitofrontal': 'rLORB',\
'ctx-rh-lingual': 'rLIN',\
'ctx-rh-medialorbitofrontal': 'rMORB',\
'ctx-rh-middletemporal': 'rMT',\
'ctx-rh-parahippocampal': 'rPARH',\
'ctx-rh-paracentral': 'rPARC',\
'ctx-rh-parsopercularis': 'rPOPE',\
'ctx-rh-parsorbitalis': 'rPORB',\
'ctx-rh-parstriangularis': 'rPTRI',\
'ctx-rh-pericalcarine': 'rPCAL',\
'ctx-rh-postcentral': 'rPSTC',\
'ctx-rh-posteriorcingulate': 'rPC',\
'ctx-rh-precentral': 'rPREC',\
'ctx-rh-precuneus': 'rPCUN',\
'ctx-rh-rostralanteriorcingulate': 'rRAC',\
'ctx-rh-rostralmiddlefrontal': 'rRMF',\
'ctx-rh-superiorfrontal': 'rSF',\
'ctx-rh-superiorparietal': 'rSP',\
'ctx-rh-superiortemporal': 'rST',\
'ctx-rh-supramarginal': 'rSMAR',\
'ctx-rh-frontalpole': 'rFP',\
'ctx-rh-temporalpole': 'rTP',\
'ctx-rh-transversetemporal': 'rTT'}

DKTLabelsWithInsula = {'ctx-lh-caudalanteriorcingulate': 'lCAC',\
'ctx-lh-caudalmiddlefrontal': 'lCMF',\
'ctx-lh-cuneus': 'lCUN',\
'ctx-lh-entorhinal': 'lENT',\
'ctx-lh-fusiform': 'lFUS',\
'ctx-lh-inferiorparietal': 'lINFP',\
'ctx-lh-insula': 'lINS',\
'ctx-lh-inferiortemporal': 'lIT',\
'ctx-lh-isthmuscingulate': 'lISTC',\
'ctx-lh-lateraloccipital': 'lLOCC',\
'ctx-lh-lateralorbitofrontal': 'lLORB',\
'ctx-lh-lingual': 'lLIN',\
'ctx-lh-medialorbitofrontal': 'lMORB',\
'ctx-lh-middletemporal': 'lMT',\
'ctx-lh-parahippocampal': 'lPARH',\
'ctx-lh-paracentral': 'lPARC',\
'ctx-lh-parsopercularis': 'lPOPE',\
'ctx-lh-parsorbitalis': 'lPORB',\
'ctx-lh-parstriangularis': 'lPTRI',\
'ctx-lh-pericalcarine': 'lPCAL',\
'ctx-lh-postcentral': 'lPSTC',\
'ctx-lh-posteriorcingulate': 'lPC',\
'ctx-lh-precentral': 'lPREC',\
'ctx-lh-precuneus': 'lPCUN',\
'ctx-lh-rostralanteriorcingulate': 'lRAC',\
'ctx-lh-rostralmiddlefrontal': 'lRMF',\
'ctx-lh-superiorfrontal': 'lSF',\
'ctx-lh-superiorparietal': 'lSP',\
'ctx-lh-superiortemporal': 'lST',\
'ctx-lh-supramarginal': 'lSMAR',\
'ctx-lh-transversetemporal': 'lTT',\
'ctx-rh-caudalanteriorcingulate': 'rCAC',\
'ctx-rh-caudalmiddlefrontal': 'rCMF',\
'ctx-rh-cuneus': 'rCUN',\
'ctx-rh-entorhinal': 'rENT',\
'ctx-rh-fusiform': 'rFUS',\
'ctx-rh-inferiorparietal': 'rINFP',\
'ctx-rh-insula': 'rINS',\
'ctx-rh-inferiortemporal': 'rIT',\
'ctx-rh-isthmuscingulate': 'rISTC',\
'ctx-rh-lateraloccipital': 'rLOCC',\
'ctx-rh-lateralorbitofrontal': 'rLORB',\
'ctx-rh-lingual': 'rLIN',\
'ctx-rh-medialorbitofrontal': 'rMORB',\
'ctx-rh-middletemporal': 'rMT',\
'ctx-rh-parahippocampal': 'rPARH',\
'ctx-rh-paracentral': 'rPARC',\
'ctx-rh-parsopercularis': 'rPOPE',\
'ctx-rh-parsorbitalis': 'rPORB',\
'ctx-rh-parstriangularis': 'rPTRI',\
'ctx-rh-pericalcarine': 'rPCAL',\
'ctx-rh-postcentral': 'rPSTC',\
'ctx-rh-posteriorcingulate': 'rPC',\
'ctx-rh-precentral': 'rPREC',\
'ctx-rh-precuneus': 'rPCUN',\
'ctx-rh-rostralanteriorcingulate': 'rRAC',\
'ctx-rh-rostralmiddlefrontal': 'rRMF',\
'ctx-rh-superiorfrontal': 'rSF',\
'ctx-rh-superiorparietal': 'rSP',\
'ctx-rh-superiortemporal': 'rST',\
'ctx-rh-supramarginal': 'rSMAR',\
'ctx-rh-transversetemporal': 'rTT'}


def readSurf(fileName):
    if not os.path.isfile(fileName):
        return None
    else:
        FID = open(fileName, 'rb')

        # read the magic number, should be 0xFFFFFE
        magicNumber = FID.read(3)
        if magicNumber != surfaceMagicNumber:
            # potentially an ascii file, check for "#!"
            FID.seek(0)
            magicNumber = FID.read(2).decode()
            if magicNumber != "#!":
                print("Checked for binary and ASCII versions of Freesurfer file, was neither, invalid file")
                fileFormat = None
            else:
                fileFormat = 'ascii'
        else:
            fileFormat = 'binary'

        if fileFormat == None:
            surfOut = None
        else:
            surfOut = dict()
            if fileFormat == 'binary':
                surfOut['description'] = FID.readline().decode().rstrip()
                fPos = FID.tell()

                T = FID.read(1).decode()
                if ord(T) != 10:
                    #print("Dummy wrong")
                    FID.seek(fPos)

                numVertices = struct.unpack('>I', FID.read(4))[0]
                numFaces = struct.unpack('>I', FID.read(4))[0]

                surfOut['vertices'] = numpy.array(struct.unpack('>' + 'f' * numVertices * 3, FID.read(numVertices * 3 * 4)))
                surfOut['faces'] = numpy.array(struct.unpack('>' + 'i' * numFaces * 3, FID.read(numFaces * 3 * 4)))

                surfOut['vertices'] = numpy.reshape(surfOut['vertices'], [3, numVertices], order = 'F')
                surfOut['faces'] = numpy.reshape(surfOut['faces'], [3, numFaces], order = 'F')
            #    except IOError e:
            #        print("Some io error occurred")
            #        return None

                #print(surfOut['description'])
                #print(numVertices)
                #print(numFaces)
                #print(surfOut['vertices'][:, :10])
                #print(surfOut['faces'][:, :10])
            elif fileFormat == 'ascii':
                #print("ascii mode")
                surfOut['description'] = FID.readline().decode().rstrip()
                #print(surfOut['description'])

                numVertices, numFaces = FID.readline().decode().rstrip().split()

                numVertices = int(numVertices)
                numFaces = int(numFaces)

                #print(numVertices)
                #print(numFaces)

                surfOut['vertices'] = numpy.zeros((3, numVertices), dtype = numpy.float)
                surfOut['faces'] = numpy.zeros((3, numFaces), dtype = numpy.int32)

                for z in range(numVertices):
                    T = FID.readline().decode().rstrip().split()[:3]

                    #T = [float(T[i]) for i in range(len(T))]
                    #surfOut['vertices'][:, z] = numpy.array(T)
                    surfOut['vertices'][0, z] = float(T[0])
                    surfOut['vertices'][1, z] = float(T[1])
                    surfOut['vertices'][2, z] = float(T[2])

                for z in range(numFaces):
                    T = FID.readline().decode().rstrip().split()[:3]
                    #T = [int(T[i]) for i in range(len(T))]
                    #surfOut['faces'][:, z] = numpy.array(T)
                    surfOut['faces'][0, z] = int(T[0])
                    surfOut['faces'][1, z] = int(T[1])
                    surfOut['faces'][2, z] = int(T[2])

                #print(surfOut['vertices'][:, :10])
                #print(surfOut['faces'][:, :10])
        FID.close()
        return surfOut

import subprocess
import tempfile
import shutil

def writeSurf(surfStruct, fileName, fileFormat = 'binary', geometryNIIFile = None, convertToTKR = True):
    if not isinstance(surfStruct, dict):
        raise("surfStruct should be a dict")

    if not 'vertices' in surfStruct or not 'faces' in surfStruct:
        raise("surfStruct must contain 'vertices' and 'faces'")

    if fileFormat == 'binary':
        try:
            FID = open(fileName, 'wb')
        except Exception:
            print(("Could not open: " + fileName))
            return

        FID.write(surfaceMagicNumber)
        if 'description' in surfStruct:
            FID.write((surfStruct['description'].rstrip() + "\n\n").encode())
        else:
            FID.write(("Freesurfer surface\n\n").encode())

        numVertices = int(surfStruct['vertices'].shape[1])
        if surfStruct['faces'] is None:
            numFaces = 0
        else:
            numFaces = int(surfStruct['faces'].shape[1])

        FID.write(struct.pack('>II', numVertices, numFaces))

        blockSize = 1000

        T = surfStruct['vertices'].flatten(order = 'F')

        for leftIDX in range(0, numpy.size(T), blockSize):
            rightIDX = numpy.minimum(leftIDX + blockSize, numpy.size(T))

            S = T[leftIDX:rightIDX].tolist()

            FID.write(struct.pack('>' + 'f' * len(S), *S))
            del S
        del T

        if not surfStruct['faces'] is None:
            T = surfStruct['faces'].flatten(order = 'F')
            for leftIDX in range(0, numpy.size(T), blockSize):
                rightIDX = numpy.minimum(leftIDX + blockSize, numpy.size(T))

                S = T[leftIDX:rightIDX].tolist()

                FID.write(struct.pack('>' + 'I' * len(S), *S))
                del S
            del T

    elif fileFormat == 'ascii':
        try:
            FID = open(fileName, 'w')
        except Exception:
            print(("Could not open: " + fileName))
            return

        if 'description' in surfStruct:
            FID.write(("#!" + surfStruct['description'].rstrip() + "\n"))
        else:
            FID.write("#!Freesurfer surface\n")

        FID.write(("%d %d\n" % (int(surfStruct['vertices'].shape[1]), int(surfStruct['faces'].shape[1]))))

        for z in range(surfStruct['vertices'].shape[1]):
            FID.write(("%f %f %f 1\n" % (surfStruct['vertices'][0, z], surfStruct['vertices'][1, z], surfStruct['vertices'][2, z])))

        for z in range(surfStruct['faces'].shape[1]):
            FID.write(("%d %d %d 1\n" % (surfStruct['faces'][0, z], surfStruct['faces'][1, z], surfStruct['faces'][2, z])))
    FID.close()
    if not geometryNIIFile is None:
        #geometryNII = nibabel.load(geometryNIIFile)
        # /* write whether vertex data was using the
        #        real RAS rather than conformed RAS */
        #     fwriteInt(TAG_OLD_USEREALRAS, fp);
        #FID.write(struct.pack('>I', 2)) # TAG_OLD_USEREALRAS

        #fwriteInt(mris->useRealRAS, fp);
        #FID.write(struct.pack('>I', 1))
        #volume info
#fwriteInt(TAG_OLD_SURF_GEOM, fp);
        #FID.write(struct.pack('>I', 20)) # TAG_OLD_USEREALRAS


        #2846     S = vg_i_to_r(&mris->vg);

        #2847     T = TkrVox2RASfromVolGeom(&mris->vg);
        #2848     Tinv = MatrixInverse(T,NULL);
        #2849     M = MatrixMultiply(S,Tinv,NULL);
        #2850     MRISmatrixMultiply(mris,M);
        #2851     mris->useRealRAS = 1;
        #2852     MatrixFree(&S);
        #2853     MatrixFree(&T);
        #2854     MatrixFree(&Tinv);
        #2855     MatrixFree(&M);


#        FID.write("valid = 1    # volume info valid\n")
#        FID.write("filename = " + geometryNIIFile + "\n")
#        FID.write("volume = " + str(geometryNII.shape[0]) + " " + str(geometryNII.shape[1]) + " " + str(geometryNII.shape[2]) + "\n")
#        #FID.write("voxelsize = " + str(geometryNII.header.get_zooms()[0]) + " " + str(geometryNII.header.get_zooms()[1]) + " " + str(geometryNII.header.get_zooms()[2]) + "\n")
#        A = geometryNII.get_affine()
#        T = numpy.matrix(A) / geometryNII.header.get_zooms()[0]
#        #A[0] = -A[0]
#        #A[1:3] = numpy.take(A, [2, 1], axis = 0)
#
#        #print T
#        # swap the second and third rows
#
#        #T[1:3] = numpy.take(T, [2, 1], axis = 0)
#        #T[1] = -T[1]
#        #print T
#        #T[1, 3] = geometryNII.shape[2] - T[1, 3]
#        #print T
#
#        IMGCentre = numpy.matrix((geometryNII.shape[1], geometryNII.shape[2], geometryNII.shape[0])).T / 2.0
#        IMGCentre = numpy.matrix(geometryNII.shape).T / 2.0
#
#        IMGCentre = numpy.concatenate((IMGCentre, numpy.matrix(1)))
#        CRAS = numpy.matrix(A) * IMGCentre
#
#        FID.write(("xras     = %.15e %.15e %.15e\n" % (T[0, 0], T[0, 1], T[0, 2])))
#        FID.write(("yras     = %.15e %.15e %.15e\n" % (T[1, 0], T[1, 1], T[1, 2])))
#        FID.write(("zras     = %.15e %.15e %.15e\n" % (T[2, 0], T[2, 1], T[2, 2])))
#        FID.write(("cras     = %.15e %.15e %.15e\n" % (CRAS[0, 0], CRAS[1, 0], CRAS[2, 0])))
#
        # use mris_convert to convert the surface to RAS format
        # put orientation information into the original file
        d1 = tempfile.NamedTemporaryFile(delete = False)
        head, tail = os.path.split(d1.name)
        d1OutFile = os.path.join(head, 'rh.' + tail)
        d2 = tempfile.NamedTemporaryFile(delete = False)
        head, tail = os.path.split(d2.name)
        d2OutFile = os.path.join(head, 'rh.' + tail)
        CMD = ['mris_convert', '--vol-geom', geometryNIIFile, fileName, d1.name]
        FNULL = open(os.devnull, 'w')
        subprocess.call(CMD, stdout = FNULL, stderr = subprocess.STDOUT)
        FNULL.close()

        if convertToTKR == True:
            CMD = ['mris_convert', '--vol-geom', geometryNIIFile, '--to-tkr', d1.name, d2.name]
            #print " ".join(CMD)

            FNULL = open(os.devnull, 'w')
            subprocess.call(CMD, stdout = FNULL, stderr = subprocess.STDOUT)
            FNULL.close()
            shutil.copyfile(d2.name, fileName)
            os.unlink(d1.name)
            os.unlink(d2.name)
        else:
            shutil.copyfile(d1.name, fileName)
            os.unlink(d1.name)

# 24-bit -1
curvMagicNumber = b'\xff\xff\xff'

def readCurv(fileName):
    if not os.path.isfile(fileName):
        return None
    else:
        FID = open(fileName, 'rb')

        magicNumber = FID.read(3)

        if magicNumber != curvMagicNumber:
            print("Wrong magic number for curv file")
            return None
        else:
            valuesOut = dict()
            valuesOut['numVertices'] = struct.unpack('>I', FID.read(4))[0]
            valuesOut['numFaces'] = struct.unpack('>I', FID.read(4))[0]

            numValuesPerVertex = struct.unpack('>I', FID.read(4))[0]
            valuesOut['values'] = numpy.array(struct.unpack('>' + 'f' * numValuesPerVertex * valuesOut['numVertices'], FID.read(numValuesPerVertex * valuesOut['numVertices'] * 4)))
            FID.close()
            return valuesOut

def writeLabel(values, fileName):
    if not isinstance(values, dict):
        raise("values should be a dict")

    if not 'index' in values or not 'RAS' in values:
        raise("values must contain 'index' and 'RAS' at minimum")

    if values['RAS'].shape[0] != 3:
        raise ValueError('RAS must have 3 rows')

    if values['index'].size != values['RAS'].shape[1]:
        raise ValueError('index and RAS have inconsistent sizes')

    try:
        FID = open(fileName, 'w')
    except Exception:
        print("Could not open: " + fileName)
        return

    if not 'comment' in values:
        FID.write("comment\n")
    else:
        FID.write("%s\n" % (values['comment'].rstrip()))

    if not 'numElements' in values:
        numElements = values['index'].size
    else:
        numElements = values['numElements']

    FID.write("%d\n" % (numElements))

    for z in range(values['index'].size):
        FID.write("%d\t%f\t%f\t%f\t%f\n" % (values['index'][z], values['RAS'][0, z], values['RAS'][1, z], values['RAS'][2, z], 0))

    FID.close()



def readLabel(fileName):

    if not os.path.isfile(fileName):
        return None
    else:
        FID = open(fileName, 'r')

        outL = dict()
        outL['comment'] = FID.readline()
        outL['numElements'] = int(FID.readline())

        outL['RAS'] = numpy.zeros((3, outL['numElements']))
        outL['index'] = numpy.zeros(outL['numElements'], dtype = numpy.int32)

        for z in range(outL['numElements']):
            T = FID.readline()
            S = T.split()
            outL['index'][z] = int(S[0])
            outL['RAS'][:, z] = numpy.array([float(S[1]), float(S[2]), float(S[3])])
        FID.close()
        return outL

def writeCurv(values, fileName):
    """
    Writes data to a freesurfer curv format file.

    Parameters
    ----------
    values: dict
        Curvature data.
    values.numVertices: int
        Number of vertices.
    values.numFaces: int
        Number of faces
    values.values: numpy.ndarray
        Vector of values
    fileName: str
        File to save to.
    """
    if not isinstance(values, dict):
        raise("values should be a dict")

    if not 'numVertices' in values or not 'numFaces' in values or not 'values' in values:
        raise("values must contain 'numVertices' and 'numFaces' and 'values'")

    try:
        FID = open(fileName, 'wb')
    except Exception:
        print(("Could not open: " + fileName))
        return

    FID.write(curvMagicNumber)

    FID.write(struct.pack('>II', int(values['numVertices']), int(values['numFaces'])))
    numValuesPerVertex = int(numpy.size(values['values']) / values['numVertices'])

    FID.write(struct.pack('>I', numValuesPerVertex))
    # make sure the values are of type single
    T = numpy.single(values['values']).tolist()

    FID.write(struct.pack('>' + 'f' * numpy.size(T), *T))
    FID.close()

# W files are deprecated, we dont need this
def readW(fileName):
    if not os.path.isfile(fileName):
        return None
    else:
        FID = open(fileName, 'rb')
        FID.close()

# emulates [~, IDX] = ismember(A, B) of matlab, but only returns the index (IDX)

def ismember(A, B):
    flatA = A.flatten()
    flatB = B.flatten()
    sortedAIDX = numpy.argsort(flatA)
    sortedBIDX = numpy.argsort(flatB)

    outIDX = numpy.zeros((numpy.size(flatA)), dtype = numpy.int64)

    outIDX.fill(-1)

    curAIDX = 0
    curBIDX = 0

    while curAIDX < numpy.size(sortedAIDX) and curBIDX < numpy.size(sortedBIDX):
        if flatA[sortedAIDX[curAIDX]] == flatB[sortedBIDX[curBIDX]]:
            outIDX[sortedAIDX[curAIDX]] = sortedBIDX[curBIDX]
            curAIDX += 1
        elif flatA[sortedAIDX[curAIDX]] < flatB[sortedBIDX[curBIDX]]:
            curAIDX += 1
        else:
            curBIDX += 1
    #outIDX = outIDX
    outIDX = numpy.reshape(outIDX, A.shape)
    return outIDX
# reads an annot file

def readAnnot(fileName):
    if not os.path.isfile(fileName):
        annotOut = None
    else:
        annotOut = dict()

        FID = open(fileName, 'rb')

        numElements = struct.unpack('>i', FID.read(4))[0]

        T = struct.unpack('>' + 'i' * numElements * 2, FID.read(numElements * 2 * 4))
        T = numpy.reshape(numpy.array(T), (2, numElements), order = 'F')

        annotOut['vertices'] = T[0]
        annotOut['label'] = T[1]

        #print(annotOut['vertices'].shape)
        #print(annotOut['label'].shape)

        del T
        isThereAColorTable = FID.read(4)

        if len(isThereAColorTable) == 0:
            annotOut['colortable'] = None
        else:
            isThereAColorTable = struct.unpack('>I', isThereAColorTable)[0]

            if isThereAColorTable == 0:
                annotOut['colortable'] = None
            else:
                numEntries = struct.unpack('>i', FID.read(4))[0]

                if numEntries > 0:
                    # original version
                    annotOut['colortable'] = dict()
                    annotOut['colortable']['numEntries'] = numEntries

                    L = struct.unpack('>i', FID.read(4))[0]
                    annotOut['colortable']['orig_tab'] = FID.read(L).decode()
                    for z in range(numEntries):
                        L = struct.unpack('>i', FID.read(4))[0]
                        annotOut['colortable']['struct_names'][structureIDX] = FID.read(L).decode()[:-1]

                        structureIndices = struct.unpack('>iiii', FID.read(4 * 4))
                        annotOut['colortable']['labels'][structureIDX] = (structureIndices[0] & 0xFF) | ((structureIndices[1] & 0xFF) << 8)| ((structureIndices[2] & 0xFF) << 16) | ((structureIndices[3] & 0xFF) << 24)
                else:
                    version = -numEntries
                    if version != 2:
                        annotOut['colortable'] = None
                    else:
                        annotOut['colortable'] = dict()

                        annotOut['colortable']['numEntries'] = struct.unpack('>i', FID.read(4))[0]
                        L = struct.unpack('>i', FID.read(4))[0]
                        annotOut['colortable']['orig_tab'] = FID.read(L).decode()[:-1]
                        del L

                        annotOut['colortable']['struct_names'] = [None] * annotOut['colortable']['numEntries']
                        annotOut['colortable']['labels'] = numpy.zeros((annotOut['colortable']['numEntries']), dtype = numpy.int32)
                        annotOut['colortable']['table'] = numpy.zeros((annotOut['colortable']['numEntries'], 4), dtype = numpy.int32)

                        numEntriesToRead = struct.unpack('>i', FID.read(4))[0]
                        for z in range(numEntriesToRead):
                            structureIDX = struct.unpack('>i', FID.read(4))[0]

                            if structureIDX < 0:
                                print("Warning: invalid structure index")
                                L = struct.unpack('>i', FID.read(4))[0]
                                FID.read(L)
                                FID.read(4 * 4)
                            else:
                                if annotOut['colortable']['struct_names'][structureIDX] != None:
                                    print("Warning: duplicate structure, overwriting")
                                L = struct.unpack('>i', FID.read(4))[0]
                                annotOut['colortable']['struct_names'][structureIDX] = FID.read(L).decode()[:-1]

                                structureIndices = struct.unpack('>iiii', FID.read(4 * 4))
                                annotOut['colortable']['labels'][structureIDX] = (structureIndices[0] & 0xFF) | ((structureIndices[1] & 0xFF) << 8)| ((structureIndices[2] & 0xFF) << 16) | ((structureIndices[3] & 0xFF) << 24)
                                annotOut['colortable']['table'][structureIDX] = numpy.array(structureIndices)

                #if numEntries > 0:
                if annotOut['colortable'] == None:
                    annotOut['seg'] = None
                else:
                    #print(annotOut['colortable']['labels'])
                    #print(annotOut['label'])
                    #annotOut['seg'] = ismember(annotOut['label'], annotOut['colortable']['labels'])
                    IDX = numpy.argsort(annotOut['colortable']['labels'])
                    T = numpy.searchsorted(annotOut['colortable']['labels'][IDX], annotOut['label'])
                    annotOut['seg'] = IDX[T]
                    #rint(annotOut['seg'][:30])
                    pass
            #if isThereAColorTable == 0:
        #if len(isThereAColorTable) == 0:
        #print(annotOut['colortable']['struct_names'])
        FID.close()
    #if not os.path.isfile(fileName):

    return annotOut

def writeAnnot(annotDict, fileName):
    if not isinstance(annotDict, dict):
        print("not a dict")
        return None

#fp = fopen(filename, 'w', 'b');
    try:
        FID = open(fileName, 'wb')
    except Exception:
        print(("Could not open: " + fileName))
        return
    FID.write(struct.pack('>I', annotDict['vertices'].size))

# write the vertices and labels interleaved
#fwrite(fp, [vertices(:), label(:)]', 'int');
    for z in range(annotDict['vertices'].size):
        FID.write(struct.pack('>II', annotDict['vertices'][z], annotDict['label'][z]))

    if not "colortable" in annotDict:
        FID.write(0)
    else:
#fwrite(fp, 1, 'int');
        FID.write(struct.pack('>I', 1))
#fwrite(fp, -2, 'int');
        FID.write(struct.pack('>i', int(-2)))

        #colortable.table(z,1) + colortable.table(z,2)*2^8 + colortable.table(z,3)*2^16 + colortable.table(z,4)*2^24;
#% writing version 2 tables only
#fwrite(fp, colortable.numEntries, 'int');
#fwrite(fp, length(colortable.orig_tab) + 1, 'int');
        FID.write(struct.pack('>I', annotDict['colortable']['numEntries']))
#fwrite(fp, [colortable.orig_tab 0], '*char');
        FID.write(struct.pack('>I', len(annotDict['colortable']['orig_tab']) + 1))
        FID.write(annotDict['colortable']['orig_tab'].encode())
        FID.write(struct.pack('>B', 0))
        #print len(annotDict['colortable']['struct_names'])
        FID.write(struct.pack('>I', len(annotDict['colortable']['struct_names'])))
        for z in range(len(annotDict['colortable']['struct_names'])):
#fwrite(fp, structure(z), 'int');
            FID.write(struct.pack('>I', z))
#fwrite(fp, length(colortable.struct_names{z}) + 1, 'int');
            if annotDict['colortable']['struct_names'][z] is None:
                T = "None"
                FID.write(struct.pack('>I', len(T) + 1))
                FID.write(T.encode())
            else:
                FID.write(struct.pack('>I', len(annotDict['colortable']['struct_names'][z]) + 1))
                FID.write(annotDict['colortable']['struct_names'][z].encode())
#fwrite(fp, [colortable.struct_names{z} 0], '*char');
            #print annotDict['colortable']['struct_names'][z]
            FID.write(struct.pack('>B', 0)) # null terminator
#fwrite(fp, colortable.table(z, 1:4), 'int');
            FID.write(struct.pack('>IIII', annotDict['colortable']['table'][z, 0], annotDict['colortable']['table'][z, 1], annotDict['colortable']['table'][z, 2], 0))

#end
#end
#
#fclose(fp);
    FID.close()

# reads the ICO file in the freesurfer home directory
def readICO(order):
    ICOFile = os.path.join(os.environ['FREESURFER_HOME'], 'lib', 'bem', 'ic' + str(order) + '.tri')
    if not os.path.isfile(ICOFile):
        print(("Couldnt find ICO file: " + ICOFile))
        quit()

    FID = open(ICOFile, 'r')

    curLine = FID.readline()

    nvertices = int(curLine)

    surface = dict()
    surface['vertices'] = numpy.zeros((nvertices, 3))
    curLine = 1

    n = 0
    while True:
        curLine = FID.readline()
        if curLine is None:
            break

        s = curLine.split()
        if len(s) == 4:
            vno = int(s[0]) - 1
            surface['vertices'][vno] = numpy.array([float(s[1]), float(s[2]), float(s[3])])
        elif len(s) == 6:
            vno = n
            surface['vertices'][vno] = numpy.array([float(s[0]), float(s[1]), float(s[2])])
        else:
            break
        n = n + 1
        if n >= nvertices:
            break

    nfacesLine = FID.readline()
    nfaces = int(nfacesLine)

    surface['faces'] = numpy.zeros((nfaces, 3), dtype = numpy.int32)
    n = 0
    while True:
        curLine = FID.readline()
        if curLine is None:
            break

        s = curLine.split()
        if len(s) == 3:
            vno3 = int(s[2])
            vno2 = int(s[1])
            vno1 = int(s[0])
            fno = n
        elif len(s) == 4:
            vno3 = int(s[3])
            vno2 = int(s[2])
            vno1 = int(s[1])
            fno = int(s[0])
        surface['faces'][fno - 1] = numpy.array([vno1 - 1, vno2 - 1, vno3 - 1])
        n = n + 1
        if n >= nfaces:
            break
    surface['vertices'] = surface['vertices'].T
    surface['faces'] = surface['faces'].T

    FID.close()

    return surface
