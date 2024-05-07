#!/usr/bin/env python

import numpy
import freesurfer

import os

# returns a Ir and Jc vector for a mesh
def meshIrJc(S):
    E = numpy.concatenate((numpy.stack((S['faces'][0], S['faces'][1])), numpy.stack((S['faces'][2], S['faces'][1])), numpy.stack((S['faces'][0], S['faces'][2]))), axis = 1)
    E = numpy.concatenate((E, numpy.take(E, [1, 0], axis = 0)), axis = 1)
    I = numpy.lexsort(E)

    ESorted = numpy.take(E, I, axis = 1)

    # extract unique
    NonUniqueIDX = numpy.concatenate((numpy.where(numpy.any(numpy.diff(ESorted, axis = 1) != 0, axis = 0))[0], numpy.array([ESorted.shape[1] - 1])))

    ESorted = numpy.take(ESorted, NonUniqueIDX, axis = 1)

    IrJc = dict()
    IrJc['Ir'] = numpy.array(ESorted[0])

    H = numpy.bincount(ESorted[1], minlength = S['vertices'].shape[1])
    IrJc['Jc'] = numpy.append(numpy.array([0]), numpy.cumsum(H))
    #print numpy.where(numpy.all(numpy.diff(ESorted, axis = 1) == 0, axis = 0))[0]
    #print S['vertices'].shape[1]
    #print IrJc['Jc'].shape
    #print IrJc['Jc']
    #print IrJc['Ir'].shape

    return IrJc

def getVertexNeighboursOldSlow(S):
    # vertexNeighbours = list()
    #
    # for z in range(S['vertices'].shape[1]):
    #     vertexNeighbours.append(list())

    vertexNeighbours = [[] for z in range(S['vertices'].shape[1])]

    for z in range(S['faces'].shape[1]):
        vertexNeighbours[S['faces'][0, z]].append(S['faces'][1, z])
        vertexNeighbours[S['faces'][1, z]].append(S['faces'][0, z])
        vertexNeighbours[S['faces'][0, z]].append(S['faces'][2, z])
        vertexNeighbours[S['faces'][2, z]].append(S['faces'][0, z])
        vertexNeighbours[S['faces'][2, z]].append(S['faces'][1, z])
        vertexNeighbours[S['faces'][1, z]].append(S['faces'][2, z])

    for z in range(S['vertices'].shape[1]):
        vertexNeighbours[z] = numpy.unique(numpy.array(vertexNeighbours[z]))

    return vertexNeighbours

def getVertexNeighbours(S):

    E = numpy.concatenate((numpy.stack((S['faces'][0], S['faces'][1])), numpy.stack((S['faces'][2], S['faces'][1])), numpy.stack((S['faces'][0], S['faces'][2]))), axis = 1)
    E = numpy.concatenate((E, numpy.take(E, [1, 0], axis = 0)), axis = 1)
    I = numpy.lexsort(E)

    ESorted = numpy.take(E, I, axis = 1)

    # extract unique
    NonUniqueIDX = numpy.concatenate((numpy.where(numpy.any(numpy.diff(ESorted, axis = 1) != 0, axis = 0))[0], numpy.array([ESorted.shape[1] - 1])))

    ESorted = numpy.take(ESorted, NonUniqueIDX, axis = 1)
    #print(ESorted)

    vertexNeighbours = [[] for z in range(S['vertices'].shape[1])]

    for z in range(ESorted.shape[1]):
        vertexNeighbours[ESorted[1, z]].append(ESorted[0, z])

    for z in range(S['vertices'].shape[1]):
        vertexNeighbours[z] = numpy.array(vertexNeighbours[z])
    return vertexNeighbours


def vertexFaceIDX(F):
    """_summary_

    Parameters
    ----------
    F : numpy.ndarray
        Face array, each column is an index vector into the array of vertices.

    Returns
    -------
    list
        Each element f[I] contains a list of face indices that vertex I belongs to.
    """    
    numVertices = numpy.max(F) + 1

    vertexFaceIDX = list()
    for z in range(numVertices):
        vertexFaceIDX.append(list())

    for z in range(F.shape[1]):
        vertexFaceIDX[F[0, z]].append(z)
        vertexFaceIDX[F[1, z]].append(z)
        vertexFaceIDX[F[2, z]].append(z)

    for z in range(numVertices):
        vertexFaceIDX[z] = numpy.array(vertexFaceIDX[z])
    return vertexFaceIDX

# for the surface S returns a tuple (FaceNormals, FaceAreas, VertexNormals, VertexAreas)

def surfaceAreasNormals(S):

    VA = numpy.take(S['vertices'], S['faces'][0], axis = 1)
    VB = numpy.take(S['vertices'], S['faces'][1], axis = 1)
    VC = numpy.take(S['vertices'], S['faces'][2], axis = 1)

    faceNormals = numpy.cross(VB - VA, VC - VA, axis = 0)
    faceAreas = numpy.sqrt(numpy.sum(faceNormals * faceNormals, axis = 0)) / 2
    faceNormals = faceNormals / faceAreas / 2

    vertexFaces = vertexFaceIDX(S['faces'])

    vertexNormals = numpy.zeros_like(S['vertices'])
    vertexAreas = numpy.zeros(S['vertices'].shape[1])

    for z in range(S['vertices'].shape[1]):
        vertexAreas[z] = numpy.sum(faceAreas[vertexFaces[z]]) / 3

        vertexNormals[:, z] = numpy.sum(numpy.atleast_2d(faceAreas[vertexFaces[z]]) * faceNormals[:, vertexFaces[z]], axis = 1)

    vertexNormals = vertexNormals / numpy.atleast_2d(numpy.sqrt(numpy.sum(vertexNormals * vertexNormals, axis = 1))).T
    return (faceNormals, faceAreas, vertexNormals, vertexAreas)

def connectedComponents(S, Mask=None):
    """Label vertices of a surface into connected components.

    Parameters
    ----------
    S : dict
        Surface data, must contain 'vertices' and 'faces'
    Mask : numpy.ndarray, optional
        Only process vertices that are True, by default None

    Returns
    -------
    numpy.ndarray
        Labels for each vertex according to their connected component.
    """    
    Neighbours = getVertexNeighbours(S)

    ConnIDX = numpy.zeros(Mask.size, dtype = numpy.int32)
    CurConnLabel = 1

    while True:
        if Mask is not None:
            Unvisited = numpy.where(numpy.logical_and(ConnIDX == 0, Mask))[0]
        else:
            Unvisited = numpy.where(ConnIDX == 0)[0]

        if Unvisited.size == 0:
            break

        CurVertices = [Unvisited[0]]
        ConnIDX[Unvisited[0]] = CurConnLabel

        while True:
            CurNeighbours = [Neighbours[x] for x in CurVertices]
            CurNeighbours = numpy.unique(numpy.concatenate(CurNeighbours))
            if Mask is not None:
               I = numpy.logical_and(Mask[CurNeighbours], ConnIDX[CurNeighbours] == 0)
            else:
               I = (ConnIDX[CurNeighbours] == 0)

            if numpy.all(I == False):
                break
            ConnIDX[CurNeighbours[I]] = CurConnLabel
            CurVertices = CurNeighbours[I]
        CurConnLabel = CurConnLabel + 1
    return ConnIDX





#if __name__ == "__main__":
#    S = freesurfer.readSurf(os.path.join('bonnie_mirtk', 'freesurfer', 'P01', 'surf', 'lh.white'))
#    G = meshIrJc(S)
