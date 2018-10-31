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

def getVertexNeighbours(S):
    vertexNeighbours = list()

    for z in range(S['vertices'].shape[1]):
            vertexNeighbours.append(list())

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

def vertexFaceIDX(F):
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

#def connectedComponents(S, Mask):





#if __name__ == "__main__":
#    S = freesurfer.readSurf(os.path.join('bonnie_mirtk', 'freesurfer', 'P01', 'surf', 'lh.white'))
#    G = meshIrJc(S)
