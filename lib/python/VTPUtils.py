import os
import vtk
import numpy

def vertexRegionId(F, FaceRegionId):
    numVertices = numpy.max(F) + 1
    VRegionId = numpy.zeros(numVertices, dtype = numpy.int8)
    (U, I, J) = numpy.unique(FaceRegionId, return_index = True, return_inverse = True)

    VertexFaceIDX = vertexFaceIDX(F)

    for z in range(VRegionId.size):
        H = numpy.bincount(J[VertexFaceIDX[z]], minlength = U.size)
        VRegionId[z] = U[numpy.argmax(H)]
    return VRegionId


def readVTPSurf(inFileName):
    
    if not os.path.isfile(inFileName):
        print("File Not Found: " + inFileName)
        return None
    Reader = vtk.vtkXMLPolyDataReader()
    Reader.SetFileName(inFileName)
    Reader.Update()

    Data = Reader.GetOutput()
    Vrts = Data.GetVerts()
    indices = [Vrts.GetData().GetValue(i) for i in range(1, Vrts.GetSize())]

    S = dict()
    S['vertices'] = [list(Data.GetPoint(point_id)) for point_id in range(Data.GetNumberOfPoints())]
    S['vertices'] = numpy.stack(S['vertices']).T

    if Data.GetNumberOfPolys() > 0:
        S['faces'] = [[int(Data.GetPolys().GetData().GetValue(j)) for j in range(i*4 + 1, i*4 + 4)] for i in range(Data.GetPolys().GetNumberOfCells())]
        S['faces'] = numpy.stack(S['faces']).T
    return S

def readVTPCellArray(inFileName, cellArrayName):
    
    if not os.path.isfile(inFileName):
        print("File Not Found: " + inFileName)
        return None
    Reader = vtk.vtkXMLPolyDataReader()
    Reader.SetFileName(inFileName)
    Reader.Update()

    Data = Reader.GetOutput()
    Vrts = Data.GetVerts()
    indices = [Vrts.GetData().GetValue(i) for i in range(1, Vrts.GetSize())]

    cellData = Data.GetCellData()
    for arrayIDX in range(cellData.GetNumberOfArrays()):
        curArrayName = cellData.GetArrayName(arrayIDX)
        curArrayName = curArrayName.replace(" ", "_")
        if curArrayName == cellArrayName:
            curArray = cellData.GetArray(arrayIDX)
            curArray = [curArray.GetTuple(arrayIDX) for arrayIDX in range(curArray.GetNumberOfTuples())]
            curArray = numpy.array(curArray)
            return curArray
