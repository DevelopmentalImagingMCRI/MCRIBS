import os
import vtk
import numpy
import vtk.util.numpy_support

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
    #D = Vrts.GetData()
    #indices = [D.GetValue(i) for i in range(1, Vrts.GetSize())]

    S = dict()
    #S['vertices'] = [list(Data.GetPoint(point_id)) for point_id in range(Data.GetNumberOfPoints())]
    #S['vertices'] = numpy.stack(S['vertices']).T
    S['vertices'] = numpy.array(vtk.util.numpy_support.vtk_to_numpy(Data.GetPoints().GetData())).T
    if Data.GetNumberOfPolys() > 0:
        S['faces'] = numpy.array(vtk.util.numpy_support.vtk_to_numpy(Data.GetPolys().GetData()))
        S['faces'] = numpy.reshape(S['faces'], (int(S['faces'].size / 4), 4)).T
        S['faces'] = S['faces'][1:]
        #
        # D = Data.GetPolys().GetData()
        # S['faces'] = [[int(D.GetValue(j)) for j in range(i*4 + 1, i*4 + 4)] for i in range(Data.GetPolys().GetNumberOfCells())]
        # S['faces'] = numpy.stack(S['faces']).T
    return S

def readVTPPointDataArray(inFileName, pointDataArrayName):

    if not os.path.isfile(inFileName):
        print("File Not Found: " + inFileName)
        return None
    Reader = vtk.vtkXMLPolyDataReader()
    Reader.SetFileName(inFileName)
    Reader.Update()

    Data = Reader.GetOutput()
    #Vrts = Data.GetVerts()
    #indices = [Vrts.GetData().GetValue(i) for i in range(1, Vrts.GetSize())]
    pointData = Data.GetPointData()

    for arrayIDX in range(pointData.GetNumberOfArrays()):
        curArrayName = printData.GetArrayName(arrayIDX)
        curArrayName = curArrayName.replace(" ", "_")
        if curArrayName == cellArrayName:
            curArray = numpy.array(vtk.util.numpy_support.vtk_to_numpy(pointData.GetArray(arrayIDX)))
            # curArray = cellData.GetArray(arrayIDX)
            # curArray = [curArray.GetTuple(arrayIDX) for arrayIDX in range(curArray.GetNumberOfTuples())]
            # curArray = numpy.array(curArray)
            return curArray

def readVTPCellArray(inFileName, cellArrayName):

    if not os.path.isfile(inFileName):
        print("File Not Found: " + inFileName)
        return None
    Reader = vtk.vtkXMLPolyDataReader()
    Reader.SetFileName(inFileName)
    Reader.Update()

    Data = Reader.GetOutput()
    #Vrts = Data.GetVerts()
    #indices = [Vrts.GetData().GetValue(i) for i in range(1, Vrts.GetSize())]

    cellData = Data.GetCellData()
    for arrayIDX in range(cellData.GetNumberOfArrays()):
        curArrayName = cellData.GetArrayName(arrayIDX)
        curArrayName = curArrayName.replace(" ", "_")
        if curArrayName == cellArrayName:
            curArray = numpy.array(vtk.util.numpy_support.vtk_to_numpy(cellData.GetArray(arrayIDX)))
            # curArray = cellData.GetArray(arrayIDX)
            # curArray = [curArray.GetTuple(arrayIDX) for arrayIDX in range(curArray.GetNumberOfTuples())]
            # curArray = numpy.array(curArray)
            return curArray
