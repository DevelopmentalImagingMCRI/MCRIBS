import os
import vtk
import numpy
import vtk.util.numpy_support

# def vertexRegionId(F, FaceRegionId):
#     numVertices = numpy.max(F) + 1
#     VRegionId = numpy.zeros(numVertices, dtype = numpy.int8)
#     (U, I, J) = numpy.unique(FaceRegionId, return_index=True, return_inverse=True)

#     VertexFaceIDX = vertexFaceIDX(F)

#     for z in range(VRegionId.size):
#         H = numpy.bincount(J[VertexFaceIDX[z]], minlength = U.size)
#         VRegionId[z] = U[numpy.argmax(H)]
#     return VRegionId


def readVTPAll(inFileName: str) -> dict:
    """Read a VTP surface file along with all cell and point data arrays

    Parameters
    ----------
    inFileName : str
        VTP file name.

    Returns
    -------
    dict
        `vertices` and `faces`.
    """
    if not os.path.isfile(inFileName):
        print("File Not Found: " + inFileName)
        return None
    Reader = vtk.vtkXMLPolyDataReader()
    Reader.SetFileName(inFileName)
    Reader.Update()

    Data = Reader.GetOutput()
    Vrts = Data.GetVerts()

    S = dict()
    S['vertices'] = numpy.array(vtk.util.numpy_support.vtk_to_numpy(Data.GetPoints().GetData())).T
    if Data.GetNumberOfPolys() > 0:
        S['faces'] = numpy.array(vtk.util.numpy_support.vtk_to_numpy(Data.GetPolys().GetData()))
        
        if S['faces'].size == Data.GetNumberOfPolys() * 5:
            S['faces'] = numpy.reshape(S['faces'], (int(S['faces'].size / 5), 5)).T
            S['faces'] = S['faces'][2:]
        else:
            S['faces'] = numpy.reshape(S['faces'], (int(S['faces'].size / 4), 4)).T
            S['faces'] = S['faces'][1:]
    
    pointData = Data.GetPointData()
    S['pointData'] = dict()

    for arrayIDX in range(pointData.GetNumberOfArrays()):
        curArrayName = pointData.GetArrayName(arrayIDX)
        curArrayName = curArrayName.replace(" ", "_")
        S['pointData'][curArrayName] = numpy.array(vtk.util.numpy_support.vtk_to_numpy(pointData.GetArray(arrayIDX)))
    
    cellData = Data.GetCellData()
    S['cellData'] = dict()

    for arrayIDX in range(cellData.GetNumberOfArrays()):
        curArrayName = cellData.GetArrayName(arrayIDX)
        curArrayName = curArrayName.replace(" ", "_")
        S['cellData'][curArrayName] = numpy.array(vtk.util.numpy_support.vtk_to_numpy(cellData.GetArray(arrayIDX)))
    return S

def readVTPSurf(inFileName: str) -> dict:
    """Read a VTP surface file.

    Parameters
    ----------
    inFileName : str
        VTP file name.

    Returns
    -------
    dict
        `vertices` and `faces`.
    """
    if not os.path.isfile(inFileName):
        print("File Not Found: " + inFileName)
        return None
    Reader = vtk.vtkXMLPolyDataReader()
    Reader.SetFileName(inFileName)
    Reader.Update()

    Data = Reader.GetOutput()
    Vrts = Data.GetVerts()

    S = dict()
    S['vertices'] = numpy.array(vtk.util.numpy_support.vtk_to_numpy(Data.GetPoints().GetData())).T
    if Data.GetNumberOfPolys() > 0:
        S['faces'] = numpy.array(vtk.util.numpy_support.vtk_to_numpy(Data.GetPolys().GetData()))
        
        if S['faces'].size == Data.GetNumberOfPolys() * 5:
            S['faces'] = numpy.reshape(S['faces'], (int(S['faces'].size / 5), 5)).T
            S['faces'] = S['faces'][2:]
        else:
            S['faces'] = numpy.reshape(S['faces'], (int(S['faces'].size / 4), 4)).T
            S['faces'] = S['faces'][1:]
    return S


def readVTPPointDataArray(inFileName: str, pointDataArrayName: str) -> numpy.array:
    """Read a PointData array from a VTP file.

    Parameters
    ----------
    inFileName : str
        File name of the VTP file.
    pointDataArrayName : str
        Name of the point data array to extract.

    Returns
    -------
    numpy.array
        The PointData array `pointDataArrayName`, None if not found.
    """
    if not os.path.isfile(inFileName):
        print("File Not Found: " + inFileName)
        return None
    Reader = vtk.vtkXMLPolyDataReader()
    Reader.SetFileName(inFileName)
    Reader.Update()

    Data = Reader.GetOutput()
    pointData = Data.GetPointData()

    for arrayIDX in range(pointData.GetNumberOfArrays()):
        curArrayName = pointData.GetArrayName(arrayIDX)
        curArrayName = curArrayName.replace(" ", "_")
        if curArrayName == pointDataArrayName:
            curArray = numpy.array(vtk.util.numpy_support.vtk_to_numpy(pointData.GetArray(arrayIDX)))
            return curArray
    return None


def readVTPCellArray(inFileName: str, cellDataArrayName: str) -> numpy.array:
    """Read a Cell array from a VTP file.

    Parameters
    ----------
    inFileName : str
        File name of the VTP file.
    cellDataArrayName : str
        Name of the cell data array to extract.

    Returns
    -------
    numpy.array
        The CellData array `pointDataArrayName`, None if not found.
    """
    if not os.path.isfile(inFileName):
        print("File Not Found: " + inFileName)
        return None
    Reader = vtk.vtkXMLPolyDataReader()
    Reader.SetFileName(inFileName)
    Reader.Update()

    Data = Reader.GetOutput()
    cellData = Data.GetCellData()

    for arrayIDX in range(cellData.GetNumberOfArrays()):
        curArrayName = cellData.GetArrayName(arrayIDX)
        curArrayName = curArrayName.replace(" ", "_")
        if curArrayName == cellDataArrayName:
            curArray = numpy.array(vtk.util.numpy_support.vtk_to_numpy(cellData.GetArray(arrayIDX)))
            return curArray
    return None
