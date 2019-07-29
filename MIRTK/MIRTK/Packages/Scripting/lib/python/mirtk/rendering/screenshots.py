##############################################################################
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2016 Imperial College London
# Copyright 2016 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############################################################################

"""Python module for the rendering of screenshots."""

import os
import vtk


default_colors = [
    (1, 0, 0),
    (0, 0, 1),
    (0, 1, 0),
    (0, 1, 1),
    (1, 0, 1),
    (1, 1, 0),
    (1, 1, 1)
]


def deep_copy(obj):
    """Make deep copy of VTK object."""
    copy = obj.NewInstance()
    copy.DeepCopy(obj)
    return copy


def iround(x):
    """Round floating point number and cast to int."""
    return int(round(x))


def nearest_voxel(index):
    """Get indices of nearest voxel."""
    i = iround(index[0])
    j = iround(index[1])
    k = iround(index[2])
    return (i, j, k)


def invert_matrix(m):
    """Get inverse of a vtkMatrix4x4."""
    inv = deep_copy(m)
    inv.Invert()
    return inv


def index_to_point(index, origin, spacing):
    """Transform voxel indices to image data point coordinates."""
    x = origin[0] + index[0] * spacing[0]
    y = origin[1] + index[1] * spacing[1]
    z = origin[2] + index[2] * spacing[2]
    return (x, y, z)


def point_to_index(point, origin, spacing):
    """Transform image data point coordinates to voxel."""
    i = (point[0] - origin[0]) / spacing[0]
    j = (point[1] - origin[1]) / spacing[1]
    k = (point[2] - origin[2]) / spacing[2]
    return (i, j, k)


def matrix_to_affine(matrix):
    """Convert vtkMatrix4x4 to NiBabel 'affine' 2D array."""
    return [[matrix.GetElement(0, 0), matrix.GetElement(0, 1),
             matrix.GetElement(0, 2), matrix.GetElement(0, 3)],
            [matrix.GetElement(1, 0), matrix.GetElement(1, 1),
             matrix.GetElement(1, 2), matrix.GetElement(1, 3)],
            [matrix.GetElement(2, 0), matrix.GetElement(2, 1),
             matrix.GetElement(2, 2), matrix.GetElement(2, 3)],
            [matrix.GetElement(3, 0), matrix.GetElement(3, 1),
             matrix.GetElement(3, 2), matrix.GetElement(3, 3)]]


def range_to_level_window(min_value, max_value):
    """Convert min/max value range to level/window parameters."""
    window = max_value - min_value
    level = min_value + .5 * window
    return (level, window)


def auto_image_range(image, percentiles=(1, 99)):
    """Compute range for color transfer function."""
    stats = vtk.vtkImageHistogramStatistics()
    stats.SetInputData(image)
    stats.AutomaticBinningOn()
    stats.SetMaximumNumberOfBins(512)
    stats.SetAutoRangePercentiles(percentiles)
    stats.UpdateWholeExtent()
    return tuple(stats.GetAutoRange())


def auto_level_window(image, percentiles=(1, 99)):
    """Compute level/window for color transfer function."""
    return range_to_level_window(*auto_image_range(image, percentiles))


def add_contour(renderer, plane, polydata,
                transform=None, line_width=3, color=(1, 0, 0)):
    """Add contour of mesh cut by given image plane to render scene."""
    if transform:
        transformer = vtk.vtkTransformPolyDataFilter()
        transformer.SetInputData(polydata)
        transformer.SetTransform(transform)
        transformer.Update()
        polydata = deep_copy(transformer.GetOutput())
        transformer = None
    cutter = vtk.vtkCutter()
    cutter.SetInputData(polydata)
    cutter.SetCutFunction(plane)
    cutter.Update()
    contour = deep_copy(cutter.GetOutput())
    cutter = None

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(contour)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.LightingOff()
    prop.SetRepresentationToWireframe()
    prop.SetLineWidth(line_width)

    if color:
        prop.SetColor(color)
        mapper.ScalarVisibilityOff()
    elif polydata.GetPointData().GetScalars():
        mapper.SetScalarModeToUsePointData()
        mapper.ScalarVisibilityOn()
    elif polydata.GetCellData().GetScalars():
        mapper.SetScalarModeToUseCellData()
        mapper.ScalarVisibilityOn()

    renderer.AddActor(actor)
    return actor


def slice_axes(zdir):
    """Get volume dimensions corresponding to x and z axes of slice."""
    if zdir == 0:
        return (1, 2)
    elif zdir == 1:
        return (0, 2)
    elif zdir == 2:
        return (0, 1)
    else:
        raise Exception("Invalid zdir argument: " + zdir)


def cropping_region(center, axes, width, height):
    bounds = [center[0], center[0], center[1], center[1], center[2], center[2]]
    i = 2 * axes[0]
    j = 2 * axes[1]
    bounds[i] -= (width - 1) / 2
    bounds[i + 1] += (width + 1) / 2
    bounds[j] -= (height - 1) / 2
    bounds[j + 1] += (height + 1) / 2
    return bounds


def slice_view(image, index, width, height, zdir=2, polydata=[], colors=default_colors,
               transform=None, line_width=3, level_window=None, image_lut=None, interpolation="nearest"):
    """Return vtkRenderer for orthogonal image slice."""

    # determine orientation of medical volume
    flip = [False, False, False]
    if transform:
        transform.Update()
        matrix = deep_copy(transform.GetMatrix())
        try:
            from nibabel import aff2axcodes
            codes = aff2axcodes(matrix_to_affine(matrix))
        except Exception:
            codes = ('L', 'A', 'S')
            if matrix.GetElement(0, 0) < 0:
                codes = 'R'
            if matrix.GetElement(1, 1) < 0:
                codes = 'P'
            if matrix.GetElement(2, 2) < 0:
                codes = 'I'
        if codes[0] == 'R':
            flip[0] = True
        if codes[1] == 'P':
            flip[1] = True
        if codes[2] == 'I':
            flip[2] = True

    dims = image.GetDimensions()
    axes = slice_axes(zdir)
    if width < 1:
        width = dims[axes[0]]
    if height < 1:
        height = dims[axes[1]]
    size = [1, 1, 1]
    size[axes[0]] = width
    size[axes[1]] = height
    if zdir == 2:
        up = (0, 1, 0)
    else:
        up = (0, 0, 1)
    spacing = image.GetSpacing()
    distance = 10. * spacing[zdir]
    focal_point = index_to_point(index, image.GetOrigin(), spacing)
    position = list(focal_point)
    if flip[zdir]:
        position[zdir] = position[zdir] - distance
    else:
        position[zdir] = position[zdir] + distance

    margin = 2
    extent = cropping_region(index, axes, size[axes[0]] + margin, size[axes[1]] + margin)

    if flip[0] or flip[1] or flip[2]:
        flip_transform = vtk.vtkTransform()
        flip_transform.Translate(+focal_point[0], +focal_point[1], +focal_point[2])
        flip_transform.Scale(-1. if flip[0] else 1.,
                             -1. if flip[1] else 1.,
                             -1. if flip[2] else 1.)
        flip_transform.Translate(-focal_point[0], -focal_point[1], -focal_point[2])
        points_transform = vtk.vtkTransform()
        points_transform.SetMatrix(matrix)
        points_transform.PostMultiply()
        points_transform.Concatenate(flip_transform)
    else:
        flip_transform = None
        points_transform = None

    mapper = vtk.vtkImageSliceMapper()
    mapper.SetInputData(image)
    mapper.SetOrientation(zdir)
    mapper.SetSliceNumber(extent[2 * zdir])
    mapper.SetCroppingRegion(extent)
    mapper.CroppingOn()
    mapper.Update()

    actor = vtk.vtkImageSlice()
    actor.SetMapper(mapper)
    if flip_transform:
        actor.SetUserTransform(flip_transform)
    prop = actor.GetProperty()
    interpolation = interpolation.lower()
    if interpolation in ("nn", "nearest"):
        prop.SetInterpolationTypeToNearest()
    elif interpolation == "linear":
        prop.SetInterpolationTypeToLinear()
    elif interpolation == "cubic":
        prop.SetInterpolationTypeToCubic()
    else:
        raise ValueError("Invalid interpolation mode: {}".format(interpolation))

    if not level_window:
        level_window = auto_level_window(image)
    prop.SetColorLevel(level_window[0])
    prop.SetColorWindow(level_window[1])
    if image_lut:
        prop.SetLookupTable(image_lut)
        prop.UseLookupTableScalarRangeOn()

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)

    camera = renderer.GetActiveCamera()
    camera.SetViewUp(up)
    camera.SetPosition(position)
    camera.SetFocalPoint(focal_point)
    camera.SetParallelScale(.5 * max((size[axes[0]] - 1) * spacing[axes[0]],
                                     (size[axes[1]] - 1) * spacing[axes[1]]))
    camera.SetClippingRange(distance - .5 * spacing[zdir],
                            distance + .5 * spacing[zdir])
    camera.ParallelProjectionOn()

    # add contours of polygonal data intersected by slice plane
    if isinstance(polydata, vtk.vtkPolyData):
        polydata = [polydata]
    if isinstance(line_width, int):
        line_width = [line_width]
    if isinstance(colors[0], float):
        colors = [colors]
    for i in xrange(len(polydata)):
        if i < len(colors):
            color = colors[i]
        else:
            color = colors[-1]
        if i < len(line_width):
            width = line_width[i]
        else:
            width = line_width[-1]
        add_contour(renderer, plane=mapper.GetSlicePlane(),
                    polydata=polydata[i], transform=points_transform,
                    line_width=width, color=color)
    return renderer


def take_screenshot(window, path=None):
    """Takes vtkRenderWindow instance and writes a screenshot of the rendering.

    window : vtkRenderWindow
        The render window from which a screenshot is taken.
    path : str
        File name path of output PNG file.
        A .png file name extension is appended if missing.

    """
    _offscreen = window.GetOffScreenRendering()
    window.OffScreenRenderingOn()
    window_to_image = vtk.vtkWindowToImageFilter()
    window_to_image.SetInput(window)
    window_to_image.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetInputConnection(window_to_image.GetOutputPort())
    if path:
        if os.path.splitext(path)[1].lower() != '.png':
            path += '.png'
        writer.SetFileName(path)
    else:
        writer.WriteToMemoryOn()
    writer.Write()
    window.SetOffScreenRendering(_offscreen)
    if writer.GetWriteToMemory():
        from IPython.display import Image
        data = str(buffer(writer.GetResult()))
        return Image(data)


def take_orthogonal_screenshots(image, center=None, length=0, offsets=None,
                                size=(512, 512), prefix='screenshot_{n}',
                                suffix=('axial', 'coronal', 'sagittal'),
                                path_format='{prefix}_{suffix}',
                                level_window=None, qform=None,
                                polydata=[], colors=[], line_width=3,
                                trim=False, overwrite=True):
    """Take three orthogonal screenshots of the given image patch.

    Arguments
    ---------

    image : vtkImageData
        Volume data.
    center : (float, float, float) or (int, int, int)
        vtkImageData coordinates or voxel indices (3-tuple of ints)
        of patch center point. Be sure to pass tuple with correct type.
        When not specified, the image center point is used.
    length : float or int
        Side length of patch either in mm (float) or number of voxels (int).
    offsets : list of float or int, optional
        One or more offsets in either mm (float) or number of voxels (int)
        from the center point along the orthogonal viewing direction.
        A screenshot of the volume of interest is taken for each combination
        of offset and orthgonal viewing direction, i.e., the number of
        screenshots is `3 * len(offset)`. For example, when the volume of
        interest is the entire image, i.e., `center=None` and `length=0`,
        this can be used to render multiple slices with one call of this function.
    size : (int, int)
        Either int or 2-tuple/-list of int values specifying the width and height of the screenshot.
    prefix : str
        Common output path prefix of screenshot files.
    suffix : (str, str, str)
        List or tuple of three strings used as suffix of the screenshot
        take from the respective orthogonal view. The order is:
        axial, coronal, sagittal.
    path_format: str
        A format string used to construct the file name of each individual screenshot.
        Make sure to use a suitable format string such that every screenshot has a
        unique file name. The allowed placeholders are:
        - `{prefix}`: Substituted by the `prefix` argument.
        - `{suffix}`: Substituted by the `suffix` argument corresponding to the orthogonal view.
        - `{n}`: The 1-based index of the screenshot.
        - `{i}`: The i (int) volume index of the screenshot center point before trimming.
        - `{j}`: The j (int) volume index of the screenshot center point before trimming.
        - `{k}`: The k (int) volume index of the screenshot center point before trimming.
        - `{x}`: The x (float) vtkImageData coordinates of the center point before trimming.
        - `{y}`: The y (float) vtkImageData coordinates of the center point before trimming.
        - `{z}`: The z (float) vtkImageData coordinates of the center point before trimming.
        - `{vi}`: The i (int) volume index of the viewport center point after trimming.
        - `{vj}`: The j (int) volume index of the viewport center point after trimming.
        - `{vk}`: The k (int) volume index of the viewport center point after trimming.
        - `{vx}`: The x (float) vtkImageData coordinates of the viewport center point after trimming.
        - `{vy}`: The y (float) vtkImageData coordinates of the viewport center point after trimming.
        - `{vz}`: The z (float) vtkImageData coordinates of the viewport center point after trimming.
    qform : vtkMatrix4x4, optional
        Homogeneous vtkImageData coordinates to world transformation matrix.
        For example, pass the vtkNIFTIImageReader.GetQFormMatrix().
    level_window : (float, float), optional
        2-tuple/-list of level and window color transfer function parameters.
        When not specified, the auto_level_window function with default percentiles
        is used to compute an intensity range that is robust to outliers.
    polydata : vtkPolyData, list, optional
        List of vtkPolyData objects to be cut by each orthogonal
        image slice plane and the contours rendered over the image.
        When a `qform` matrix is given, the points are transformed
        to image coordinates using the inverse of the `qform` matrix.
    colors : list, optional
        List of colors (3-tuples of float RGB values in [0, 1]) to use for each `polydata` contour.
    trim : bool, optional
        Whether to trim cropping region such that screenshot is contained within image bounds.
    overwrite : bool, optional
        Whether to overwrite existing output files. When this option is False,
        a screenshot is only taken when the output file does not exist.

    Returns
    -------

    paths: list
        List of 6-tuples with parameters used to render each screenshot and the
        absolute file path of the written PNG image files. The number of screenshots
        is divisable by three, where the first `len(paths)/3` screenshots are taken from
        axial image slices, the next `len(paths)/3` screenshots are taken from coronal
        slices, and the last `len(paths)/3` screenshots are taken from sagittal slices.

        For each screenshot, the returned 6-tuple contains the following values:
        - path:   Absolute path of PNG file.
        - zdir:   Volume dimension corresponding to viewing direction.
        - center: 3-tuple of screenshot center voxel indices before trimming.
        - origin: 3-tuple of viewport center voxel indices after trimming.
        - size:   2-tuple of extracted slice width and height.
        - isnew:  Whether image file was newly written or existed already.
                  This value is always True when `overwrite=True`.

    """
    if isinstance(size, int):
        size = (size, size)
    if not level_window:
        level_window = auto_level_window(image)
    if not suffix or len(suffix) != 3:
        raise Exception("suffix must be a tuple/list of three elements")
    suffix = list(suffix)
    suffix.reverse()
    origin = image.GetOrigin()
    spacing = image.GetSpacing()
    dims = image.GetDimensions()
    if os.path.splitext(path_format)[1] != '.png':
        path_format += '.png'
    if center and (isinstance(center[0], int) and
                   isinstance(center[1], int) and
                   isinstance(center[2], int)):
        index = center
    else:
        if not center:
            center = image.GetCenter()
        index = nearest_voxel(point_to_index(center, origin, spacing))
    if not offsets:
        offsets = [0]
    elif isinstance(offsets, (int, float)):
        offsets = [offsets]
    if qform:
        linear_transform = vtk.vtkMatrixToLinearTransform()
        linear_transform.SetInput(invert_matrix(qform))
        linear_transform.Update()
    else:
        linear_transform = None
    args = dict(
        polydata=polydata,
        transform=linear_transform,
        level_window=level_window,
        colors=colors,
        line_width=line_width
    )
    n = 0
    screenshots = []
    for zdir in (2, 1, 0):
        xdir, ydir = slice_axes(zdir)
        if isinstance(length, int):
            width = length
            height = length
        else:
            width = iround(length / spacing[xdir])
            height = iround(length / spacing[ydir])
        for offset in offsets:
            ijk = list(index)
            if isinstance(offset, int):
                ijk[zdir] += offset
            else:
                ijk[zdir] += iround(offset / spacing[zdir])
            if ijk[zdir] < 0 or ijk[zdir] >= dims[zdir]:
                continue
            vox = list(ijk)
            if trim:
                bounds = cropping_region(vox, (xdir, ydir, zdir), width, height)
                if bounds[0] < 0:
                    bounds[0] = 0
                if bounds[1] >= dims[0]:
                    bounds[1] = dims[0] - 1
                if bounds[2] < 0:
                    bounds[2] = 0
                if bounds[3] >= dims[1]:
                    bounds[3] = dims[1] - 1
                if bounds[4] < 0:
                    bounds[4] = 0
                if bounds[5] >= dims[2]:
                    bounds[5] = dims[2] - 1
                dim = 2 * xdir
                width = (bounds[dim + 1] - bounds[dim])
                vox[xdir] = bounds[dim] + width / 2
                dim = 2 * ydir
                height = (bounds[dim + 1] - bounds[dim])
                vox[ydir] = bounds[dim] + height / 2
            renderer = slice_view(image, zdir=zdir, index=vox, width=width, height=height, **args)
            window = vtk.vtkRenderWindow()
            window.SetSize(size)
            window.AddRenderer(renderer)
            n += 1
            xyz = index_to_point(ijk, origin, spacing)
            pos = index_to_point(vox, origin, spacing)
            path = os.path.abspath(path_format.format(
                n=n, prefix=prefix, suffix=suffix[zdir],
                vi=vox[0], vj=vox[1], vk=vox[2],
                vx=pos[0], vy=pos[1], vz=pos[2],
                i=ijk[0], j=ijk[1], k=ijk[2],
                x=xyz[0], y=xyz[1], z=xyz[2],
            ))
            if overwrite or not os.path.isfile(path):
                directory = os.path.dirname(path)
                if not os.path.isdir(directory):
                    os.makedirs(directory)
                take_screenshot(window, path=path)
                isnew = True
            else:
                isnew = False
            screenshots.append((path, zdir, ijk, vox, (width, height), isnew))
            renderer = None
            window = None
    return screenshots
