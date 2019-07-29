/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/NumericsConfig.h"
#include "mirtk/DeformableConfig.h"
#include "mirtk/TransformationConfig.h"

#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/SurfaceBoundary.h"
#include "mirtk/Transformation.h"
#include "mirtk/BSplineFreeFormTransformation3D.h"
#include "mirtk/BSplineFreeFormTransformationSV.h"
#include "mirtk/NearestNeighborInterpolateImageFunction.h"
#include "mirtk/LinearInterpolateImageFunction.h"

// Deformable surface model / parameterization
#include "mirtk/DeformableSurfaceModel.h"
#include "mirtk/DeformableSurfaceLogger.h"
#include "mirtk/DeformableSurfaceDebugger.h"

// Optimization method
#include "mirtk/LocalOptimizer.h"
#include "mirtk/EulerMethod.h"
#include "mirtk/EulerMethodWithMomentum.h"
#include "mirtk/GradientDescent.h"
#include "mirtk/InexactLineSearch.h"
#include "mirtk/BrentLineSearch.h"

// Stopping criteria
#include "mirtk/EnergyThreshold.h"
#include "mirtk/MinActiveStoppingCriterion.h"
#include "mirtk/InflationStoppingCriterion.h"

// External forces
#include "mirtk/BalloonForce.h"
#include "mirtk/ImageEdgeForce.h"
#include "mirtk/ImageEdgeDistance.h"
#include "mirtk/ImplicitSurfaceDistance.h"

// Internal forces
#include "mirtk/SpringForce.h"
#include "mirtk/InflationForce.h"
#include "mirtk/CurvatureConstraint.h"
#include "mirtk/GaussCurvatureConstraint.h"
#include "mirtk/MeanCurvatureConstraint.h"
#include "mirtk/MaximumCurvatureConstraint.h"
#include "mirtk/QuadraticCurvatureConstraint.h"
#include "mirtk/MetricDistortion.h"
#include "mirtk/StretchingForce.h"
#include "mirtk/RepulsiveForce.h"
#include "mirtk/NonSelfIntersectionConstraint.h"
#include "mirtk/NormalForce.h"

// Transformation constraints
#include "mirtk/SmoothnessConstraint.h"

// VTK
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkCellTreeLocator.h"
#include "vtkSortDataArray.h"
#include "vtkPolyDataNormals.h"


using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  DeformableSurfaceModel model; // with default parameters
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Iteratively minimizes a deformable surface model energy functional. The gradient of" << endl;
  cout << "  the energy terms are the internal and external forces of the deformable surface model." << endl;
  cout << endl;
  cout << "Input options:" << endl;
  cout << "  -initial <file>" << endl;
  cout << "      Point set used to initialize the deformed output mesh. Usually the output of a" << endl;
  cout << "      previous optimization with possibly saved node status (see :option:`-save-status`)." << endl;
  cout << "      (default: input)" << endl;
  cout << "  -dof <type> [<dx> [<dy> <dz>]]" << endl;
  cout << "      Optimize spatial transformation of named <type> to deform the mesh points." << endl;
  cout << "      The optional <dx>, <dy>, and <dz> arguments specify the control point spacing" << endl;
  cout << "      of free-form deformation (FFD) transformations. Common transformation types are:" << endl;
  cout << "      - ``FFD``:   Cubic B-spline FFD." << endl;
  cout << "      - ``SVFFD``: Stationary velocity (SV) cubic B-spline FFD." << endl;
  cout << "  -image <file>" << endl;
  cout << "      Intensity image on which external forces are based. (default: none)" << endl;
  cout << "  -distance-image, -dmap <file>" << endl;
  cout << "      Euclidean distance image on which implicit surface forces are based. (default: none)" << endl;
  cout << "  -distance-offset, -dmap-offset <value>" << endl;
  cout << "      Implicit surface isovalue of :option:`-distance-image`. (default: 0)" << endl;
  cout << "  -mask <file>" << endl;
  cout << "      Mask defining region in which external forces are non-zero. (default: none)" << endl;
  cout << "  -padding <value>" << endl;
  cout << "      Padding/Background value of input :option:`-image`. (default: none)" << endl;
  cout << "  -inflate-brain" << endl;
  cout << "      Set default parameters of cortical surface inflation process equivalent" << endl;
  cout << "      to FreeSurfer's mris_inflate command." << endl;
  cout << endl;
  cout << "Optimization options:" << endl;
  cout << "  -optimizer <name>" << endl;
  cout << "      Optimization method used to minimize the energy of the deformable surface model:" << endl;
  cout << "      - ``EulerMethod``:              Forward Euler integration (default)" << endl;
  cout << "      - ``EulerMethodWithDamping``:   Forward Euler integration with momentum." << endl;
  cout << "      - ``EulerMethodWithMomentum``:  Forward Euler integration with momentum." << endl;
  cout << "      - ``GradientDescent``:          Gradient descent optimizer." << endl;
  cout << "      - ``ConjugateGradientDescent``: Conjugate gradient descent." << endl;
  cout << "  -line-search, -linesearch <name>" << endl;
  cout << "      Line search method used by gradient descent optimizers:" << endl;
  cout << "      - ``Adaptive``: Line search with adaptive step length. (default)" << endl;
  cout << "      - ``Brent``: Brent's line search method." << endl;
  cout << "  -damping <value>" << endl;
  cout << "      Damping ratio used by Euler method with momentum modelling the effect" << endl;
  cout << "      of dissipation of kinetic energy." << endl;
  cout << "  -momentum <value>" << endl;
  cout << "      Momentum of Euler method with momentum, i.e., :math:`1 - damping` (see :option:`-damping`)" << endl;
  cout << "  -mass <value>" << endl;
  cout << "      Node mass used by Euler methods with momentum. (default: 1)" << endl;
  cout << "  -levels <max> | <min> <max>" << endl;
  cout << "      Perform optimization on starting at level <max> until level <min> (> 0)." << endl;
  cout << "      When only the <max> level argument is given, the <min> level is set to 1." << endl;
  cout << "      On each level, the node forces are averaged :math:`2^{level-1}` times which" << endl;
  cout << "      is similar to computing the forces on a coarser mesh. See :option:`-force-averaging`. (default: 0 0)" << endl;
  cout << "  -force-averaging <n>..." << endl;
  cout << "      Number of force averaging steps. (default: 0)" << endl;
  cout << "      Cannot be combined with :option:`-magnitude-averaging`." << endl;
  cout << "  -magnitude-averaging <n>..." << endl;
  cout << "      Number of force magnitude averaging steps. (default: 0)" << endl;
  cout << "      Cannot be combined with :option:`-force-averaging`." << endl;
  cout << "  -distance-averaging <n>..." << endl;
  cout << "      Number of :option:`-distance` force averaging steps. (default: 0)" << endl;
  cout << "  -steps, -max-steps, -iterations, -max-iterations <n>..." << endl;
  cout << "      Maximum number of iterations. (default: 100)" << endl;
  cout << "  -step, -dt <value>..." << endl;
  cout << "      Length of integration/gradient steps. (default: 1)" << endl;
  cout << "  -max-dx, -maxdx, -dx <value>..." << endl;
  cout << "      Maximum displacement of a node at each iteration. By default, the node displacements" << endl;
  cout << "      are normalized by the maximum node displacement. When this option is used, the node" << endl;
  cout << "      displacements are clamped to the specified maximum length instead. (default: :option:`-step`)" << endl;
  cout << "  -max-displacement <value>" << endl;
  cout << "      Maximum distance from input surface. (default: +inf)" << endl;
  cout << "  -remesh <n>" << endl;
  cout << "      Remesh surface mesh every n-th iteration. (default: " << model.RemeshInterval() << ")" << endl;
  cout << "  -remesh-adaptively" << endl;
  cout << "      Remesh surface mesh using an adaptive edge length interval based on local curvature" << endl;
  cout << "      of the deformed surface mesh or input implicit surface (:option:`-distance-image`)." << endl;
  cout << "  -[no]triangle-inversion" << endl;
  cout << "      Whether to allow inversion of pair of triangles during surface remeshing. (default: on)" << endl;
  cout << "  -min-edgelength <value>..." << endl;
  cout << "      Minimum edge length used for local adaptive remeshing. (default: " << model.MinEdgeLength() << ")" << endl;
  cout << "  -max-edgelength <value>..." << endl;
  cout << "      Maximum edge length used for local adaptive remeshing. (default: " << model.MaxEdgeLength() << ")" << endl;
  cout << "  -min-angle <degrees>..." << endl;
  cout << "      Minimum angle between edge node normals for an edge be excluded from collapsing during" << endl;
  cout << "      iterative :option:`-remesh` operations. (default: " << model.MinFeatureAngle() << ")" << endl;
  cout << "  -max-angle <degrees>..." << endl;
  cout << "      Maximum angle between edge node normals for an edge be excluded from splitting during" << endl;
  cout << "      iterative :option:`-remesh` operations. (default: " << model.MaxFeatureAngle() << ")" << endl;
  cout << "  -lowpass <n>" << endl;
  cout << "      Low-pass filter surface mesh every n-th iteration. (default: " << model.LowPassInterval() << ")" << endl;
  cout << "  -lowpass-iterations <n>" << endl;
  cout << "      Number of :option:`-lowpass` filter iterations. (default: " << model.LowPassIterations() << ")" << endl;
  cout << "  -lowpass-band <band>" << endl;
  cout << "      Low-pass filtering band argument, usually in the range [0, 2]. (default: " << model.LowPassBand() << ")" << endl;
  cout << "  -nointersection" << endl;
  cout << "      Hard non-self-intersection constraint for surface meshes. (default: off)" << endl;
  cout << "  -mind, -min-distance <value>" << endl;
  cout << "      Minimum distance to other triangles in front of a given triangle." << endl;
  cout << "  -minw, -min-width <value>" << endl;
  cout << "      Minimum distance to other triangles in the back of a given triangle." << endl;
  cout << "  -max-collision-angle <degrees>" << endl;
  cout << "      Maximum angle between vector connecting centers of nearby triangles and the face normal" << endl;
  cout << "      of the reference triangle for a collision to be detected. When the triangles are within" << endl;
  cout << "      the same flat neighborhood of the surface mesh, this angle will be close to 90 degrees." << endl;
  cout << "      This parameter reduces false collision detection between neighboring triangles. (default: " << model.MaxCollisionAngle() << ")" << endl;
  cout << "  -fast-collision-test\n";
  cout << "      Use fast approximate triangle-triangle collision test based on distance of their centers only. (default: off)" << endl;
  cout << "  -reset-status" << endl;
  cout << "      Set status of all mesh nodes to active again after each level (see :option:`-levels`). (default: off)" << endl;
  cout << endl;
  cout << "Deformable model options:" << endl;
  cout << "  -neighborhood <n>" << endl;
  cout << "      Size of node neighborhoods used by internal force terms that consider more" << endl;
  cout << "      than only the adjacent nodes, but also up to n-connected nodes. (default: " << model.NeighborhoodRadius() << ")" << endl;
  cout << "  -distance <w>" << endl;
  cout << "      Weight of implicit surface distance. (default: 0)" << endl;
  cout << "  -distance-measure <name>" << endl;
  cout << "      Implicit surface distance measure used by :option:`-distance`:" << endl;
  cout << "      - ``minimum``: Minimum surface distance (see :option:`-distance-image`, default)" << endl;
  cout << "      - ``normal``:  Estimate distance by casting rays along normal direction." << endl;
  cout << "  -balloon-inflation, -balloon <w>" << endl;
  cout << "      Weight of inflation force based on local intensity statistics. (default: 0)" << endl;
  cout << "  -balloon-deflation <w>" << endl;
  cout << "      Weight of deflation force based on local intensity statistics. (default: 0)" << endl;
  cout << "  -balloon-min <intensity>" << endl;
  cout << "      Global lower intensity threshold for :option:`-balloon-inflation` or :option:`-balloon-deflation`. (default: -inf)" << endl;
  cout << "  -balloon-max <intensity>" << endl;
  cout << "      Global upper intensity threshold for :option:`-balloon-inflation` or :option:`-balloon-deflation`. (default: +inf)" << endl;
  cout << "  -balloon-range <min> <max>" << endl;
  cout << "      Global intensity thresholds for :option:`-balloon-inflation` or :option:`-balloon-deflation`. (default: [-inf +inf])" << endl;
  cout << "  -balloon-radius <r>" << endl;
  cout << "      Radius for local intensity statistics of :option:`-balloon-inflation` or :option:`-balloon-deflation`. (default: 7 times voxel size)" << endl;
  cout << "  -balloon-sigma <sigma>" << endl;
  cout << "      Local intensity standard deviation scaling factor of :option:`-balloon-inflation` or :option:`-balloon-deflation`. (default: 5)" << endl;
  cout << "  -balloon-mask <file>" << endl;
  cout << "      Image mask used for local intensity statistics for :option:`-balloon-inflation` or :option:`-balloon-deflation`." << endl;
  cout << "      (default: interior of deformed surface)" << endl;
  cout << "  -edges <w>" << endl;
  cout << "      Weight of image edge force. (default: 0)" << endl;
  cout << "  -edge-distance <w>" << endl;
  cout << "      Weight of closest image edge distance force. (default: 0)" << endl;
  cout << "  -inflation <w>" << endl;
  cout << "      Weight of surface inflation force used for cortical surface inflation. (default: 0)" << endl;
  cout << "  -bending-energy <w>" << endl;
  cout << "      Weight of bending energy of :option:`-dof` transformation. (default: 0)" << endl;
  cout << "  -spring <w>" << endl;
  cout << "      Weight of internal spring force. (default: 0)" << endl;
  cout << "  -normal-spring, -nspring <w>" << endl;
  cout << "      Weight of internal spring force in normal direction. (default: 0)" << endl;
  cout << "  -tangential-spring, -tspring <w>" << endl;
  cout << "      Weight of internal spring force in tangent plane. (default: 0)" << endl;
  cout << "  -normalized-spring <w>" << endl;
  cout << "      Weight of internal spring force normalized w.r.t. force in normal direction. (default: 0)" << endl;
  cout << "  -curvature <w>" << endl;
  cout << "      Weight of surface curvature. (default: 0)" << endl;
  cout << "  -quadratic-curvature, -qcurvature <w>" << endl;
  cout << "      Weight of surface curvature estimated by quadratic fit of node neighbor to tangent plane distance. (default: 0)" << endl;
  cout << "  -distant-quadratic-curvature, -distant-qcurvature <w>" << endl;
  cout << "      Weight of :option:`-quadratic-curvature` proportional to the :option:`-distance` magnitude. (default: 0)" << endl;
  cout << "  -mean-curvature, -mcurvature <w>" << endl;
  cout << "      Weight of mean curvature constraint. (default: 0)" << endl;
  cout << "  -distant-mean-curvature, -distant-mcurvature <w>" << endl;
  cout << "      Weight of :option:`-mean-curvature` proportional to the :option:`-distance` magnitude. (default: 0)" << endl;
  cout << "  -gauss-curvature, -gcurvature <w>" << endl;
  cout << "      Weight of Gauss curvature constraint. (default: 0)" << endl;
  cout << "  -distant-gauss-curvature, -distant-gcurvature <w>" << endl;
  cout << "      Weight of :option:`-gauss-curvature` proportional to the :option:`-distance` magnitude. (default: 0)" << endl;
  cout << "  -distortion <w>" << endl;
  cout << "      Weight of metric distortion." << endl;
  cout << "  -stretching <w>" << endl;
  cout << "      Weight of spring force based on difference of neighbor distance compared to" << endl;
  cout << "      initial distance. (default: 0)" << endl;
  cout << "  -repulsion <w>" << endl;
  cout << "      Weight of node repulsion force. (default: 0)" << endl;
  cout << "  -repulsion-radius <r>" << endl;
  cout << "      Radius of node repulsion force. (default: average edge length)" << endl;
  cout << "  -repulsion-distance <r>" << endl;
  cout << "      Frontface radius of node repulsion force. (default: average edge length)" << endl;
  cout << "  -repulsion-width <r>" << endl;
  cout << "      Backface radius of node repulsion force. (default: average edge length)" << endl;
  cout << "  -collision <w>" << endl;
  cout << "      Weight of triangle repulsion force." << endl;
  cout << "  -normal <w>" << endl;
  cout << "      Constant force along outwards normal direction (positive weight)" << endl;
  cout << "      or inwards normal direction (negative weight)." << endl;
  cout << endl;
  cout << "Stopping criterion options:" << endl;
  cout << "  -extrinsic-energy" << endl;
  cout << "      Consider only sum of external energy terms as total energy value of deformable model functional." << endl;
  cout << "      Internal forces still contribute to the gradient of the functional, but are excluded from the" << endl;
  cout << "      energy function value (see :option:`-epsilon` and :option:`-min-energy`). (default: off)" << endl;
  cout << "  -epsilon <value>" << endl;
  cout << "      Minimum change of deformable surface energy convergence criterion." << endl;
  cout << "  -delta <value>..." << endl;
  cout << "      Minimum maximum node displacement or :option:`-dof` parameter value. (default: 1e-6)" << endl;
  cout << "  -min-energy <value>" << endl;
  cout << "      Target deformable surface energy value. (default: 0)" << endl;
  cout << "  -min-active <ratio>..." << endl;
  cout << "      Minimum ratio of active nodes in [0, 1]. (default: 0)" << endl;
  cout << "  -inflation-error <threshold>" << endl;
  cout << "      Threshold of surface inflation RMS measure. (default: off)" << endl;
  cout << endl;
  cout << "Output options:" << endl;
  cout << "  -track [<name>]" << endl;
  cout << "      Record sum of node displacements along normal direction. The integrated" << endl;
  cout << "      displacements are stored in the point data array named \"NormalDisplacement\"" << endl;
  cout << "      by default or with the specified <name>. (default: off)" << endl;
  cout << "  -track-zero-mean [<name>]" << endl;
  cout << "      Same as :option:`-track`, but subtract mean from tracked values." << endl;
  cout << "      This option is implicit when :option:`-inflate-brain` is given to inflate a cortical" << endl;
  cout << "      brain surface. The default output array name is \"SulcalDepth\". Otherwise, the default" << endl;
  cout << "      point data array name is \"NormalDisplacementZeroMean\". (default: off)" << endl;
  cout << "  -track-zero-median [<name>]" << endl;
  cout << "      Same as :option:`-track`, but subtract median from tracked values." << endl;
  cout << "      The default point data array name is \"NormalDisplacementZeroMedian\". (default: off)" << endl;
  cout << "  -track-unit-variance [<name>]" << endl;
  cout << "      Same as :option:`-track`, but divide tracked values by their standard deviation." << endl;
  cout << "      The default point data array name is \"NormalDisplacementUnitVariance\". (default: off)" << endl;
  cout << "  -track-zvalues [<name>]" << endl;
  cout << "      Same as :option:`-track`, but subtract mean from tracked values and divide by standard deviation." << endl;
  cout << "      The resulting values are the Z-score normalized standard scores of the tracked displacements." << endl;
  cout << "      The default point data array name is \"NormalDisplacementZValues\". (default: off)" << endl;
  cout << "  -track-zero-median-zvalues [<name>]" << endl;
  cout << "      Same as :option:`-track`, but subtract median from tracked values and divide by standard deviation." << endl;
  cout << "      It can be used with :option:`-inflate-brain` to obtain a normalized curvature measure." << endl;
  cout << "      The default point data array name is \"NormalDisplacementZeroMedianZValues\". (default: off)" << endl;
  cout << "  -track-without-momentum" << endl;
  cout << "      When tracking the total displacement of a node in normal direction using the EulerMethodWithMomentum" << endl;
  cout << "      :option:`-optimizer` as used in particular by :option:`-inflate-brain`, exclude the momentum from the." << endl;
  cout << "      tracked displacement. This is idential to FreeSurfer's mrisTrackTotalDisplacement used for the curvature" << endl;
  cout << "      output of mris_inflate. The correct curvature value is, however, obtained by including the momentum" << endl;
  cout << "      component as it integrates the actual amount by which each node is displaced during the Euler steps." << endl;
  cout << "      Not using this option corresponds to the mrisTrackTotalDisplacementNew function. (default: off)" << endl;
  cout << "  -notrack" << endl;
  cout << "      Do not track node displacements along normal direction." << endl;
  cout << "  -center-output" << endl;
  cout << "      Center output mesh such that center is at origin. (default: off)" << endl;
  cout << "  -match-area" << endl;
  cout << "      Scale output mesh by ratio of input and output surface area. (default: off)" << endl;
  cout << "  -match-sampling" << endl;
  cout << "      Resample output mesh at corresponding positions of input mesh." << endl;
  cout << "      This option is only useful in conjunction with :option:`-remesh`. (default: off)" << endl;
  cout << "  -save-status" << endl;
  cout << "      Save node status (active/passive) to output file. (default: off)" << endl;
  cout << "  -ascii | -nobinary" << endl;
  cout << "      Write legacy VTK in ASCII format. (default: off)" << endl;
  cout << "  -binary | -noascii" << endl;
  cout << "      Write legacy VTK in binary format. (default: on)" << endl;
  cout << "  -[no]compress" << endl;
  cout << "      Write XML VTK file with or without compression. (default: on)" << endl;
  cout << "  -debug-prefix <prefix>" << endl;
  cout << "      File name prefix for :option:`-debug` output. (default: deform_mesh\\_)" << endl;
  cout << "  -debug-interval <n>" << endl;
  cout << "      Write :option:`-debug` output every n-th iteration. (default: 10)" << endl;
  cout << "  -[no]level-prefix" << endl;
  cout << "      Write :option:`-debug` output without level prefix in file names. (default: on)" << endl;
  cout << endl;
  cout << "Advanced options:" << endl;
  cout << "  -par <name> <value>" << endl;
  cout << "      Advanced option to set surface model or optimizer parameter." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Resample image
template <class VoxelType>
void ResampleImage(GenericImage<VoxelType> &image, const ImageAttributes &attr)
{
  const GenericImage<VoxelType> input(image);
  GenericLinearInterpolateImageFunction<GenericImage<VoxelType> > func;
  func.Input(&input);
  func.Initialize();
  image.Initialize(attr, 1);
  double x, y, z;
  for (int k = 0; k < image.Z(); ++k)
  for (int j = 0; j < image.Y(); ++j)
  for (int i = 0; i < image.X(); ++i) {
    x = i, y = j, z = k;
    image.ImageToWorld(x, y, z);
    func.WorldToImage(x, y, z);
    image(i, j, k) = static_cast<VoxelType>(func.Evaluate(x, y, z));
  }
}

// -----------------------------------------------------------------------------
/// Resample mask
void ResampleMask(BinaryImage &mask, const ImageAttributes &attr)
{
  RealImage resampled(mask);
  ResampleImage(resampled, attr);
  mask.Initialize(resampled.Attributes());
  const int nvox = resampled.NumberOfVoxels();
  for (int vox = 0; vox < nvox; ++vox) {
    mask(vox) = BinaryPixel(resampled(vox) >= .5 ? 1 : 0);
  }
}

// -----------------------------------------------------------------------------
/// Get parameter for current level
///
/// Some parameters can be set per level. The number of levels is defined by
/// the maximum number of values for each parameter. When the user did not
/// specify a parameter for each level, the first value is used for the initial
/// levels for which a parameter value is missing.
template <class T>
T ParameterValue(int level, int nlevels, const Array<T> &values, T default_value)
{
  mirtkAssert(nlevels > 0,                   "at least one level");
  mirtkAssert(0 <= level && level < nlevels, "valid level index");
  const int nvalues = static_cast<int>(values.size());
  if (nvalues == 0) return default_value;
  const int idx = level - (nlevels - nvalues);
  return values[idx < 0 ? 0 : idx];
}

// -----------------------------------------------------------------------------
/// Initialize array of total energy gradient averaging steps
Array<int> GradientAveraging(int min_level, int max_level)
{
  Array<int> navgs;
  navgs.reserve(max_level - min_level + 1);
  for (int level = max_level; level >= min_level; --level) {
    navgs.push_back(level > 1 ? static_cast<int>(pow(2, level - 2)) : 0);
  }
  return navgs;
}

// -----------------------------------------------------------------------------
/// Subract mean of tracked normal displacements
void DemeanValues(vtkDataArray *values, bool use_median = false)
{
  double mu;

  // Compute mean (or median)
  if (use_median) {
    vtkSmartPointer<vtkDataArray> sorted;
    sorted.TakeReference(values->NewInstance());
    sorted->DeepCopy(values);
    vtkSortDataArray::Sort(sorted);
    mu = sorted->GetComponent(sorted->GetNumberOfTuples() / 2, 0);
  } else {
    mu = .0;
    for (vtkIdType id = 0; id < values->GetNumberOfTuples(); ++id) {
      mu += values->GetComponent(id, 0);
    }
    mu /= values->GetNumberOfTuples();
  }

  // Subtract mean from curvature measures
  for (vtkIdType id = 0; id < values->GetNumberOfTuples(); ++id) {
    values->SetComponent(id, 0, values->GetComponent(id, 0) - mu);
  }
}

// -----------------------------------------------------------------------------
/// Normalize variance of tracked normal displacements
void NormalizeVariance(vtkDataArray *values)
{
  // Compute mean
  double mu = .0;
  for (vtkIdType id = 0; id < values->GetNumberOfTuples(); ++id) {
    mu += values->GetComponent(id, 0);
  }
  mu /= values->GetNumberOfTuples();

  // Compute variance
  double var = .0;
  for (vtkIdType id = 0; id < values->GetNumberOfTuples(); ++id) {
    var += pow(values->GetComponent(id, 0) - mu, 2);
  }
  var /= values->GetNumberOfTuples();

  // Normalize variance
  double sigma = (var == .0 ? 1.0 : sqrt(var));
  for (vtkIdType id = 0; id < values->GetNumberOfTuples(); ++id) {
    values->SetComponent(id, 0, mu + (values->GetComponent(id, 0) - mu) / sigma);
  }
}

// -----------------------------------------------------------------------------
/// Normalize tracked normal displacements
void NormalizeValues(vtkDataArray *values, bool use_median = false)
{
  double mu, var;

  // Compute mean (or median)
  if (use_median) {
    vtkSmartPointer<vtkDataArray> sorted;
    sorted.TakeReference(values->NewInstance());
    sorted->DeepCopy(values);
    vtkSortDataArray::Sort(sorted);
    mu = sorted->GetComponent(sorted->GetNumberOfTuples() / 2, 0);
  } else {
    mu = .0;
    for (vtkIdType id = 0; id < values->GetNumberOfTuples(); ++id) {
      mu += values->GetComponent(id, 0);
    }
    mu /= values->GetNumberOfTuples();
  }

  // Compute variance
  var = .0;
  for (vtkIdType id = 0; id < values->GetNumberOfTuples(); ++id) {
    var += pow(values->GetComponent(id, 0) - mu, 2);
  }
  var /= values->GetNumberOfTuples();

  // Z-score normalize curvature measures
  double sigma = (var == .0 ? 1.0 : sqrt(var));
  for (vtkIdType id = 0; id < values->GetNumberOfTuples(); ++id) {
    values->SetComponent(id, 0, (values->GetComponent(id, 0) - mu) / sigma);
  }
}

// -----------------------------------------------------------------------------
/// Resample surface mesh at corresponding positions of the initial points
///
/// \todo Resample also output point data.
vtkSmartPointer<vtkPointSet>
ResampleAtInitialPoints(vtkSmartPointer<vtkPointSet> input, vtkSmartPointer<vtkPointSet> output)
{
  double    p[3], x[3], pcoords[3], dist2, *weights;
  vtkIdType cellId;
  int       subId;

  vtkDataArray *initial_position = output->GetPointData()->GetArray("InitialPoints");
  if (!initial_position) {
    Warning("Cannot resample surface mesh at points corresponding to points of input mesh:"
            " deformed mesh has not point data array named \"InitialPoints\".");
    return output;
  }

  vtkSmartPointer<vtkPoints> initial_points = vtkSmartPointer<vtkPoints>::New();
  initial_points->SetNumberOfPoints(output->GetNumberOfPoints());
  for (vtkIdType ptId = 0; ptId < output->GetNumberOfPoints(); ++ptId) {
    initial_points->SetPoint(ptId, initial_position->GetTuple(ptId));
  }

  vtkSmartPointer<vtkPointSet> initial;
  initial.TakeReference(output->NewInstance());
  initial->ShallowCopy(output);
  initial->SetPoints(initial_points);

  vtkSmartPointer<vtkAbstractCellLocator> locator;
  locator = vtkSmartPointer<vtkCellTreeLocator>::New();
  locator->SetDataSet(initial);
  locator->BuildLocator();

  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  weights = new double[output->GetMaxCellSize()];

  vtkSmartPointer<vtkPoints> resampled_points = vtkSmartPointer<vtkPoints>::New();
  resampled_points->SetNumberOfPoints(input->GetNumberOfPoints());

  for (vtkIdType ptId = 0; ptId < input->GetNumberOfPoints(); ++ptId) {
    input->GetPoint(ptId, p);
    locator->FindClosestPoint(p, x, cell, cellId, subId, dist2);
    cell->EvaluatePosition(x, NULL, subId, pcoords, dist2, weights);
    p[0] = p[1] = p[2];
    for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); ++i) {
      output->GetPoint(cell->GetPointId(i), x);
      p[0] += weights[i] * x[0];
      p[1] += weights[i] * x[1];
      p[2] += weights[i] * x[2];
    }
    resampled_points->SetPoint(ptId, p);
  }

  delete[] weights;

  vtkSmartPointer<vtkPointSet> resampled;
  resampled.TakeReference(input->NewInstance());
  resampled->ShallowCopy(input);
  resampled->SetPoints(resampled_points);
  return resampled;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
/// Helper macro to parse multiple parameter arguments
#define PARSE_ARGUMENTS(T, name) \
  do { \
    T value; \
    PARSE_ARGUMENT(value); \
    (name).push_back(value); \
  } while (HAS_ARGUMENT); \
  nlevels = max(nlevels, static_cast<int>((name).size()))

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  FileOption output_fopt = FO_Default;

  verbose = 1; // default verbosity level
  EXPECTS_POSARGS(2);

  // Initialize libraries / object factories
  InitializeIOLibrary();
  InitializeNumericsLibrary();
  InitializeDeformableLibrary();
  InitializeTransformationLibrary();

  // Deformable surface model and default optimizer
  UniquePtr<Transformation> dof;
  DeformableSurfaceModel    model;
  DeformableSurfaceLogger   logger;
  DeformableSurfaceDebugger debugger(&model);
  UniquePtr<LocalOptimizer> optimizer(new EulerMethod(&model));
  ParameterList             params;

  // Read input point set
  vtkSmartPointer<vtkPointSet> input = ReadPointSet(POSARG(1), output_fopt);
  vtkPointData * const inputPD = input->GetPointData();
  vtkCellData  * const inputCD = input->GetCellData();
  ImageAttributes domain = PointSetDomain(input);
  model.Input(input);

  // External forces (inactive if weight == 0)
  BalloonForce            balloon ("Balloon force", .0);
  ImageEdgeForce          edges   ("Edge force",    .0);
  ImageEdgeDistance       dedges  ("Edge distance", .0);
  ImplicitSurfaceDistance distance("Distance",      .0);

  // Internal forces (inactive if weight == 0)
  SpringForce                   spring    ("Bending",        .0);
  InflationForce                normspring("Bending",        .0);
  InflationForce                inflation ("Inflation",      .0);
  CurvatureConstraint           curvature ("Curvature",      .0);
  GaussCurvatureConstraint      gcurvature("Gauss curv.",    .0);
  MeanCurvatureConstraint       mcurvature("Mean curv.",     .0);
  MaximumCurvatureConstraint    pcurvature("Max. curv.",     .0);
  QuadraticCurvatureConstraint  qcurvature("Quad. curv.",    .0);
  MetricDistortion              distortion("Distortion",     .0);
  StretchingForce               stretching("Stretching",     .0);
  RepulsiveForce                repulsion ("Repulsion",      .0);
  NonSelfIntersectionConstraint collision ("Collision",      .0);
  SmoothnessConstraint          dofbending("Bending energy", .0);
  NormalForce                   nforce    ("Normal",         .0);

  // Default force settings, also set after option parsing if unchanged
  gcurvature.WeightInside (0.);
  gcurvature.WeightOutside(0.);

  mcurvature.WeightInside (0.);
  mcurvature.WeightOutside(0.);

  pcurvature.WeightInside (0.);
  pcurvature.WeightOutside(0.);

  qcurvature.WeightInside (0.);
  qcurvature.WeightOutside(0.);

  spring.InwardNormalWeight (0.);
  spring.OutwardNormalWeight(0.);
  spring.TangentialWeight   (0.);

  repulsion.FrontfaceRadius(0.);
  repulsion.BackfaceRadius (0.);

  // Stopping criteria (disabled by default)
  EnergyThreshold            min_energy(&model);
  MinActiveStoppingCriterion min_active(&model);
  InflationStoppingCriterion inflation_error(&model);

  min_active.Threshold(0.);
  inflation_error.Threshold(NaN);

  // Optional arguments
  const char *image_name        = nullptr;
  const char *balloon_mask_name = nullptr;
  const char *dmap_name         = nullptr;
  const char *dmag_name         = nullptr;
  const char *mask_name         = nullptr;
  const char *track_name        = nullptr; // track normal movement of nodes
  bool        track_zero_mean   = false;   // subtract mean from tracked normal movements
  bool        track_unit_var    = false;   // Z-score normalize tracked normal movements
  bool        track_use_median  = false;   // use median instead of mean for normalization
  const char *initial_name      = nullptr;
  const char *debug_prefix      = "deform-mesh_";
  double      padding           = NaN;
  bool        level_prefix      = true;
  bool        reset_status      = false;
  bool        fix_boundary      = false;
  bool        center_output     = false;
  bool        match_area        = false;
  bool        match_sampling    = false;
  bool        signed_gradient   = false; // average gradient vectors with positive dot product
  bool        average_magnitude = false; // average gradient magnitude only
  bool        save_status       = false;
  bool        inflate_brain     = false; // mimick mris_inflate
  int         nlevels           = 1;     // no. of levels

  const char *t1w_image_name       = nullptr;
  const char *t2w_image_name       = nullptr;
  const char *wm_mask_name         = nullptr;
  const char *gm_mask_name         = nullptr;
  const char *cortex_dmap_name     = nullptr;
  const char *vents_dmap_name      = nullptr;
  const char *cerebellum_dmap_name = nullptr;

  Array<int>    navgs;           // no. of total gradient averaging steps
  Array<int>    distance_navgs;  // no. of distance gradient averaging steps
  Array<int>    dedges_navgs;    // no. of edge distance gradient averaging steps
  Array<int>    balloon_navgs;   // no. of balloon force gradient averaging steps
  Array<int>    nsteps;          // maximum no. of integration steps
  Array<double> delta;           // max node change convergence criterion
  Array<double> max_dt;          // maximum integration step length
  Array<double> max_dx;          // maximum node displacement at each integration step
  Array<double> min_edge_length; // minimum average edge length
  Array<double> max_edge_length; // maximum average edge length
  Array<double> min_edge_angle;  // minimum angle between edge end points to allow melting
  Array<double> max_edge_angle;  // maximum angle between edge end points before enforcing subdivision
  Array<double> min_distance;    // minimum front facing distance
  Array<double> min_width;       // minimum back facing distance
  Array<double> dmap_offsets;    // implicit surface distance offset
  Array<double> min_active_thres; // minimum active stopping criterion

  bool   barg;
  int    iarg;
  double farg;

  for (ALL_OPTIONS) {

    // Note: Need to split if-else if blocks to not reach compiler limit
    //       of maximum number of nested blocks.
    bool unknown_option;

    // Input
    unknown_option = false;
    if (OPTION("-image")) {
      image_name = ARGUMENT;
    }
    else if (OPTION("-t1-weighted-image") || OPTION("-t1w-image") || OPTION("-t1-image") || OPTION("-t1w")) {
      t1w_image_name = ARGUMENT;
    }
    else if (OPTION("-t2-weighted-image") || OPTION("-t2w-image") || OPTION("-t2-image") || OPTION("-t2w")) {
      t2w_image_name = ARGUMENT;
    }
    else if (OPTION("-dmap") || OPTION("-distance-map") || OPTION("-distance-image") || OPTION("-implicit-surface")) {
      dmap_name = ARGUMENT;
    }
    else if (OPTION("-dmap-offset") || OPTION("-distance-offset") || OPTION("-implicit-surface-offset")) {
      PARSE_ARGUMENTS(double, dmap_offsets);
    }
    else if (OPTION("-dmag") || OPTION("-distance-magnitude")) {
      dmag_name = ARGUMENT;
    }
    else if (OPTION("-mask")) {
      mask_name = ARGUMENT;
    }
    else if (OPTION("-white-matter-mask") || OPTION("-wm-mask")) {
      wm_mask_name = ARGUMENT;
    }
    else if (OPTION("-grey-matter-mask") || OPTION("-gm-mask")) {
      gm_mask_name = ARGUMENT;
    }
    else if (OPTION("-inner-cortical-dmap") || OPTION("-inner-cortical-distance-map") || OPTION("-inner-cortical-distance-image")) {
      cortex_dmap_name = ARGUMENT;
    }
    else if (OPTION("-ventricles-dmap") || OPTION("-ventricles-distance-map") || OPTION("-ventricles-distance-image")) {
      vents_dmap_name = ARGUMENT;
    }
    else if (OPTION("-cerebellum-dmap") || OPTION("-cerebellum-distance-map") || OPTION("-cerebellum-distance-image")) {
      cerebellum_dmap_name = ARGUMENT;
    }
    else if (OPTION("-initial")) {
      initial_name = ARGUMENT;
    }
    else if (OPTION("-padding")) {
      PARSE_ARGUMENT(padding);
    }
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Presets
    unknown_option = false;
    if (OPTION("-inflate-brain")) { // cf. FreeSurfer's mris_inflate
      inflate_brain = true;
      nlevels = 6;
      navgs = GradientAveraging(1, 6);
      nsteps.resize(1);
      nsteps[0] = 10;
      inflation.Weight(.5);   //  1 / 2 b/c InflationForce   gradient weight incl. factor 2
      distortion.Weight(.05); // .1 / 2 b/c MetricDistortion gradient weight incl. factor 2
      inflation_error.Threshold(.015);
      max_dt.resize(1);
      max_dt[0] = .9;
      model.NeighborhoodRadius(2);
      UniquePtr<EulerMethodWithMomentum> euler(new EulerMethodWithMomentum());
      euler->Momentum(.9);
      euler->NormalizeStepLength(false);
      euler->MaximumDisplacement(1.0);
      optimizer.reset(euler.release());
      center_output    = true;
      match_area       = true;
      signed_gradient  = false;
      track_name       = "SulcalDepth";
      track_zero_mean  = true;
      track_unit_var   = false;
      track_use_median = false;
    }
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Optimization method
    unknown_option = false;
    if (OPTION("-optimizer") || OPTION("-optimiser")) {
      OptimizationMethod m;
      PARSE_ARGUMENT(m);
      optimizer.reset(LocalOptimizer::New(m, &model));
    }
    else if (OPTION("-line-search") || OPTION("-linesearch")) {
      Insert(params, "Line search strategy", ARGUMENT);
    }
    else if (OPTION("-dof")) {
      string arg = ARGUMENT;
      double dx = 1.0, dy = 1.0, dz = 1.0;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(dx);
        if (HAS_ARGUMENT) {
          PARSE_ARGUMENT(dy);
          PARSE_ARGUMENT(dz);
        } else {
          dy = dz = dx;
        }
      }
      string larg = ToLower(arg);
      if (larg == "none") {
        dof.reset(NULL);
      } else if (larg == "ffd") {
        dof.reset(new BSplineFreeFormTransformation3D(domain, dx, dy, dz));
      } else if (larg == "svffd") {
        dof.reset(new BSplineFreeFormTransformationSV(domain, dx, dy, dz));
      } else {
        TransformationType type = Transformation::TypeOfClass(arg.c_str());
        if (type == TRANSFORMATION_UNKNOWN) {
          FatalError("Invalid -dof transformation type argument: " << arg);
        }
        dof.reset(Transformation::New(type));
      }
      model.Transformation(dof.get());
    }
    else if (OPTION("-levels")) {
      int min_level, max_level;
      PARSE_ARGUMENT(min_level);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(max_level);
      } else {
        max_level = min_level;
        min_level = 1;
      }
      if (min_level < 1 || max_level < 1) {
        FatalError("Invalid -levels argument");
      }
      navgs   = GradientAveraging(min_level, max_level);
      nlevels = max(nlevels, static_cast<int>(navgs.size()));
    }
    else if (OPTION("-force-averaging")) {
      PARSE_ARGUMENTS(int, navgs);
      average_magnitude = false;
    }
    else if (OPTION("-magnitude-averaging")) {
      PARSE_ARGUMENTS(int, navgs);
      average_magnitude = true;
    }
    else HANDLE_BOOLEAN_OPTION("signed-averaging", signed_gradient);
    else if (OPTION("-max-steps")      || OPTION("-steps")      ||
             OPTION("-max-iterations") || OPTION("-iterations") ||
             OPTION("-max-iter")       || OPTION("-iter")) {
      PARSE_ARGUMENTS(int, nsteps);
    }
    else if (OPTION("-step") || OPTION("-dt") || OPTION("-h")) {
      PARSE_ARGUMENTS(double, max_dt);
    }
    else if (OPTION("-max-dx") || OPTION("-maxdx") || OPTION("-dx") || OPTION("-maxd") || OPTION("-d")) {
      PARSE_ARGUMENTS(double, max_dx);
    }
    else if (OPTION("-max-displacement") || OPTION("-maxdisplacement")) {
      PARSE_ARGUMENT(model.MaxInputDistance());
    }
    else if (OPTION("-normalize-forces") || OPTION("-normalise-forces")) {
      Insert(params, "Normalize gradient vectors", true);
    }
    else if (OPTION("-damping"))   Insert(params, "Deformable surface damping", ARGUMENT);
    else if (OPTION("-momentum"))  Insert(params, "Deformable surface momentum", ARGUMENT);
    else if (OPTION("-mass"))      Insert(params, "Deformable surface mass", ARGUMENT);
    else if (OPTION("-epsilon"))   Insert(params, "Epsilon", ARGUMENT);
    else if (OPTION("-delta")) {
      PARSE_ARGUMENTS(double, delta);
    }
    else if (OPTION("-min-energy") || OPTION("-minenergy")) {
      PARSE_ARGUMENT(min_energy.Threshold());
    }
    else if (OPTION("-min-active") || OPTION("-minactive")) {
      do {
        string value, units = ValueUnits(ARGUMENT, &value);
        if (!FromString(value, farg)) {
          FatalError("Invalid -min-active value: " << value);
        }
        if (units == "%") {
          farg /= 100.;
        } else if (!units.empty()) {
          FatalError("Invalid -min-active units: " << units);
        }
        min_active_thres.push_back(farg);
      } while (HAS_ARGUMENT);
      nlevels = max(nlevels, static_cast<int>(min_active_thres.size()));
    }
    else if (OPTION("-extrinsic-energy")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(barg);
      else barg = true;
      model.MinimizeExtrinsicEnergy(barg);
    }
    else HANDLE_BOOLEAN_OPTION("reset-status", reset_status);
    else HANDLE_BOOLEAN_OPTION("fix-boundary", fix_boundary);
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // External forces
    unknown_option = false;
    if (OPTION("-distance")) {
      PARSE_ARGUMENT(distance.Weight());
    }
    else if (OPTION("-distance-maximum") || OPTION("-distance-max") || OPTION("-distance-max-depth")) {
      PARSE_ARGUMENT(distance.MaxDistance());
    }
    else if (OPTION("-distance-maximum-threshold") || OPTION("-distance-max-threshold")) {
      PARSE_ARGUMENT(distance.MaxThreshold());
    }
    else if (OPTION("-distance-minimum-threshold") || OPTION("-distance-min-threshold") || OPTION("-distance-threshold")) {
      PARSE_ARGUMENT(distance.MinThreshold());
    }
    else if (OPTION("-distance-averaging")) {
      PARSE_ARGUMENTS(int, distance_navgs);
    }
    else if (OPTION("-distance-smoothing")) {
      PARSE_ARGUMENT(distance.DistanceSmoothing());
    }
    else if (OPTION("-distance-measure")) {
      PARSE_ARGUMENT(distance.DistanceMeasure());
    }
    else HANDLE_BOOLEAN_OPTION("distance-hole-filling", distance.FillInHoles());
    else if (OPTION("-balloon-inflation") || OPTION("-balloon")) {
      PARSE_ARGUMENT(farg);
      balloon.Weight(farg);
      balloon.DeflateSurface(false);
    }
    else if (OPTION("-balloon-deflation")) {
      PARSE_ARGUMENT(farg);
      balloon.Weight(farg);
      balloon.DeflateSurface(true);
    }
    else if (OPTION("-balloon-min")) {
      PARSE_ARGUMENT(farg);
      balloon.LowerIntensity(farg);
    }
    else if (OPTION("-balloon-max")) {
      PARSE_ARGUMENT(farg);
      balloon.UpperIntensity(farg);
    }
    else if (OPTION("-balloon-range")) {
      PARSE_ARGUMENT(farg);
      balloon.LowerIntensity(farg);
      PARSE_ARGUMENT(farg);
      balloon.UpperIntensity(farg);
    }
    else if (OPTION("-balloon-radius")) {
      PARSE_ARGUMENT(farg);
      balloon.Radius(farg);
    }
    else if (OPTION("-balloon-sigma")) {
      PARSE_ARGUMENT(farg);
      balloon.LowerIntensitySigma(farg);
      balloon.UpperIntensitySigma(farg);
    }
    else if (OPTION("-balloon-lower-sigma")) {
      PARSE_ARGUMENT(farg);
      balloon.LowerIntensitySigma(farg);
    }
    else if (OPTION("-balloon-upper-sigma")) {
      PARSE_ARGUMENT(farg);
      balloon.UpperIntensitySigma(farg);
    }
    else if (OPTION("-balloon-averaging")) {
      PARSE_ARGUMENTS(int, balloon_navgs);
    }
    else if (OPTION("-balloon-mask")) {
      balloon_mask_name = ARGUMENT;
    }
    else if (OPTION("-edges")) {
      PARSE_ARGUMENT(farg);
      edges.Weight(farg);
    }
    else if (OPTION("-edge-distance")) {
      PARSE_ARGUMENT(farg);
      dedges.Weight(farg);
    }
    else if (OPTION("-edge-distance-type")) {
      PARSE_ARGUMENT(dedges.EdgeType());
    }
    else if (OPTION("-edge-distance-padding")) {
      PARSE_ARGUMENT(dedges.Padding());
    }
    else if (OPTION("-edge-distance-white-matter-window") ||
             OPTION("-edge-distance-wm-window")) {
      PARSE_ARGUMENT(dedges.WhiteMatterWindowWidth());
    }
    else if (OPTION("-edge-distance-grey-matter-window") ||
             OPTION("-edge-distance-gm-window")) {
      PARSE_ARGUMENT(dedges.GreyMatterWindowWidth());
    }
    else if (OPTION("-edge-distance-min-intensity")) {
      PARSE_ARGUMENT(dedges.MinIntensity());
    }
    else if (OPTION("-edge-distance-max-intensity")) {
      PARSE_ARGUMENT(dedges.MaxIntensity());
    }
    else if (OPTION("-edge-distance-min-gradient")) {
      PARSE_ARGUMENT(dedges.MinGradient());
    }
    else if (OPTION("-edge-distance-max-depth") || OPTION("-edge-distance-maximum")) {
      PARSE_ARGUMENT(dedges.MaxDistance());
    }
    else if (OPTION("-edge-distance-threshold")) {
      PARSE_ARGUMENT(dedges.DistanceThreshold());
    }
    else if (OPTION("-edge-distance-median")) {
      PARSE_ARGUMENT(dedges.MedianFilterRadius());
    }
    else if (OPTION("-edge-distance-smoothing")) {
      PARSE_ARGUMENT(dedges.DistanceSmoothing());
    }
    else if (OPTION("-edge-distance-averaging")) {
      PARSE_ARGUMENTS(int, dedges_navgs);
    }
    else if (OPTION("-inflation")) {
      PARSE_ARGUMENT(farg);
      inflation.Name("Inflation");
      inflation.Weight(farg);
    }
    else if (OPTION("-normal") || OPTION("-normal-force") || OPTION("-nforce")) {
      PARSE_ARGUMENT(nforce.Weight());
    }
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Internal forces
    unknown_option = false;
    if (OPTION("-bending-energy")) {
      PARSE_ARGUMENT(farg);
      dofbending.Weight(farg);
    }
    else if (OPTION("-neighborhood") || OPTION("-neighbourhood")) {
      model.NeighborhoodRadius(atoi(ARGUMENT));
    }
    else if (OPTION("-spring") || OPTION("-bending")) {
      PARSE_ARGUMENT(farg);
      spring.Weight(farg);
    }
    else if (OPTION("-normal-spring") || OPTION("-nspring")) {
      PARSE_ARGUMENT(farg);
      spring.InwardNormalWeight(farg);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(farg);
        spring.OutwardNormalWeight(farg);
      } else {
        spring.OutwardNormalWeight(spring.InwardNormalWeight());
      }
    }
    else if (OPTION("-tangential-spring") || OPTION("-tspring")) {
      PARSE_ARGUMENT(farg);
      spring.TangentialWeight(farg);
    }
    else if (OPTION("-normalized-spring")) {
      PARSE_ARGUMENT(farg);
      normspring.Weight(farg);
    }
    else if (OPTION("-curvature")) {
      PARSE_ARGUMENT(farg);
      curvature.Weight(farg);
    }
    else if (OPTION("-gauss-curvature") || OPTION("-gaussian-curvature") || OPTION("-gcurvature")) {
      PARSE_ARGUMENT(gcurvature.Weight());
    }
    else if (OPTION("-gauss-curvature-inside") || OPTION("-gaussian-curvature-inside") || OPTION("-gcurvature-inside")) {
      PARSE_ARGUMENT(gcurvature.WeightInside());
    }
    else if (OPTION("-gauss-curvature-outside") || OPTION("-gaussian-curvature-outside") || OPTION("-gcurvature-outside")) {
      PARSE_ARGUMENT(gcurvature.WeightOutside());
    }
    else if (OPTION("-gauss-curvature-action") || OPTION("-gaussian-curvature-action") || OPTION("-gcurvature-action")) {
      PARSE_ARGUMENT(gcurvature.NegativeGaussCurvatureAction());
      gcurvature.PositiveGaussCurvatureAction(gcurvature.NegativeGaussCurvatureAction());
    }
    else if (OPTION("-negative-gauss-curvature-action") || OPTION("-negative-gaussian-curvature-action") || OPTION("-negative-gcurvature-action")) {
      PARSE_ARGUMENT(gcurvature.NegativeGaussCurvatureAction());
    }
    else if (OPTION("-positive-gauss-curvature-action") || OPTION("-positive-gaussian-curvature-action") || OPTION("-positive-gcurvature-action")) {
      PARSE_ARGUMENT(gcurvature.PositiveGaussCurvatureAction());
    }
    else if (OPTION("-gauss-curvature-minimum") || OPTION("-gaussian-curvature-minimum") || OPTION("-gcurvature-minimum")) {
      PARSE_ARGUMENT(gcurvature.MinGaussCurvature());
    }
    else if (OPTION("-gauss-curvature-maximum") || OPTION("-gaussian-curvature-maximum") || OPTION("-gcurvature-maximum")) {
      PARSE_ARGUMENT(gcurvature.MaxGaussCurvature());
    }
    else if (OPTION("-mean-curvature") || OPTION("-mcurvature")) {
      PARSE_ARGUMENT(mcurvature.Weight());
    }
    else if (OPTION("-mean-curvature-inside") || OPTION("-mcurvature-inside")) {
      PARSE_ARGUMENT(mcurvature.WeightInside());
    }
    else if (OPTION("-mean-curvature-outside") || OPTION("-mcurvature-outside")) {
      PARSE_ARGUMENT(mcurvature.WeightOutside());
    }
    else if (OPTION("-max-curvature") || OPTION("-maximum-curvature")) {
      PARSE_ARGUMENT(pcurvature.Weight());
    }
    else if (OPTION("-max-curvature-inside") || OPTION("-maximum-curvature-inside")) {
      PARSE_ARGUMENT(pcurvature.WeightInside());
    }
    else if (OPTION("-max-curvature-outside") || OPTION("-maximum-curvature-outside")) {
      PARSE_ARGUMENT(pcurvature.WeightOutside());
    }
    else if (OPTION("-max-curvature-threshold") || OPTION("-maximum-curvature-threshold")) {
      PARSE_ARGUMENT(pcurvature.Threshold());
    }
    else if (OPTION("-quadratic-curvature") || OPTION("-qcurvature")) {
      PARSE_ARGUMENT(qcurvature.Weight());
    }
    else if (OPTION("-quadratic-curvature-inside") || OPTION("-qcurvature-inside")) {
      PARSE_ARGUMENT(qcurvature.WeightInside());
    }
    else if (OPTION("-quadratic-curvature-outside") || OPTION("-qcurvature-outside")) {
      PARSE_ARGUMENT(qcurvature.WeightOutside());
    }
    else if (OPTION("-distortion")) {
      PARSE_ARGUMENT(farg);
      distortion.Weight(farg);
    }
    else if (OPTION("-stretching")) {
      PARSE_ARGUMENT(farg);
      stretching.Weight(farg);
    }
    else if (OPTION("-stretching-rest-length")) {
      const char *arg = ARGUMENT;
      if (strcmp(arg, "avg") == 0) {
        stretching.RestLength(-1.);
        stretching.UseCurrentAverageLength(false);
      } else if (strcmp(arg, "curavg") == 0) {
        stretching.RestLength(-1.);
        stretching.UseCurrentAverageLength(true);
      } else if (FromString(arg, farg)) {
        stretching.RestLength(farg);
        stretching.UseCurrentAverageLength(false);
      } else {
        FatalError("Invalid -stretching-rest-length argument: " << arg);
      }
    }
    else if (OPTION("-repulsion")) {
      PARSE_ARGUMENT(repulsion.Weight());
    }
    else if (OPTION("-repulsion-radius")) {
      PARSE_ARGUMENT(repulsion.FrontfaceRadius());
      repulsion.BackfaceRadius(repulsion.FrontfaceRadius());
    }
    else if (OPTION("-repulsion-distance")) {
      PARSE_ARGUMENT(repulsion.FrontfaceRadius());
    }
    else if (OPTION("-repulsion-width")) {
      PARSE_ARGUMENT(repulsion.BackfaceRadius());
    }
    else if (OPTION("-collision")) {
      PARSE_ARGUMENT(farg);
      collision.Weight(farg);
    }
    else if (OPTION("-collision-radius")) {
      PARSE_ARGUMENT(collision.MinDistance());
    }
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Stopping criteria
    unknown_option = false;
    if (OPTION("-inflation-error")) {
      PARSE_ARGUMENT(farg);
      inflation_error.Threshold(farg);
    }
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Iterative local remeshing
    unknown_option = false;
    if (OPTION("-remesh")) {
      PARSE_ARGUMENT(iarg);
      model.RemeshInterval(iarg);
    }
    else if (OPTION("-remesh-adaptively")) {
      model.RemeshAdaptively(true);
    }
    else if (OPTION("-min-edge-length") || OPTION("-min-edgelength") || OPTION("-minedgelength")) {
      PARSE_ARGUMENTS(double, min_edge_length);
    }
    else if (OPTION("-max-edge-length") || OPTION("-max-edgelength") || OPTION("-maxedgelength")) {
      PARSE_ARGUMENTS(double, max_edge_length);
    }
    else if (OPTION("-min-angle") || OPTION("-minangle")) {
      PARSE_ARGUMENTS(double, min_edge_angle);
    }
    else if (OPTION("-max-angle") || OPTION("-maxangle")) {
      PARSE_ARGUMENTS(double, max_edge_angle);
    }
    else if (OPTION("-triangle-inversion")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(model.AllowTriangleInversion());
      else model.AllowTriangleInversion(true);
    }
    else if (OPTION("-notriangle-inversion")) {
      model.AllowTriangleInversion(false);
    }
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Iterative low-pass filtering
    unknown_option = false;
    if (OPTION("-lowpass")) {
      PARSE_ARGUMENT(iarg);
      model.LowPassInterval(iarg);
    }
    else if (OPTION("-lowpass-iterations")) {
      PARSE_ARGUMENT(iarg);
      model.LowPassIterations(iarg);
    }
    else if (OPTION("-lowpass-band")) {
      PARSE_ARGUMENT(farg);
      model.LowPassBand(farg);
    }
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Non-self-intersection / collision detection
    unknown_option = false;
    if (OPTION("-non-self-intersection")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(barg);
      else barg = true;
      model.HardNonSelfIntersection(barg);
    }
    else if (OPTION("-nointersection")) {
      model.HardNonSelfIntersection(true);
    }
    else if (OPTION("-min-distance") || OPTION("-mindistance") || OPTION("-mind")) {
      PARSE_ARGUMENTS(double, min_distance);
    }
    else if (OPTION("-min-width") || OPTION("-minwidth") || OPTION("-minw")) {
      PARSE_ARGUMENTS(double, min_width);
    }
    else if (OPTION("-max-collision-angle")) {
      PARSE_ARGUMENT(farg);
      model.MaxCollisionAngle(farg);
      collision.MaxAngle(farg);
    }
    else if (OPTION("-fast-collision-test")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(barg);
      else barg = true;
      model.FastCollisionTest(barg);
    }
    else if (OPTION("-nofast-collision-test")) {
      model.FastCollisionTest(false);
    }
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Output format
    unknown_option = false;
    if (OPTION("-center-output"))  center_output  = true;
    else if (OPTION("-match-area"))     match_area     = true;
    else if (OPTION("-match-sampling")) match_sampling = true;
    else if (OPTION("-track")) {
      if (HAS_ARGUMENT) track_name = ARGUMENT;
      else              track_name = "NormalDisplacement";
      track_zero_mean  = false;
      track_unit_var   = false;
      track_use_median = false;
    }
    else if (OPTION("-notrack")) track_name = nullptr;
    else if (OPTION("-track-zvalues")) {
      if (HAS_ARGUMENT)     track_name = ARGUMENT;
      else if (!track_name) track_name = "NormalDisplacementZValues";
      track_zero_mean  = true;
      track_unit_var   = true;
      track_use_median = false;
    }
    else if (OPTION("-track-zero-median-zvalues")) {
      if (HAS_ARGUMENT)     track_name = ARGUMENT;
      else if (!track_name) track_name = "NormalDisplacementZeroMedianZValues";
      track_zero_mean  = true;
      track_unit_var   = true;
      track_use_median = true;
    }
    else if (OPTION("-track-zero-mean")) {
      if (HAS_ARGUMENT)     track_name = ARGUMENT;
      else if (!track_name) track_name = "NormalDisplacementZeroMean";
      track_zero_mean  = true;
      track_use_median = false;
    }
    else if (OPTION("-track-zero-median")) {
      if (HAS_ARGUMENT)     track_name = ARGUMENT;
      else if (!track_name) track_name = "NormalDisplacementZeroMedian";
      track_zero_mean  = true;
      track_use_median = true;
    }
    else if (OPTION("-track-unit-variance")) {
      if (HAS_ARGUMENT)     track_name = ARGUMENT;
      else if (!track_name) track_name = "NormalDisplacementUnitVariance";
      track_zero_mean  = false;
      track_use_median = false;
      track_unit_var   = true;
    }
    else if (OPTION("-track-without-momentum")) {
      Insert(params, "Exclude momentum from tracked normal displacement", true);
    }
    else HANDLE_BOOLEAN_OPTION("save-status", save_status);
    else if (OPTION("-level-prefix") || OPTION("-levelprefix")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(level_prefix);
      else level_prefix = true;
    }
    else if (OPTION("-nolevel-prefix") || OPTION("-nolevelprefix")) {
      level_prefix = false;
    }
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Debugging and other common/advanced options
    unknown_option = false;
    if (OPTION("-par")) {
      const char *name  = ARGUMENT;
      const char *value = ARGUMENT;
      Insert(params, name, value);
    }
    else if (OPTION("-debug-prefix") || OPTION("-debugprefix")) {
      debug_prefix = ARGUMENT;
    }
    else if (OPTION("-debug-interval") || OPTION("-debuginterval")) {
      PARSE_ARGUMENT(iarg);
      debugger.Interval(iarg);
    }
    else HANDLE_POINTSETIO_OPTION(output_fopt);
    else {
      unknown_option = true;
    }
    if (!unknown_option) continue;

    // Common or unknown option
    HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (debug   < 0) debug   = 0;
  if (verbose < 0) verbose = 0;

  if (!image_name) {
    if (dedges.Weight() != 0.) {
      image_name = t2w_image_name;
    } else if (t1w_image_name && t2w_image_name) {
      FatalError("Not both T1-w and T2-w images used, specify only one or use -image option");
    } else {
      image_name = (t1w_image_name ? t1w_image_name : t2w_image_name);
    }
  }

  double nspring = .5 * (spring.InwardNormalWeight() + spring.OutwardNormalWeight());
  if (spring.Weight() == .0) { // no -spring, but -nspring and/or -tspring
    spring.Weight(nspring + spring.TangentialWeight());
  } else if (nspring + spring.TangentialWeight() == .0) {
    // no -nspring and -tspring, but -spring
    spring.InwardNormalWeight(.5);
    spring.OutwardNormalWeight(.5);
    spring.TangentialWeight(.5);
  }
  if (nspring + spring.TangentialWeight() <= .0) {
    spring.Weight(.0);
  }
  if ((balloon.Weight() || edges.Weight()) && !image_name) {
    FatalError("Input -image required by external forces!");
  }
  if (distance.Weight() && !dmap_name) {
    FatalError("Input -distance-map required by external -distance forces!");
  }

  string dmagnitude;
  if      (distance.Weight() > 0.) dmagnitude = distance.Name() + " magnitude";
  else if (dedges  .Weight() > 0.) dmagnitude = dedges  .Name() + " magnitude";
  if (!dmagnitude.empty()) {
    if (qcurvature.WeightInside() > 0. || qcurvature.WeightOutside() > 0.) {
      qcurvature.ExternalMagnitudeArrayName(dmagnitude);
      qcurvature.WeightMinimum(qcurvature.Weight());
      qcurvature.Weight(1.);
    }
    if (mcurvature.WeightInside() > 0. || mcurvature.WeightOutside() > 0.) {
      mcurvature.ExternalMagnitudeArrayName(dmagnitude);
      mcurvature.WeightMinimum(mcurvature.Weight());
      mcurvature.Weight(1.);
    }
    if (gcurvature.WeightInside() > 0. || gcurvature.WeightOutside() > 0.) {
      gcurvature.ExternalMagnitudeArrayName(dmagnitude);
      gcurvature.WeightMinimum(gcurvature.Weight());
      gcurvature.Weight(1.);
    }
    if (pcurvature.WeightInside() > 0. || pcurvature.WeightOutside() > 0.) {
      pcurvature.ExternalMagnitudeArrayName(dmagnitude);
      pcurvature.WeightMinimum(pcurvature.Weight());
      pcurvature.Weight(1.);
    }
  }

  // Common image attributes
  const bool force_update = true; // named variable for better readability
  ImageAttributes attr;

  // Read input image
  RegisteredImage::InputImageType input_image;
  BinaryImage image_mask;
  RegisteredImage image;
  if (image_name) {
    input_image.Read(image_name);
    attr = input_image.Attributes();
    input_image.PutBackgroundValueAsDouble(padding, true);
    if (mask_name) {
      image_mask.Read(mask_name);
      if (!image_mask.Attributes().EqualInSpace(attr)) {
        ResampleMask(image_mask, attr);
      }
      input_image.PutMask(&image_mask);
    }
    image.InputImage(&input_image);
    image.Initialize(attr);
    image.Update(true, false, false, force_update);
    image.SelfUpdate(false);
    model.Image(&image);
  }

  // Read implicit surface distance map
  RegisteredImage::InputImageType input_dmap;
  RegisteredImage dmap;
  if (dmap_name) {
    input_dmap.Read(dmap_name);
    if (attr) {
      if (!input_dmap.Attributes().EqualInSpace(attr)) {
        ResampleImage(input_dmap, attr);
      }
    } else {
      attr = input_dmap.Attributes();
    }
    dmap.InputImage(&input_dmap);
    dmap.Initialize(attr);
    dmap.Update(true, false, false, force_update);
    dmap.SelfUpdate(false);
    model.ImplicitSurface(&dmap);
  }

  // Read implicit surface distance force magnitude map
  RegisteredImage::InputImageType input_dmag;
  RegisteredImage dmag;
  if (dmag_name) {
    input_dmag.Read(dmag_name);
    if (attr && !input_dmag.Attributes().EqualInSpace(attr)) {
      ResampleImage(input_dmag, attr);
    }
    dmag.InputImage(&input_dmag);
    dmag.Initialize(input_dmag.Attributes());
    dmag.Update(true, false, false, force_update);
    dmag.SelfUpdate(false);
    distance.MagnitudeImage(&dmag);
    distance.InvertMagnitude(false);
    distance.NormalizeMagnitude(false);
  }

  // Read foreground mask of balloon force
  BinaryImage balloon_mask;
  if (balloon_mask_name) {
    balloon_mask.Read(balloon_mask_name);
    if (attr && !balloon_mask.Attributes().EqualInSpace(attr)) {
      ResampleMask(balloon_mask, attr);
    }
    balloon.ForegroundMask(&balloon_mask);
  }

  // Read tissue masks of image edge distance force
  BinaryImage wm_mask, gm_mask;
  RealImage t1w_image, cortex_dmap, vents_dmap, cerebellum_dmap;
  if (dedges.Weight() != 0.) {
    if (t1w_image_name) {
      t1w_image.Read(t1w_image_name);
      if (attr && !t1w_image.Attributes().EqualInSpace(attr)) {
        ResampleImage(t1w_image, attr);
      }
      dedges.T1WeightedImage(&t1w_image);
    }
    if (wm_mask_name) {
      wm_mask.Read(wm_mask_name);
      if (attr && !wm_mask.Attributes().EqualInSpace(attr)) {
        ResampleMask(wm_mask, attr);
      }
      dedges.WhiteMatterMask(&wm_mask);
    }
    if (gm_mask_name) {
      gm_mask.Read(gm_mask_name);
      if (attr && !gm_mask.Attributes().EqualInSpace(attr)) {
        ResampleMask(gm_mask, attr);
      }
      dedges.GreyMatterMask(&gm_mask);
    }
    if (cortex_dmap_name) {
      cortex_dmap.Read(cortex_dmap_name);
      if (attr && !cortex_dmap.Attributes().EqualInSpace(attr)) {
        ResampleImage(cortex_dmap, attr);
      }
      dedges.CorticalHullDistance(&cortex_dmap);
    }
    if (vents_dmap_name) {
      vents_dmap.Read(vents_dmap_name);
      if (attr && !vents_dmap.Attributes().EqualInSpace(attr)) {
        ResampleImage(vents_dmap, attr);
      }
      dedges.VentriclesDistance(&vents_dmap);
    }
    if (cerebellum_dmap_name) {
      cerebellum_dmap.Read(cerebellum_dmap_name);
      if (attr && !cerebellum_dmap.Attributes().EqualInSpace(attr)) {
        ResampleImage(cerebellum_dmap, attr);
      }
      dedges.CerebellumDistance(&cerebellum_dmap);
    }
  }

  // Add energy terms
  model.Add(&nforce,      false);
  model.Add(&distance,    false);
  model.Add(&balloon,     false);
  model.Add(&edges,       false);
  model.Add(&dedges,      false);
  model.Add(&spring,      false);
  model.Add(&normspring,  false);
  model.Add(&inflation,   false);
  model.Add(&curvature,   false);
  model.Add(&gcurvature,  false);
  model.Add(&mcurvature,  false);
  model.Add(&pcurvature,  false);
  model.Add(&qcurvature,  false);
  model.Add(&distortion,  false);
  model.Add(&stretching,  false);
  model.Add(&repulsion,   false);
  model.Add(&collision,   false);
  model.Add(&dofbending,  false);

  // Add stopping criteria
  if (min_energy.Threshold() > .0) {
    optimizer->AddStoppingCriterion(&min_energy);
  }
  if (min_active.Threshold() >= .0) {
    // Can be disabled with -minactive -1, otherwise, even when the stopping
    // criterion will never be fulfilled, it is needed to label nodes as passive
    optimizer->AddStoppingCriterion(&min_active);
  }
  if (!IsNaN(inflation_error.Threshold())) {
    optimizer->AddStoppingCriterion(&inflation_error);
  }

  // Set parameters
  bool ok = true;
  for (ParameterConstIterator it = params.begin(); it != params.end(); ++it) {
    if (!model     .Set(it->first.c_str(), it->second.c_str()) &&
        !optimizer->Set(it->first.c_str(), it->second.c_str())) {
      Warning("Unused/invalid parameter: " << it->first << " = " << it->second);
      ok = false;
    }
  }
  if (!ok) cout << endl;

  // Rename spring terms (after setting of parameters!)
  if (spring.Weight()) {
    if (spring.InwardNormalWeight() + spring.OutwardNormalWeight() == .0) {
      spring.Name("Tang. spring");
    }
    if (spring.TangentialWeight() == .0) {
      spring.Name("Normal spring");
    }
  }

  // Initialize deformable surface model
  const double front_repulsion_radius = repulsion.FrontfaceRadius();
  const double back_repulsion_radius  = repulsion.BackfaceRadius();
  if (front_repulsion_radius == 0.) repulsion.FrontfaceRadius(1.);

  model.GradientAveraging(0);
  model.AverageSignedGradients(signed_gradient);
  model.AverageGradientMagnitude(average_magnitude);
  model.Initialize();

  vtkPointSet  *output   = model.Output();
  vtkPointData *outputPD = output->GetPointData();
  vtkCellData  *outputCD = output->GetCellData();

  // Set output points to initial positions from previous execution
  //
  // This is necessary because some energy terms are based on the properties
  // of the original surface mesh, such as the original surface area.
  // Therefore, the input surface mesh must be identical between executions.
  // To continue optimizing a given deformable model, only replace the points
  // of the output by those of the previous output mesh (-initial argument).
  if (initial_name) {
    if (model.Transformation()) {
      FatalError("Option -initial not allowed when optimizing a parametric deformation!");
    }
    vtkSmartPointer<vtkPointSet> initial = ReadPointSet(initial_name);
    if (initial->GetNumberOfPoints() != output->GetNumberOfPoints()) {
      FatalError("Point set with initial deformed mesh points has differing number of points");
    }
    output->GetPoints()->DeepCopy(initial->GetPoints());
  }

  // Remember input point status and initialize first level status
  vtkSmartPointer<vtkDataArray> current_status = outputPD->GetArray("Status");
  vtkSmartPointer<vtkDataArray> initial_status = outputPD->GetArray("InitialStatus");
  if (current_status) {
    if (!initial_status) {
      initial_status.TakeReference(current_status->NewInstance());
      initial_status->DeepCopy(current_status);
      initial_status->SetName("InitialStatus");
      outputPD->AddArray(initial_status);
    }
  } else {
    if (!initial_status) {
      initial_status = NewVtkDataArray(VTK_UNSIGNED_CHAR, output->GetNumberOfPoints(), 1, "InitialStatus");
      initial_status->FillComponent(0, 1.);
      outputPD->AddArray(initial_status);
    }
    current_status.TakeReference(initial_status->NewInstance());
    current_status->DeepCopy(initial_status);
    current_status->SetName("Status");
    outputPD->AddArray(current_status);
  }
  if (fix_boundary) {
    vtkPolyData  * const surface = vtkPolyData::SafeDownCast(output);
    if (surface) {
      SurfaceBoundary boundary(surface);
      for (auto ptId : boundary.PointIds()) {
        initial_status->SetComponent(ptId, 0, 0.);
        current_status->SetComponent(ptId, 0, 0.);
      }
    } else {
      FatalError("Option -fix-boundary currenly only supported for surface meshes!");
    }
  }

  // Initialize optimizer
  GradientDescent *gd    = dynamic_cast<GradientDescent *>(optimizer.get());
  EulerMethod     *euler = dynamic_cast<EulerMethod     *>(optimizer.get());

  optimizer->Function(&model);
  optimizer->Initialize();

  if (gd) {
    InexactLineSearch *linesearch;
    BrentLineSearch   *brentls;
    linesearch = dynamic_cast<InexactLineSearch *>(gd->LineSearch());
    brentls    = dynamic_cast<BrentLineSearch   *>(gd->LineSearch());
    if (linesearch) {
      if (!Contains(params, "Strict total step length range")) {
        int strict = (gd->LineSearchStrategy() == LS_Brent ? 0 : 2);
        linesearch->StrictStepLengthRange(strict);
      }
      if (!Contains(params, "Maximum streak of rejected steps")) {
        int maxrejected = 1;
        if (gd->LineSearchStrategy() == LS_Brent) maxrejected = -1;
        linesearch->MaxRejectedStreak(maxrejected);
      }
      if (!Contains(params, "Minimum length of steps")) {
        linesearch->MinStepLength(.01);
      }
      if (!Contains(params, "Maximum length of steps")) {
        linesearch->MaxStepLength(1.0);
      }
      if (!Contains(params, "Maximum no. of line search iterations")) {
        linesearch->NumberOfIterations(12);
      }
    }
    if (brentls) {
      if (!Contains(params, "Brent's line search tolerance")) {
        brentls->Tolerance(.1);
      }
    }
  }

  // Add point data array to keep track of node displacment in normal direction
  // (i.e., sulcal depth measure in case of surface -inflation)
  if (track_name) {
    if (euler == NULL || model.Transformation() != NULL || !IsSurfaceMesh(output)) {
      FatalError("Option -track can currently only be used with an Euler method as -optimizer to\n"
          "       directly deform a surface mesh without a parametric transformation (no input -dof).");
    }
    vtkSmartPointer<vtkDataArray> track_array;
    track_array = outputPD->GetArray(track_name);
    if (!track_array) {
      track_array = vtkSmartPointer<vtkFloatArray>::New();
      track_array->SetName(track_name);
      track_array->SetNumberOfComponents(1);
      track_array->SetNumberOfTuples(input->GetNumberOfPoints());
      track_array->FillComponent(0, .0);
      outputPD->AddArray(track_array);
    }
    euler->NormalDisplacement(track_array);
  }

  // Deform surface until either local minimum of energy function is reached
  // or the internal and external forces of the model are in equilibrium
  const double distortion_weight = distortion.Weight();

  if (verbose > 0) {
    cout << endl;
    logger.Verbosity(verbose - 1);
    optimizer->AddObserver(logger);
  }
  if (debug > 0) {
    debugger.Prefix(debug_prefix);
    optimizer->AddObserver(debugger);
  }

  for (int level = 0; level < nlevels; ++level) {

    // Apply current distance-offset
    if (dmap_name && !dmap_offsets.empty()) {
      input_dmap.Read(dmap_name);
      input_dmap -= ParameterValue(level, nlevels, dmap_offsets, 0.);
      dmap.Update(true, false, false, force_update);
    }

    // Set number of integration steps and length of each step
    const auto dt = ParameterValue(level, nlevels, max_dt, 1.);
    optimizer->Set("Maximum length of steps", ToString(dt).c_str());
    optimizer->NumberOfSteps(ParameterValue(level, nlevels, nsteps, 100));

    optimizer->Delta(ParameterValue(level, nlevels, delta, 1e-6));

    // Set maximum node displacement at each step
    if (!max_dx.empty()) {
      const auto dx = ParameterValue(level, nlevels, max_dx, 0.);
      optimizer->Set("Normalize length of steps", "No");
      optimizer->Set("Maximum node displacement", ToString(dx).c_str());
    }
    model.MinBackfaceDistance (ParameterValue(level, nlevels, min_width,    0.));
    model.MinFrontfaceDistance(ParameterValue(level, nlevels, min_distance, 0.));

    // Set parameters of iterative remeshing step
    model.MinEdgeLength  (ParameterValue(level, nlevels, min_edge_length, 0.));
    model.MaxEdgeLength  (ParameterValue(level, nlevels, max_edge_length, inf));
    model.MinFeatureAngle(ParameterValue(level, nlevels, min_edge_angle,  180.));
    model.MaxFeatureAngle(ParameterValue(level, nlevels, max_edge_angle,  180.));

    // Set radius of repulsion force based on average edge length range
    if (front_repulsion_radius == 0. || back_repulsion_radius == 0.) {
      const auto r = model.MinEdgeLength() + .5 * (model.MaxEdgeLength() - model.MinEdgeLength());
      if (front_repulsion_radius == 0.) repulsion.FrontfaceRadius(r);
      if (back_repulsion_radius  == 0.) repulsion.BackfaceRadius (r);
    }

    // Set number of gradient averaging iterations and adjust metric distortion
    // weight for current level (cf. FreeSurfer's MRISinflateBrain function)
    const auto navg = ParameterValue(level, nlevels, navgs, 0);
    if (inflate_brain) {
      distortion.Weight(distortion_weight * sqrt(double(navg)));
      distortion.GradientAveraging(navg);
    } else {
      model.GradientAveraging(navg);
      distance.GradientAveraging(ParameterValue(level, nlevels, distance_navgs, 0));
      dedges  .GradientAveraging(ParameterValue(level, nlevels, dedges_navgs,   0));
      balloon .GradientAveraging(ParameterValue(level, nlevels, balloon_navgs,  0));
    }

    // Stopping criteria
    if (level > 0 && reset_status) {
      vtkPointSet  * const output = model.Output();
      vtkDataArray * const status = output->GetPointData()->GetArray("Status");
      if (status) {
        vtkDataArray * const initial = output->GetPointData()->GetArray("InitialStatus");
        if (initial) {
          status->CopyComponent(0, initial, 0);
        } else {
          status->FillComponent(0, 1.);
        }
      }
    }
    min_active.Threshold(ParameterValue(level, nlevels, min_active_thres, 0.));

    // Initialize optimizer
    optimizer->Initialize();

    // Debug/log output
    if (verbose > 0) {
      cout << "Level " << (level + 1) << " out of " << nlevels << "\n";
    }
    if (verbose > 1) {
      cout << "\n";
      PrintParameter(cout, "Maximum no. of steps", optimizer->NumberOfSteps());
      PrintParameter(cout, "Maximum length of steps", dt);
      PrintParameter(cout, "No. of gradient averaging steps", navg);
      if (model.RemeshInterval() > 0) {
        PrintParameter(cout, "Minimum edge length", model.MinEdgeLength());
        PrintParameter(cout, "Maximum edge length", model.MaxEdgeLength());
        PrintParameter(cout, "Minimum edge angle",  model.MinFeatureAngle());
        PrintParameter(cout, "Maximum edge angle",  model.MaxFeatureAngle());
      }
      if (inflate_brain) {
        PrintParameter(cout, "Distortion weight", distortion.Weight());
      }
      if (repulsion.Weight()) {
        PrintParameter(cout, "Repulsion frontface radius", repulsion.FrontfaceRadius());
        PrintParameter(cout, "Repulsion backface radius",  repulsion.BackfaceRadius());
      }
    }
    cout << endl;
    if (level_prefix) {
      char prefix[64];
      snprintf(prefix, 64, "%slevel_%d_", debug_prefix, level + 1);
      debugger.Prefix(prefix);
      debugger.Iteration(0);
    }

    // Perform optimization at current level
    optimizer->Run();
    if (verbose > 0) cout << endl;
  }

  optimizer->ClearObservers();

  // Remove stopping criteria to avoid their deletion by the optimizer
  optimizer->RemoveStoppingCriterion(&min_energy);
  optimizer->RemoveStoppingCriterion(&min_active);
  optimizer->RemoveStoppingCriterion(&inflation_error);

  // Get final output mesh
  output   = model.Output();
  outputPD = output->GetPointData();
  outputCD = output->GetCellData();

  // Remove data arrays used by optimizer
  if (!save_status) inputPD->RemoveArray("Status");
  for (int i = 0; i < outputPD->GetNumberOfArrays(); ++i) {
    const char *name = outputPD->GetArrayName(i);
    if (name) {
      if ((!track_name  || strcmp(name, track_name) != 0) &&
          (!save_status || strcmp(name, "Status")   != 0)) {
        if (inputPD->HasArray(name) == 0) {
          outputPD->RemoveArray(name);
          --i;
        }
      }
    }
  }
  for (int i = 0; i < outputCD->GetNumberOfArrays(); ++i) {
    const char *name = outputCD->GetArrayName(i);
    if (name && inputCD->HasArray(name) == 0) {
      outputCD->RemoveArray(name);
      --i;
    }
  }

  // Normalize sulcal depth measure tracked during inflation process
  if (track_name) {
    vtkDataArray *values = outputPD->GetArray(track_name);
    if (track_zero_mean && track_unit_var) NormalizeValues  (values, track_use_median);
    else if (track_zero_mean)              DemeanValues     (values, track_use_median);
    else if (track_unit_var)               NormalizeVariance(values);
  }

  // Resample remeshed output mesh at corresponding initial points
  if (match_sampling && model.RemeshInterval() > 0) {
    output = ResampleAtInitialPoints(input, output);
  }

  // Center output point set
  if (center_output) Center(output);

  // Scale output surface to match input area
  if (match_area) Scale(output, sqrt(Area(input) / Area(output)));

  // Write deformed output surface
  if (!WritePointSet(POSARG(2), output, output_fopt)) {
    FatalError("Failed to write output to file " << output);
  }

  return 0;
}
