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

#include "mirtk/Triangle.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/GenericImage.h"
#include "mirtk/GradientImageFilter.h"
#include "mirtk/PiecewiseLinearMap.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkImplicitPolyDataDistance.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char* name)
{
  cout << "\n";
  cout << "usage: " << name << " <input> [options]\n";
  cout << "\n";
  cout << "Evaluates quantitative measures of a surface map or volumetric map\n";
  cout << "such as its harmonic energy.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input   Surface map or volumetric map.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -target <file>              Point set representing the map domain. (default: none)\n";
  cout << "  -source <file>              Point set representing the map codomain. (default: none)\n";
  cout << "  -harmonic-energy [<name>]   Evaluate harmonic energy of the map. Optionally write harmonic\n";
  cout << "                              energy at each point to the named image file. (default: off)\n";
  cout << "  -lattice <file>             Lattice attributes used to discretize the map domain\n";
  cout << "                              on a regular grid are read from the given image file.\n";
  cout << "                              (default: derived from map domain)\n";
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Read surface or volumetric mesh
vtkSmartPointer<vtkPointSet> ReadMesh(const char *fname)
{
  vtkSmartPointer<vtkPointSet> pointset;
  if (verbose) cout << "Reading point set from " << fname << "...", cout.flush();
  const bool exit_on_failure = false;
  pointset = ReadPointSet(fname, exit_on_failure);
  if (pointset->GetNumberOfPoints() == 0) {
    if (verbose) cout << " failed" << endl;
    cerr << "Error: Failed to read point set from " << fname << " or it has no points" << endl;
    exit(1);
  }
  if (verbose) cout << " done" << endl;
  if (pointset->GetNumberOfCells() == 0) {
    if (verbose) cout << "Computing convex hull...", cout.flush();
    pointset = ConvexHull(pointset);
    if (verbose) cout << " done" << endl;
  }
  return pointset;
}

// -----------------------------------------------------------------------------
/// Count number of triangles with normal direction inconsistent with majority
///
/// \param[in] map Piecewise linear map from surface in 3D to 2D.
///
/// \returns Number of flipped triangles.
int NumberOfFlippedTriangles(const PiecewiseLinearMap *map)
{
  vtkDataArray * const values  = map->Values();
  vtkPolyData  * const surface = vtkPolyData::SafeDownCast(map->Domain());
  if (surface == nullptr || values->GetNumberOfComponents() != 2) return 0;

  vtkIdType npts, *pts;
  double    a[2], b[2], c[2];
  double    area;

  int n_negative = 0;
  int n_positive = 0;

  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
    surface->GetCellPoints(cellId, npts, pts);
    if (npts == 3) {
      values->GetTuple(pts[0], a);
      values->GetTuple(pts[1], b);
      values->GetTuple(pts[2], c);
      area = Triangle::DoubleSignedArea2D(a, b, c);
      if (area < 0.) {
        ++n_negative;
      } else if (area > 0.) {
        ++n_positive;
      }
    }
  }

  return min(n_positive, n_negative);
}

// -----------------------------------------------------------------------------
/// Normalize target and source models to fit into the unit box
template <class Real>
void NormalizeCoordinates(GenericImage<Real>          &map,
                          vtkSmartPointer<vtkPointSet> target,
                          vtkSmartPointer<vtkPointSet> source)
{
  double bounds[6];

  // Get homogeneous transformation matrix to map target points to unit box
  Matrix T(4, 4);
  target->GetBounds(bounds);
  T(0, 0) = 1.0 / (bounds[1] - bounds[0]);
  T(1, 1) = 1.0 / (bounds[3] - bounds[2]);
  T(2, 2) = 1.0 / (bounds[5] - bounds[4]);
  T(0, 3) = - bounds[0] / (bounds[1] - bounds[0]);
  T(1, 3) = - bounds[2] / (bounds[3] - bounds[2]);
  T(2, 3) = - bounds[4] / (bounds[5] - bounds[4]);
  T(3, 3) = 1.0;

  // Get homogeneous transformation matrix to map source points to unit box
  Matrix S(4, 4);
  source->GetBounds(bounds);
  S(0, 0) = 1.0 / (bounds[1] - bounds[0]);
  S(1, 1) = 1.0 / (bounds[3] - bounds[2]);
  S(2, 2) = 1.0 / (bounds[5] - bounds[4]);
  S(0, 3) = - bounds[0] / (bounds[1] - bounds[0]);
  S(1, 3) = - bounds[2] / (bounds[3] - bounds[2]);
  S(2, 3) = - bounds[4] / (bounds[5] - bounds[4]);
  S(3, 3) = 1.0;

  // Adjust attributes of discrete target domain lattice
  map.PutAffineMatrix(T, true);

  // Map source domain to unit box
  const int nvox = map.NumberOfSpatialVoxels();
  Real *x = map.Data(), *y = x + nvox, *z = y + nvox;
  for (int vox = 0; vox < nvox; ++vox, ++x, ++y, ++z) {
    Transform(S, *x, *y, *z);
  }
}

// -----------------------------------------------------------------------------
/// Evaluate harmonic energy of discretized volumetric map
///
/// Li et al. (2009). Meshless harmonic volumetric mapping using fundamental solution methods.
/// IEEE Transactions on Automation Science and Engineering, 6(3), 409–422.
template <class Real>
double EvaluateHarmonicEnergy(const GenericImage<Real> &map,
                              GenericImage<Real>       *energy = NULL)
{
  typedef GradientImageFilter<Real> GradientFilter;

  const ImageAttributes &lattice = map.Attributes();
  if (energy) energy->Initialize(lattice, 1);

  // Compute squared norm of gradient vectors for each component
  GenericImage<Real> squared_gradient_norm;
  GradientFilter gradient(GradientFilter::GRADIENT_DOT_PRODUCT);
  gradient.Input(&map);
  gradient.Output(&squared_gradient_norm);
  gradient.UseVoxelSize(true);
  gradient.UseOrientation(false); // R = I
  gradient.PaddingValue(numeric_limits<Real>::quiet_NaN());
  gradient.Run();

  // Evaluate harmonic energy
  const double vol  = lattice._dx * lattice._dy * lattice._dz;
  const int    nvox = lattice.NumberOfSpatialPoints();
  double value, harmonic_energy = .0;
  for (int vox = 0; vox < nvox; ++vox) {
    value = .0;
    for (int j = 0; j < lattice._t; ++j) {
      value += squared_gradient_norm(vox + j * nvox);
    }
    value *= vol;
    if (energy) energy->Put(vox, value);
    harmonic_energy += value;
  }

  return harmonic_energy;
}

// -----------------------------------------------------------------------------
/// Evaluate deformation energy of discretized volumetric map
///
/// Li et al. (2009). Meshless harmonic volumetric mapping using fundamental solution methods.
/// IEEE Transactions on Automation Science and Engineering, 6(3), 409–422.
template <class Real>
double EvaluateDeformationEnergy(const GenericImage<Real> &map,
                                 GenericImage<Real>       *energy = NULL,
                                 double lambda = .0335, double mu = .0224)
{
  typedef GradientImageFilter<Real> GradientFilter;

  if (map.T() != 3) {
    cerr << "Error: Can compute deformation energy only from 3D -> 3D volumetric map" << endl;
    exit(1);
  }

  const ImageAttributes &lattice = map.Attributes();
  if (energy) energy->Initialize(lattice, 1);

  // Compute squared norm of gradient vectors for each component
  GenericImage<Real> jac;
  GradientFilter gradient(GradientFilter::GRADIENT_VECTOR);
  gradient.Input(&map);
  gradient.Output(&jac);
  gradient.UseVoxelSize(true);
  gradient.UseOrientation(false); // R = I
  gradient.PaddingValue(numeric_limits<Real>::quiet_NaN());
  gradient.Run();

  // Evaluate deformation energy
  const double vol = lattice._dx * lattice._dy * lattice._dz;

  double stress, value, deformation_energy = .0;
  Vector3D<Real> ga, gb;
  Matrix strain(3, 3);

  for (int k = 0; k < lattice._z; ++k)
  for (int j = 0; j < lattice._y; ++j)
  for (int i = 0; i < lattice._x; ++i) {
    // Compute strain tensor, \epsilon
    for (int b = 0; b < 3; ++b) {
      gb._x = jac(i, j, k, b    ); // dq_x / dp_b
      gb._y = jac(i, j, k, b + 3); // dq_y / dp_b
      gb._z = jac(i, j, k, b + 6); // dq_z / dp_b
      for (int a = 0; a < 3; ++a) {
        ga._x = jac(i, j, k, a    ); // dq_x / dp_a
        ga._y = jac(i, j, k, a + 3); // dq_y / dp_a
        ga._z = jac(i, j, k, a + 6); // dq_z / dp_a
        strain(a, b) = ga.DotProduct(gb);
        if (a == b) strain(a, b) -= 1.0;
      }
    }
    // Compute elastic potential, \eta
    value = .0;
    for (int b = 0; b < 3; ++b)
    for (int a = 0; a < 3; ++a) {
      stress = 2.0 * mu * strain(a, b);
      if (a == b) {
        stress += lambda * strain(0, 0);
        stress += lambda * strain(1, 1);
        stress += lambda * strain(2, 2);
      }
      value += stress * strain(a, b);
    }
    value *= .5;
    // Multiply by cube volume
    value *= vol;
    // Add to total deformation energy
    if (energy) energy->Put(i, j, k, value);
    deformation_energy += value;
  }

  return deformation_energy;
}

// -----------------------------------------------------------------------------
/// Count number of lattice points mapped outside the output domain
template <class Real>
int NumberOfPointsOutside(const GenericImage<Real>     &map,
                          vtkSmartPointer<vtkPointSet>  source,
                          RealImage                    *dfield = NULL)
{
  const int nvox = map.NumberOfSpatialVoxels();

  int    n = 0;
  double p[3], d;

  vtkSmartPointer<vtkPolyData> surface = DataSetSurface(source);

  vtkSmartPointer<vtkImplicitPolyDataDistance> dist;
  dist = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
  dist->SetInput(surface);

  if (dfield) dfield->Initialize(map.Attributes(), 1);

  const Real *x = map.Data(), *y = x + nvox, *z = y + nvox;
  for (int vox = 0; vox < nvox; ++vox, ++x, ++y, ++z) {
    p[0] = *x, p[1] = *y, p[2] = *z;
    d = dist->EvaluateFunction(p);
    if (dfield) dfield->Put(vox, d);
    if (d > 0) ++n;
  }

  return n;
}

// -----------------------------------------------------------------------------
/// Compute distance of each mapped point to the output domain boundary
vtkSmartPointer<vtkDataSet> DistanceField(const PiecewiseLinearMap    *map,
                                          vtkSmartPointer<vtkPointSet> source)
{
  vtkSmartPointer<vtkPolyData> surface = DataSetSurface(source);

  vtkSmartPointer<vtkImplicitPolyDataDistance> dist;
  dist = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
  dist->SetInput(surface);

  vtkSmartPointer<vtkDataArray> darray = vtkSmartPointer<vtkFloatArray>::New();
  darray->SetName("Distance");
  darray->SetNumberOfComponents(1);
  darray->SetNumberOfTuples(map->Domain()->GetNumberOfPoints());

  vtkSmartPointer<vtkDataSet> dfield;
  dfield.TakeReference(map->Domain()->NewInstance());
  dfield->ShallowCopy(map->Domain());
  dfield->GetPointData()->Initialize();
  dfield->GetPointData()->SetScalars(darray);

  double p[3], d;
  for (vtkIdType ptId = 0; ptId < dfield->GetNumberOfPoints(); ++ptId) {
    map->Values()->GetTuple(ptId, p);
    d = dist->EvaluateFunction(p);
    darray->SetTuple(ptId, &d);
  }

  return dfield;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  double harmonic_energy;
  double deformation_energy;
  int    noutside;

  REQUIRES_POSARGS(1);

  const char *input_name  = POSARG(1);
  const char *output_name = (NUM_POSARGS >= 2 ? POSARG(2) : nullptr);

  const char *target_name  = nullptr;
  const char *source_name  = nullptr;
  const char *lattice_name = nullptr;

  bool eval_harmonic_energy    = false;
  bool eval_deformation_energy = false;
  bool eval_outside            = false;

  const char *harmonic_energy_name    = nullptr;
  const char *deformation_energy_name = nullptr;
  const char *outside_name            = nullptr;
  const char *distance_name           = nullptr;

  for (ALL_OPTIONS) {
    if      (OPTION("-target") || OPTION("-domain"))   target_name = ARGUMENT;
    else if (OPTION("-source") || OPTION("-codomain")) source_name = ARGUMENT;
    else if (OPTION("-lattice")) lattice_name  = ARGUMENT;
    else if (OPTION("-harmonic-energy")) {
      eval_harmonic_energy = true;
      if (HAS_ARGUMENT) harmonic_energy_name = ARGUMENT;
    }
    else if (OPTION("-deformation-energy")) {
      eval_deformation_energy = true;
      if (HAS_ARGUMENT) deformation_energy_name = ARGUMENT;
    }
    else if (OPTION("-outside")) {
      eval_outside = true;
      if (HAS_ARGUMENT) outside_name = ARGUMENT;
    }
    else if (OPTION("-distance")) {
      distance_name = ARGUMENT;
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (eval_outside && !source_name) {
    FatalError("Input -source required by -outside option");
  }

  // Read input point sets
  vtkSmartPointer<vtkPointSet> target, source;
  if (source_name) source = ReadMesh(source_name);
  if (target_name) target = ReadMesh(target_name);

  // Read input map
  if (verbose) cout << "Read map from " << input_name << "...", cout.flush();
  UniquePtr<Mapping> map(Mapping::New(input_name));
  PiecewiseLinearMap *dmap = dynamic_cast<PiecewiseLinearMap *>(map.get());
  if (verbose) cout << " done" << endl;

  if (target && output_name) {
    double p[3];

    // Evaluate volumetric map at points of target point set
    if (verbose) cout << "Evaluate map...", cout.flush();
    vtkSmartPointer<vtkDataArray> discrete_map = vtkSmartPointer<vtkFloatArray>::New();
    discrete_map->SetName("Map");
    discrete_map->SetNumberOfComponents(map->NumberOfComponents());
    discrete_map->SetNumberOfTuples(target->GetNumberOfPoints());
    double *v = new double[map->NumberOfComponents()];
    for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
      target->GetPoint(ptId, p);
      if (!map->Evaluate(v, p)) {
        cerr << "Warning: Map undefined at point ("
             << p[0] << ", " << p[1] << ", " << p[2]
             << ") with ID " << ptId << endl;
      }
      discrete_map->SetTuple(ptId, v);
    }
    delete[] v;
    if (discrete_map->GetNumberOfComponents() == 1) {
      target->GetPointData()->SetScalars(discrete_map);
    } else if (discrete_map->GetNumberOfComponents() == 3) {
      target->GetPointData()->SetTCoords(discrete_map);
    } else {
      target->GetPointData()->AddArray(discrete_map);
    }
    if (verbose) cout << " done" << endl;
  }

  // Write output point set
  if (!WritePointSet(output_name, target)) {
    FatalError("Failed to write point set to " << output_name);
  }

  // Evaluate volumetric map at lattice points
  RealImage discrete_map;
  if (eval_harmonic_energy || eval_deformation_energy || eval_outside || (distance_name && !dmap)) {
    if (verbose) cout << "Discretize volumetric map", cout.flush();
    if (lattice_name) {
      discrete_map.Read(lattice_name);
      discrete_map.Initialize(discrete_map.Attributes(), map->NumberOfComponents());
    } else {
      discrete_map.Initialize(map->Attributes(128, 128, 128), map->NumberOfComponents());
    }
    if (verbose > 1) {
      cout << " (N = " << discrete_map.X()
           <<    " x " << discrete_map.Y()
           <<    " x " << discrete_map.Z() << ")";
    }
    cout << "...";
    cout.flush();
    map->Evaluate(discrete_map, 0, target);
    if (debug) discrete_map.Write("discrete_map.nii.gz");
    if (verbose) cout << " done" << endl;
  }

  // Evaluate distance of mapped points to output domain boundary
  if (distance_name) {
    vtkSmartPointer<vtkPointSet> dfield;
    if (dmap) {
      vtkSmartPointer<vtkDataSet> dists = DistanceField(dmap, source);
      dfield = vtkPointSet::SafeDownCast(dists);
    }
    if (dfield) {
      if (!WritePointSet(distance_name, dfield)) {
        cerr << "Error: Failed to write mapped distance field to " << distance_name << endl;
        exit(1);
      }
    } else {
      // TODO: Modify NumberOfPointsOutside function
      cerr << "Error: Can compute -distance field currently only for tetrahedral volumetric maps" << endl;
      exit(1);
    }
  }

  // Evaluate error of volumetric map at boundary
//  double rms_boundary_error, min_boundary_error, max_boundary_error;
//  if (boundary_map_name) {
//    if (verbose) {
//      cout << "Evaluate boundary map error...";
//      if (verbose > 1) cout << "\n";
//      cout.flush();
//    }
//    vtkDataArray *boundary_map = GetBoundaryMap(target, boundary_map_name);
//    if (map->NumberOfComponents() != boundary_map->GetNumberOfComponents()) {
//      cerr << "Error: Boundary map and volumetric map have differing output domoin dimension" << endl;
//      exit(1);
//    }
//    int    num = 0;
//    double p[3], error, sum = .0;
//    min_boundary_error = +numeric_limits<double>::infinity();
//    max_boundary_error = -numeric_limits<double>::infinity();
//    double *b = new double[map->NumberOfComponents()];
//    double *v = new double[map->NumberOfComponents()];
//    for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
//      target->GetPoint(ptId, p);
//      if (!map->Evaluate(v, p)) {
//        cerr << "Warning: Volumetric map undefined at point ("
//               << p[0] << ", " << p[1] << ", " << p[2]
//               << ") with ID " << ptId << endl;
//        continue;
//      }
//      error = .0;
//      for (int j = 0; j < map->NumberOfComponents(); ++j) {
//        error += pow(v[j] - boundary_map->GetComponent(ptId, j), 2);
//      }
//      if (verbose > 1) {
//        cout << "Point " << setw(7) << (ptId + 1)
//             << ": Boundary map error = " << sqrt(error) << endl;
//      }
//      if (error < min_boundary_error) min_boundary_error = error;
//      if (error > max_boundary_error) max_boundary_error = error;
//      sum += error;
//      num += 1;
//    }
//    min_boundary_error = sqrt(min_boundary_error);
//    max_boundary_error = sqrt(max_boundary_error);
//    rms_boundary_error = (num > 0 ? sqrt(sum / num) : .0);
//    delete[] v;
//    delete[] b;
//    if (verbose) {
//      if (verbose > 1) cout << "Evaluate boundary map error...";
//      cout << " done";
//      if (verbose > 1) cout << "\n";
//      cout << endl;
//    }
//  }

  // Determine how many points are mapped outside the output domain
  if (eval_outside || (distance_name && !dmap)) {
    if (verbose) cout << "Evaluate distance of mapped points to codomain boundary...", cout.flush();
    if (outside_name || distance_name) {
      RealImage dfield;
      noutside = NumberOfPointsOutside(discrete_map, source, &dfield);
      if (outside_name ) dfield.Write(outside_name);
      if (distance_name) dfield.Write(distance_name);
    } else {
      noutside = NumberOfPointsOutside(discrete_map, source);
    }
    if (verbose) cout << " done" << endl;
  }

  // Normalize models to unit box before evaluating energy measures
  if ((eval_harmonic_energy || eval_deformation_energy) && target && source) {
    NormalizeCoordinates(discrete_map, target, source);
  }

  // Evaluate harmonic energy
  if (eval_harmonic_energy) {
    if (verbose) cout << "Evaluate harmonic energy...", cout.flush();
    if (harmonic_energy_name) {
      RealImage energy;
      harmonic_energy = EvaluateHarmonicEnergy(discrete_map, &energy);
      energy.Write(harmonic_energy_name);
    } else {
      harmonic_energy = EvaluateHarmonicEnergy(discrete_map);
    }
    if (verbose) cout << " done" << endl;
  }

  // Evaluate deformation energy
  if (eval_deformation_energy) {
    if (verbose) cout << "Evaluate deformation energy...", cout.flush();
    if (deformation_energy_name) {
      RealImage energy;
      deformation_energy = EvaluateDeformationEnergy(discrete_map, &energy);
      energy.Write(deformation_energy_name);
    } else {
      deformation_energy = EvaluateDeformationEnergy(discrete_map);
    }
    if (verbose) cout << " done" << endl;
  }

//  if (boundary_map_name) {
//    cout << "Boundary RMS error     = " << rms_boundary_error << endl;
//    cout << "Minimum boundary error = " << min_boundary_error << endl;
//    cout << "Maximum boundary error = " << max_boundary_error << endl;
//  }
  if (eval_harmonic_energy   ) cout << "Harmonic energy        = " << harmonic_energy << endl;
  if (eval_deformation_energy) cout << "Deformation energy     = " << deformation_energy << endl;
  if (eval_outside           ) cout << "No. of points outside  = " << noutside << endl;

  return 0;
}
