/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/PointSetIO.h"
#include "mirtk/PiecewiseLinearMap.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char* name)
{
  cout << "\n";
  cout << "Usage: " << name << " <f1> <f2>... <g>\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Computes the composition g of the given input maps such that\n";
  cout << "\n";
  cout << "  .. math::\n";
  cout << "\n";
  cout << "     g(x) = fn o ... o f2 o f1(x)\n";
  cout << "\n";
  cout << "  The first input map, f1, must be a piecewise linear map defined at a discrete set of points.\n";
  cout << "  Values of the composite output map are computed for these input points only and the resulting\n";
  cout << "  map is also a piecewise linear map with identical discretized domain. All other maps, except\n";
  cout << "  the last one, fn, must have codomain dimension 2 or 3. Each input map can be either a surface\n";
  cout << "  map computed with the map-surface command, a volumetric map computed with the map-volume\n";
  cout << "  command, or one of the following analytic maps. Parameters of the analytic maps can be passed\n";
  cout << "  in parenthesis after the map name and must be separated by comma. See Examples below.\n";
  cout << "\n";
  cout << "  +--------------------------------+---------------------------------------------------------------------+\n";
  cout << "  | Mapping                        | Description                                                         |\n";
  cout << "  +================================+=====================================================================+\n";
  cout << "  | SquareToDisk                   | Map 2D points in a square to a disk bounded by the incircle.        |\n";
  cout << "  +--------------------------------+---------------------------------------------------------------------+\n";
  cout << "  | DiskToSquare                   | Inverse of SquareToDisk.                                            |\n";
  cout << "  +--------------------------------+---------------------------------------------------------------------+\n";
  cout << "  | StereographicProjection        | Stereographic projection from the sphere to the plane. Use argument |\n";
  cout << "  |                                | 'N' or 'S' for a projection from the north pole or south pole,      |\n";
  cout << "  |                                | respectively. The default projection is to the plane at 'z=0'.      |\n";
  cout << "  |                                | When the input points are within a disk, all points are mapped to   |\n";
  cout << "  |                                | hemisphere opposite of the projection pole.                         |\n";
  cout << "  +--------------------------------+---------------------------------------------------------------------+\n";
  cout << "  | InverseStereographicProjection | Inverse of StereographicProjection.                                 |\n";
  cout << "  +--------------------------------+---------------------------------------------------------------------+\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  f   Surface map, volumetric map, or predefined map.\n";
  cout << "  g   Composite output map.\n";
  cout << "\n";
  cout << "Examples:\n";
  cout << "  " << name << " surface-to-disk.vtp 'InverseStereographicProjection(pole=N,r=1)' surface-to-southern-hemisphere.vtp\n";
  cout << "\n";
  cout << "      Composes the given input surface map which maps each point on a surface embedded in 3D\n";
  cout << "      Euclidean space to the unit disk with an inverse stereographic projection from the north pole\n";
  cout << "      to the plane containing the disk. The resulting piecewise linear surface map projects the\n";
  cout << "      points onto the southern hemisphere of the unit sphere. The default radius, 'r', is equal\n";
  cout << "      to the maximum distance of the points from the origin along the x and y axis, respectively.\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

enum MapType
{
  MAP_Unknown,
  MAP_SquareToDisk,
  MAP_DiskToSquare,
  MAP_Stereographic,
  MAP_InvStereographic,
  MAP_Other
};

// -----------------------------------------------------------------------------
/// Split analytic map into map name and "<name>=<value>" parameter list
MapType ParseMapArgument(const char *arg, ParameterList &params)
{
  MapType type = MAP_Other;
  params.clear();
  string map = Trim(arg);
  auto pos = map.find('(');
  if (pos != map.npos && pos > 0 && map.find('(', pos + 1) == map.npos &&
      map[map.size()-1] == ')' && map.rfind(')', map.size() - 2) == map.npos) {
    auto args = Split(map.substr(pos + 1, map.size() - pos - 2), ",");
    for (auto it = args.begin(); it != args.end(); ++it) {
      auto name_value_pair = Split(*it, "=");
      if (name_value_pair.size() != 2) {
        FatalError("Failed to parse map arguments list: " << arg);
      }
      Insert(params, Trim(name_value_pair[0]), Trim(name_value_pair[1]));
    }
    map = Trim(map.substr(0, pos));
    type = MAP_Unknown;
  }
  if (map == "SquareToDisk") return MAP_SquareToDisk;
  if (map == "DiskToSquare") return MAP_DiskToSquare;
  if (map == "StereographicProjection") return MAP_Stereographic;
  if (map == "InverseStereographicProjection") return MAP_InvStereographic;
  return type;
}

// =============================================================================
// Square <-> Disk
// =============================================================================

// -----------------------------------------------------------------------------
/// Compose surface to square map with map to disk
void SquareToDisk(PiecewiseLinearMap &map, const ParameterList &params)
{
  vtkDataArray * const f = map.Values();
  Vector cdisk;
  double rdisk = .0;
  for (auto it = params.begin(); it != params.end(); ++it) {
    if (it->first == "r") {
      if (!FromString(it->second, rdisk)) {
        FatalError("SquareToDisk: Failed to parse disk radius: " << it->second);
      }
    } else if (it->first == "c") {
      if (!FromString(it->second, cdisk)) {
        FatalError("SquareToDisk: Failed to parse disk center: " << it->second);
      }
      if (cdisk.Rows() == 1) {
        cdisk.Initialize(2, cdisk(0));
      } else if (cdisk.Rows() != 2) {
        FatalError("SquareToDisk: Center point must be a single scalar or 2D vector!");
      }
    } else {
      FatalError("SquareToDisk: Unknown parameter: " << it->first);
    }
  }
  double xrange[2], yrange[2];
  f->GetRange(xrange, 0);
  f->GetRange(yrange, 1);
  double csquare[2] = { xrange[0] + (xrange[1] - xrange[0]) / 2.0,
                        yrange[0] + (yrange[1] - yrange[0]) / 2.0 };
  double rsquare = max(xrange[1] - xrange[0], yrange[1] - yrange[0]) / 2.0;
  if (cdisk.Rows() ==  0) cdisk.Initialize(2, csquare);
  if (rdisk        <= .0) rdisk = rsquare;
  double p[3], s, r;
  for (vtkIdType ptId = 0; ptId < f->GetNumberOfTuples(); ++ptId) {
    f->GetTuple(ptId, p);
    p[0] -= csquare[0], p[1] -= csquare[1];
    r = sqrt(p[0]*p[0] + p[1]*p[1]);
    if (r > .0) {
      s = (rdisk / rsquare) * (max(abs(p[0]), abs(p[1])) / r);
      p[0] = cdisk(0) + s * p[0];
      p[1] = cdisk(1) + s * p[1];
      f->SetTuple(ptId, p);
    }
  }
}

// -----------------------------------------------------------------------------
/// Compose surface to disk map with map to square
void DiskToSquare(PiecewiseLinearMap &map, const ParameterList &params)
{
  vtkDataArray * const f = map.Values();
  Vector csquare;
  double rsquare = .0;
  for (auto it = params.begin(); it != params.end(); ++it) {
    if (it->first == "r") {
      if (!FromString(it->second, rsquare)) {
        FatalError("DiskToSquare: Failed to parse radius of square incircle: " << it->second);
      }
    } else if (it->first == "c") {
      if (!FromString(it->second, csquare)) {
        FatalError("DiskToSquare: Failed to parse center point of square: " << it->second);
      }
      if (csquare.Rows() == 1) {
        csquare.Initialize(2, csquare(0));
      } else if (csquare.Rows() != 2) {
        FatalError("DiskToSquare: Center point must be a single scalar or 2D vector!");
      }
    } else {
      FatalError("DiskToSquare: Unknown parameter: " << it->first);
    }
  }
  double xrange[2], yrange[2], p[3], r, r2, cosa, sina;
  f->GetRange(xrange, 0);
  f->GetRange(yrange, 1);
  double cdisk[2] = { xrange[0] + (xrange[1] - xrange[0]) / 2.0,
                      yrange[0] + (yrange[1] - yrange[0]) / 2.0 };
  double rdisk = .5 * max(xrange[1] - xrange[0],
                          yrange[1] - yrange[0]);
  if (csquare.Rows() ==  0) csquare.Initialize(2, cdisk);
  if (rsquare        <= .0) rsquare = rdisk;
  double scale = (rsquare / rdisk);
  for (vtkIdType ptId = 0; ptId < f->GetNumberOfTuples(); ++ptId) {
    f->GetTuple(ptId, p);
    p[0] -= cdisk[0], p[1] -= cdisk[1];
    r = sqrt(p[0]*p[0] + p[1]*p[1]);
    if (r > .0) {
      r2 = r * r;
      if (fequal(p[0], p[1], 1e-6)) {
        p[0] = copysign(r, p[0]);
        p[1] = copysign(r, p[1]);
      } else if (abs(p[0]) > abs(p[1])) {
        cosa = p[0] / r;
        p[0] = copysign(r,                           p[0]);
        p[1] = copysign(sqrt(r2 / (cosa*cosa) - r2), p[1]);
      } else {
        sina = p[1] / r;
        p[0] = copysign(sqrt(r2 / (sina*sina) - r2), p[0]);
        p[1] = copysign(r,                           p[1]);
      }
      p[0] = csquare(0) + scale * p[0];
      p[1] = csquare(1) + scale * p[1];
      f->SetTuple(ptId, p);
    }
  }
}

// =============================================================================
// Spherical projection
// =============================================================================

/// Enumeration of pole used for stereographic projection
enum SpherePole { NorthPole, SouthPole };

// -----------------------------------------------------------------------------
/// Compose surface to sphere map with stereographic projection to plane
void StereographicProjection(PiecewiseLinearMap &map, const ParameterList &params)
{
  vtkDataArray *f = map.Values();
  // Check input map codomain dimension
  if (f->GetNumberOfComponents() != 3) {
    FatalError("StereographicProjection: Input must be a surface map with codomain dimension 3!");
  }
  // Parse parameters
  double c[3], p[3], z = .0, r = -1.0, s;
  SpherePole pole = NorthPole;
  for (auto it = params.begin(); it != params.end(); ++it) {
    if (it->first == "pole") {
      if (it->second == "N" || it->second == "North" || it->second == "north") {
        pole = NorthPole;
      } else if (it->second == "S" || it->second == "South" || it->second == "south") {
        pole = SouthPole;
      } else {
        FatalError("Failed to parse \"pole\" parameter value: " << it->second);
      }
    } else if (it->first == "r") {
      if (!FromString(it->second, r)) {
        FatalError("Failed to parse stereographic projection radius: " << it->second);
      }
    } else if (it->first == "z") {
      if (!FromString(it->second, z)) {
        FatalError("Failed to parse stereographic projection z value: " << it->second);
      }
    } else {
      FatalError("Invalid parameter for stereographic projection: " << it->first << "=" << it->second);
    }
  }
  // Compute center of mass
  c[0] = c[1] = c[2] = .0;
  for (vtkIdType ptId = 0; ptId < f->GetNumberOfTuples(); ++ptId) {
    f->GetTuple(ptId, p);
    c[0] += p[0], c[1] += p[1], c[2] += p[2];
  }
  c[0] /= f->GetNumberOfTuples();
  c[1] /= f->GetNumberOfTuples();
  c[2] /= f->GetNumberOfTuples();
  if (fequal(c[0], .0, 1e-6)) c[0] = .0;
  if (fequal(c[1], .0, 1e-6)) c[1] = .0;
  if (fequal(c[2], .0, 1e-6)) c[2] = .0;
  // Compute (average) radius of sphere on which points are located
  if (r < .0) {
    r = .0;
    for (vtkIdType ptId = 0; ptId < f->GetNumberOfTuples(); ++ptId) {
      f->GetTuple(ptId, p);
      p[0] -= c[0], p[1] -= c[1], p[2] -= c[2];
      r += sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    }
    r /= f->GetNumberOfTuples();
  }
  // Project sphere to plane from either the north or south pole
  const double zpole = (pole == NorthPole ? r : -r);;
  vtkSmartPointer<vtkDataArray> g;
  g.TakeReference(f->NewInstance());
  g->SetNumberOfComponents(fequal(c[2] + z, .0, 1e-6) ? 2 : 3);
  g->SetNumberOfTuples(f->GetNumberOfTuples());
  for (vtkIdType ptId = 0; ptId < f->GetNumberOfTuples(); ++ptId) {
    f->GetTuple(ptId, p);
    p[0] -= c[0], p[1] -= c[1], p[2] -= c[2];
    s = (z - zpole) / (p[2] - zpole);
    p[0] = c[0] + s * p[0];
    p[1] = c[1] + s * p[1];
    p[2] = c[2] + z;
    g->SetTuple(ptId, p);
  }
  map.Values(g);
}

// -----------------------------------------------------------------------------
/// Compose planar surface map with inverse stereographic projection to sphere
void InverseStereographicProjection(PiecewiseLinearMap &map, const ParameterList &params)
{
  vtkDataArray *f = map.Values();
  // Check input map codomain dimension
  if (f->GetNumberOfComponents() != 2) {
    FatalError("InverseStereographicProjection: Input must be a surface map with codomain dimension 2!");
  }
  // Parse parameters
  double p[3], s, r = -1.0;
  SpherePole pole = NorthPole;
  for (auto it = params.begin(); it != params.end(); ++it) {
    if (it->first == "pole") {
      if (it->second == "N" || it->second == "North" || it->second == "north") {
        pole = NorthPole;
      } else if (it->second == "S" || it->second == "South" || it->second == "south") {
        pole = SouthPole;
      } else {
        FatalError("Failed to parse \"pole\" parameter value: " << it->second);
      }
    } else if (it->first == "r") {
      if (!FromString(it->second, r)) {
        FatalError("Failed to parse stereographic projection radius: " << it->second);
      }
    } else {
      FatalError("Invalid parameter for inverse stereographic projection: " << it->first << "=" << it->second);
    }
  }
  // Determine radius of disk
  if (r < .0) {
    double xrange[2], yrange[2];
    f->GetRange(xrange, 0);
    f->GetRange(yrange, 1);
    r = max(max(abs(xrange[0]), abs(xrange[1])),
            max(abs(yrange[0]), abs(yrange[1])));
  }
  // Compose input map with inverse of stereographic projection
  const double z = (pole == NorthPole ? -r : r);
  vtkSmartPointer<vtkDataArray> g;
  g.TakeReference(f->NewInstance());
  g->SetName(f->GetName());
  g->SetNumberOfComponents(3);
  g->SetNumberOfTuples(f->GetNumberOfTuples());
  for (vtkIdType ptId = 0; ptId < f->GetNumberOfTuples(); ++ptId) {
    f->GetTuple(ptId, p);
    s = 2.0 * r * r / (p[0]*p[0] + p[1]*p[1] + r*r);
    p[0] = s * p[0];
    p[1] = s * p[1];
    p[2] = (1.0 - s) * z;
    g->SetTuple(ptId, p);
  }
  map.Values(g);
}

// =============================================================================
// Input mapping
// =============================================================================

// -----------------------------------------------------------------------------
/// Compose piecewise linear map with another Mapping
void Compose(PiecewiseLinearMap &map, const char *input_name)
{
  UniquePtr<Mapping> other(Mapping::New(input_name));
  vtkDataArray *f = map.Values();
  vtkSmartPointer<vtkDataArray> values;
  values.TakeReference(f->NewInstance());
  values->SetName(f->GetName());
  values->SetNumberOfComponents(f->GetNumberOfComponents());
  values->SetNumberOfTuples(f->GetNumberOfTuples());
  double p[3] = {.0};
  double *v = new double[f->GetNumberOfComponents()];
  for (vtkIdType i = 0; i < f->GetNumberOfTuples(); ++i) {
    f->GetTuple(i, v);
    other->Evaluate(v, p);
    values->SetTuple(static_cast<vtkIdType>(i), v);
  }
  delete[] v;
  map.Values(values);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  MapType       map;
  ParameterList params;

  REQUIRES_POSARGS(3);
  const int N = NUM_POSARGS - 1;
  for (ALL_OPTIONS) {
    HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  PiecewiseLinearMap g;
  if (!g.Read(POSARG(1))) {
    FatalError("First input map must be a piecewise linear map!");
  }
  for (int n = 2; n <= N; ++n) {
    if (g.NumberOfComponents() != 2 && g.NumberOfComponents() != 3) {
      FatalError("Codomain dimension of intermediate maps must be 2 or 3!");
    }
    map = ParseMapArgument(POSARG(n), params);
    switch (map) {
      case MAP_SquareToDisk:     SquareToDisk(g, params); break;
      case MAP_DiskToSquare:     DiskToSquare(g, params); break;
      case MAP_Stereographic:    StereographicProjection(g, params); break;
      case MAP_InvStereographic: InverseStereographicProjection(g, params); break;
      case MAP_Other:            Compose(g, POSARG(n)); break;
      default:
        FatalError("Unknown analytic mapping: " << POSARG(n));
    }
  }
  if (!g.Write(POSARG(NUM_POSARGS))) {
    FatalError("Failed to write composite map to " << POSARG(NUM_POSARGS));
  }
  return 0;
}
