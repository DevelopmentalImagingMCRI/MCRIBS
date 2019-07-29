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
#include "mirtk/PointSetUtils.h"

#include "mirtk/AsConformalAsPossibleMapper.h"
#include "mirtk/HarmonicTetrahedralMeshMapper.h"
#include "mirtk/MeshlessHarmonicVolumeMapper.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char* name)
{
  cout << "\n";
  cout << "usage: " << name << " <input> <output> [options]\n";
  cout << "\n";
  cout << "This tool computes a mapping for each point of the volume of a given input point set.\n";
  cout << "The input is either a piecewise linear complex (PLC), i.e., a tesselation of the surface,\n";
  cout << "or a tesselation of the shape's volume such as a tetrahedral mesh generated from a PLC.\n";
  cout << "The output is a volumetric map which assigns points of the volume one or more values.\n";
  cout << "The volumetric map can in general interpolate any values given on the surface of the map\n";
  cout << "domain at the interior of the volume. More common use cases are to compute a bijective\n";
  cout << "map from one volumetric shape to another with identical topology. The resulting map is a\n";
  cout << "re-parameterization of the volume of the input shape. Such parameterization can be used for\n";
  cout << "texturing, object deformation (cf. \"cage deformation\"), object morphing, and surface-\n";
  cout << "constraint image registration.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Boundary surface mesh or volumetric mesh.\n";
  cout << "  output   File path of output map. A piecewise linear map is stored as VTK file.\n";
  cout << "           Other maps are stored in a custom binary format.\n";
  cout << "\n";
  cout << "Output options:\n";
  cout << "  -acap         As-conformal-as-possible volumetric map.\n";
  cout << "  -harmonic     Harmonic volumetric map.\n";
  cout << "  -meshless     Use meshless mapping method if possible.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  PrintCommonOptions(cout);
  cout << "\n";
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of implemented volumetric mapping methods
enum MapVolumeMethod
{
  MAP_Barycentric,   ///< Compute volumetric map using generalized barycentric coordinates
  MAP_MeanValue,     ///< Compute volumetric map using mean value coordinates
  MAP_ACAP,          ///< As-conformal-as-possible (ACAP) volumetric map
  MAP_Harmonic,      ///< Compute harmonic volumetric map
  MAP_HarmonicFEM,   ///< Use FEM to compute harmonic volumetric map
  MAP_HarmonicMFS,   ///< Use MFS to compute harmonic volumetric map
  MAP_Biharmonic,    ///< Compute biharmonic volumetric map
  MAP_BiharmonicMFS, ///< Use MFS to compute biharmonic volumetric map
  MAP_Spectral       ///< Compute volumetric map using spectral coordinates
};

// -----------------------------------------------------------------------------
/// Read surface mesh or volumetric mesh
vtkSmartPointer<vtkPointSet> ReadMesh(const char *fname)
{
  vtkSmartPointer<vtkPointSet> pointset;
  if (verbose) cout << "Reading point set from " << fname << "...", cout.flush();
  const bool exit_on_failure = false;
  pointset = ReadPointSet(fname, exit_on_failure);
  if (pointset->GetNumberOfPoints() == 0) {
    if (verbose) cout << " failed" << endl;
    FatalError("Failed to read point set from " << fname << " or it has no points");
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
/// Compute volumetric map for free interior points
SharedPtr<Mapping> SolveVolumetricMap(vtkSmartPointer<vtkPointSet>  domain,
                                      vtkSmartPointer<vtkDataArray> values,
                                      vtkSmartPointer<vtkDataArray> mask,
                                      MapVolumeMethod               method,
                                      int                           niterations)
{
  SharedPtr<Mapping> map;
  if (method == MAP_Harmonic) {
    if (IsTetrahedralMesh(domain)) method = MAP_HarmonicFEM;
    else                           method = MAP_HarmonicMFS;
  } else if (method == MAP_Biharmonic) {
    method = MAP_BiharmonicMFS;
  }
  switch (method) {
    case MAP_ACAP: {
      if (verbose) cout << "Computing as-conformal-as-possible map...", cout.flush();
      AsConformalAsPossibleMapper mapper;
      mapper.InputSet(domain);
      mapper.InputMap(values);
      mapper.Run();
      map = mapper.Output();
    } break;
    case MAP_HarmonicFEM: {
      if (verbose) cout << "Computing piecewise linear harmonic map...", cout.flush();
      HarmonicTetrahedralMeshMapper mapper;
      mapper.NumberOfIterations(niterations);
      mapper.InputSet(domain);
      mapper.InputMap(values);
      mapper.InputMask(mask);
      mapper.Run();
      map = mapper.Output();
    } break;
    case MAP_HarmonicMFS: {
      if (verbose) cout << "Computing harmonic map using MFS...", cout.flush();
      MeshlessHarmonicVolumeMapper mapper;
      mapper.InputSet(domain);
      mapper.InputMap(values);
      mapper.Run();
      map = mapper.Output();
    } break;
    case MAP_BiharmonicMFS: {
      FatalError("Biharmonic mapping using MFS not implemented");
    } break;
    default:
      FatalError("Invalid volumetric map type: " << method);
  }
  if (verbose) cout << " done" << endl;
  return map;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(2);

  const char *input_name  = POSARG(1); // Input PLC or volumetric mesh
  const char *output_name = POSARG(2); // File name of output map
  const char *values_name = nullptr;   // Name of point data array with fixed point values
  const char *mask_name   = nullptr;   // Name of point data array with fixed point mask

  MapVolumeMethod method   = MAP_Harmonic;
  bool            meshless = false;
  int             niter    = 0;

  for (ALL_OPTIONS) {
    if      (OPTION("-name")) values_name = ARGUMENT;
    else if (OPTION("-mask")) mask_name   = ARGUMENT;
    // Mapping method
    else if (OPTION("-acap"))        method = MAP_ACAP;
    else if (OPTION("-barycentric")) method = MAP_Barycentric;
    else if (OPTION("-mean-value"))  method = MAP_MeanValue;
    else if (OPTION("-harmonic"))    method = MAP_Harmonic;
    else if (OPTION("-biharmonic"))  method = MAP_Biharmonic;
    else if (OPTION("-meshless"))    meshless = true;
    // Parameters of mapping method
    else if (OPTION("-max-iterations") || OPTION("-max-iter") || OPTION("-iterations") || OPTION("-iter")) {
      PARSE_ARGUMENT(niter);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }
  if (meshless) {
    if (method == MAP_Harmonic) method = MAP_HarmonicMFS;
    else                        method = MAP_HarmonicFEM;
  }

  // Read input point set
  vtkSmartPointer<vtkPointSet> domain = ReadMesh(input_name);
  vtkPointData *pd = domain->GetPointData();

  // Get boundary map
  vtkSmartPointer<vtkDataArray> values;
  const string name = ToLower(values_name);
  if      (name == "tcoords") values = pd->GetTCoords();
  else if (name == "vectors") values = pd->GetVectors();
  else if (name == "scalars") values = pd->GetScalars();
  else values = GetArrayByCaseInsensitiveName(pd, values_name);
  if (!values) {
    FatalError("Input has no point data array named " << values_name);
  }

  // Get boundary mask
  vtkSmartPointer<vtkDataArray> mask;
  if (mask_name) {
    mask = GetArrayByCaseInsensitiveName(pd, mask_name);
    if (!mask) {
      FatalError("Input has no point data array named " << mask_name);
    }
  }

  // Compute volumetric map given boundary surface map
  SharedPtr<Mapping> map(SolveVolumetricMap(domain, values, mask, method, niter));
  if (!map->Write(output_name)) {
    FatalError("Failed to write volumetric map to " << output_name);
  }

  return 0;
}
