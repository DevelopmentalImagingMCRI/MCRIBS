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

#include "mirtk/UniformSurfaceMapper.h"                   // Tutte (1964)
#include "mirtk/ChordLengthSurfaceMapper.h"               // Kent et al. (1991), Floater (1997)
#include "mirtk/HarmonicSurfaceMapper.h"                  // Pinker & Polthier (1993)
#include "mirtk/ShapePreservingSurfaceMapper.h"           // Floater (1997)
#include "mirtk/AuthalicSurfaceMapper.h"                  // Meyer et al. (2002)
#include "mirtk/IntrinsicSurfaceMapper.h"                 // Meyer et al. (2002)
#include "mirtk/IntrinsicLeastAreaDistortionSurfaceMapper.h"   // Meyer et al. (2002)
#include "mirtk/IntrinsicLeastEdgeLengthDistortionSurfaceMapper.h"   // Meyer et al. (2002)
#include "mirtk/MeanValueSurfaceMapper.h"                 // Floater (2003)
#include "mirtk/ConformalSurfaceFlattening.h"             // Angenent (1999), Haker (2000)
#include "mirtk/LeastSquaresConformalSurfaceMapper.h"     // Levy (2002), Desbrun et al. (2002)

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
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
  cout << "This tool computes a mapping for each point on the surface of a given input shape\n";
  cout << "embedded in 3D space. The output is a (piecewise linear) function which assigns each\n";
  cout << "point on the surface of the input shape one or more values. In case of non-closed surfaces,\n";
  cout << "the output map can interpolate any values given on the boundary of the surface at the\n";
  cout << "interior points of the surface. More common use cases are to compute a bijective mapping\n";
  cout << "from one geometric shape to another geometric shape with identical topology. The resulting\n";
  cout << "map is a parameterization of the surface of the input shape. Such parameterization can be\n";
  cout << "used for texturing, object morphing, and surface registration.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Point set delineating the boundary of the map domain or name\n";
  cout << "           of primitive shape such as: \"disk\", \"square\", or \"sphere\".\n";
  cout << "  output   File path of output map. A piecewise linear map is stored as VTK file.\n";
  cout << "           Other maps are stored in a custom binary format.\n";
  cout << "\n";
  cout << "Output options:\n";
  cout << "  -barycentric   Use spring constants based on generalized barycentric coordiantes.\n";
  cout << "  -mean-value    Use spring constants based on mean value coordinates.\n";
  cout << "  -conformal     Conformal surface map or as-conformal-as-possible volumetric map.\n";
  cout << "  -harmonic      Harmonic volumetric map.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -p <n>                Exponent of harmonic energy term. When non-positive, solve for an\n";
  cout << "                        approximate harmonic surface map using a spring network. (default: 0)\n";
  cout << "  -name <string>        Name of point data array used as fixed point map.  (default: tcoords)\n";
  cout << "  -mask <string>        Name of point data array used as fixed point mask. (default: boundary)\n";
  cout << "  -max-iterations <n>   Maximum no. of linear solver iterations. (default: 1 or size of problem)\n";
  PrintCommonOptions(cout);
  cout << "\n";
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of implemented surface mapping methods
enum SurfaceMappingMethod
{
  MAP_Uniform,                      ///< Uniform edge weights, Tutte's planar graph mapping
  MAP_ChordLength,                  ///< Edge weights inverse proportional to edge length
  MAP_ShapePreserving,              ///< Floater's shape-preserving weights
  MAP_MeanValue,                    ///< Floater's mean value convex map
  MAP_Harmonic,                     ///< Pinker and Polthier's cotangent weights
                                    ///< also known as discrete conformal parameterization (DCP)
  MAP_PHarmonic,                    ///< Joshi's p-harmonic map
  MAP_Authalic,                     ///< Desbrun and Meyer's area preserving weights
                                    ///< referred to as discrete authalic parameterization (DAP)
  MAP_Intrinsic,                    ///< Meyer's intrinsic parameterization based on
                                    ///< generalized Barycentric coordinates
  MAP_IntrinsicLeastAreaDistortion, ///< Meyer's intrinsic parameterization with
                                    ///< non-linear minimization of area distortion
  MAP_IntrinsicLeastEdgeDistortion, ///< Meyer's intrinsic parameterization with
                                    ///< non-linear minimization of edge length distortion
  MAP_LeastSquaresConformal,        ///< Levy's least squares conformal map (LSCM) which
                                    ///< is identical to Meyer's discrete natural conformal
                                    ///< parameterization (DNCP)
  MAP_ConformalFlattening,          ///< Angenent and Haker's conformal map to the sphere
  MAP_Spectral,                     ///< Spectral surface map w/o boundary constraints
  MAP_Spherical                     ///< Spherical surface map w/o boundary constraints
};

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(2);

  const char *input_name        = POSARG(1);
  const char *output_name       = POSARG(2);
  const char *boundary_map_name = nullptr;

  SurfaceMappingMethod method = MAP_MeanValue;

  int    niters                = -1; // Number of iterations
  int    p_harmonic_exponent   = 2;  // Exponent of p-harmonic energy
  int    chord_length_exponent = 1;  // Weighted least squares exponent
  double intrinsic_lambda      = .5; // Conformal vs. authalic energy weight
  Array<int> selection;              // Selected (boundary) points

  for (ALL_OPTIONS) {
    // Fixed boundary map
    if (OPTION("-boundary-map")) boundary_map_name = ARGUMENT;
    else if (OPTION("-select")) {
      int i;
      do {
        PARSE_ARGUMENT(i);
        selection.push_back(i);
      } while (HAS_ARGUMENT);
    }
    // Surface mapping method
    else if (OPTION("-uniform")) {
      method = MAP_Uniform;
    }
    else if (OPTION("-chord-length")) {
      method = MAP_ChordLength;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(chord_length_exponent);
      else chord_length_exponent = 1;
    }
    else if (OPTION("-shape-preserving") || OPTION("-average-barycenter")) {
      method = MAP_ShapePreserving;
    }
    else if (OPTION("-mean-value") || OPTION("-mean-value-coordinates") || OPTION("-mvc")) {
      method = MAP_MeanValue;
    }
    else if (OPTION("-harmonic")) {
      method = MAP_Harmonic;
    }
    else if (OPTION("-p-harmonic")) {
      method = MAP_PHarmonic;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(p_harmonic_exponent);
      else p_harmonic_exponent = 2;
    }
    else if (OPTION("-authalic")) {
      method = MAP_Authalic;
    }
    else if (OPTION("-intrinsic")) {
      method = MAP_Intrinsic;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(intrinsic_lambda);
      else intrinsic_lambda = .5;
    }
    else if (OPTION("-intrinsic-least-area-distortion") ||
             OPTION("-intrinsic-min-area-distortion")) {
      method = MAP_IntrinsicLeastAreaDistortion;
    }
    else if (OPTION("-intrinsic-least-edge-length-distortion") ||
             OPTION("-intrinsic-min-edge-length-distortion") ||
             OPTION("-intrinsic-least-edge-distortion") ||
             OPTION("-intrinsic-min-edge-distortion")) {
      method = MAP_IntrinsicLeastEdgeDistortion;
    }
    else if (OPTION("-conformal-flattening")) {
      method = MAP_ConformalFlattening;
    }
    else if (OPTION("-least-squares-conformal") || OPTION("-lscm") ||
             OPTION("-natural-conformal") || OPTION("-discrete-natural-conformal") || OPTION("-dncp")) {
      method = MAP_LeastSquaresConformal;
    }
    // Linear solver parameters
    else if (OPTION("-max-iterations") || OPTION("-max-iter") || OPTION("-iterations") || OPTION("-iter")) {
      PARSE_ARGUMENT(niters);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  vtkSmartPointer<vtkPolyData>  surface;
  SharedPtr<PiecewiseLinearMap> boundary_map;
  SharedPtr<Mapping>            surface_map;

  surface = ReadPolyData(input_name);
  if (boundary_map_name) {
    boundary_map = NewShared<PiecewiseLinearMap>();
    if (!boundary_map->Read(boundary_map_name)) {
      FatalError("Failed to read boundary map from " << boundary_map_name);
    }
    boundary_map->OutsideValue(0.);
    boundary_map->Initialize();
  }

  switch (method) {
    case MAP_Uniform: {
      const char *msg = "Computing uniform surface map...";
      if (verbose) cout << msg, cout.flush();
      UniformSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Input(boundary_map);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_ChordLength: {
      const char *msg = "Computing chord length weighted surface map...";
      if (verbose) cout << msg, cout.flush();
      ChordLengthSurfaceMapper mapper(chord_length_exponent);
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Input(boundary_map);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_ShapePreserving: {
      const char *msg = "Computing shape preserving surface map...";
      if (verbose) cout << msg, cout.flush();
      ShapePreservingSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Input(boundary_map);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_MeanValue: {
      const char *msg = "Computing mean value convex map...";
      if (verbose) cout << msg, cout.flush();
      MeanValueSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Input(boundary_map);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_Harmonic: {
      const char *msg = "Computing harmonic surface map...";
      if (verbose) cout << msg, cout.flush();
      HarmonicSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Input(boundary_map);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_PHarmonic: {
      FatalError("p-harmonic mapping using finite element method (FEM) not implemented");
      if (verbose) cout << "Computing p=" << p_harmonic_exponent << " harmonic surface map...", cout.flush();
    } break;

    case MAP_Authalic: {
      const char *msg = "Computing discrete authalic surface map...";
      if (verbose) cout << msg, cout.flush();
      AuthalicSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Input(boundary_map);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_Intrinsic: {
      const char *msg = "Computing surface map using intrinsic parameterization...";
      if (verbose) cout << msg, cout.flush();
      IntrinsicSurfaceMapper mapper(intrinsic_lambda);
      if (verbose) {
        cout << "\n  Conformal energy weight      = " << mapper.Lambda();
        cout << "\n  Authalic  energy weight      = " << 1. - mapper.Lambda();
        cout.flush();
      }
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Input(boundary_map);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_IntrinsicLeastAreaDistortion: {
      const char *msg = "Computing intrinsic surface map with least area distortion...";
      if (verbose) cout << msg, cout.flush();
      IntrinsicLeastAreaDistortionSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Input(boundary_map);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_IntrinsicLeastEdgeDistortion: {
      const char *msg = "Computing intrinsic surface map with least edge length distortion...";
      if (verbose) cout << msg, cout.flush();
      IntrinsicLeastEdgeLengthDistortionSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Input(boundary_map);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_ConformalFlattening: {
      if (boundary_map) {
        FatalError("Conformal flattening requires closed genus-0 input surface mesh! Input -boundary-map makes no sense.");
      }
      const char *msg = "Computing conformal flattening...";
      if (verbose) cout << msg, cout.flush();
      ConformalSurfaceFlattening mapper;
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    case MAP_LeastSquaresConformal: {
      if (boundary_map) {
        Warning("Input -boundary-map ignored by least squares conformal mapping.");
      }
      const char *msg = "Computing least squares conformal map...";
      if (verbose) cout << msg, cout.flush();
      LeastSquaresConformalSurfaceMapper mapper;
      mapper.NumberOfIterations(niters);
      mapper.Surface(surface);
      if (selection.size() > 0) {
        mapper.AddFixedPoint(selection[0], 0., 0.);
        if (selection.size() > 1) {
          mapper.AddFixedPoint(selection[1], 1., 0.);
        }
      }
      mapper.Run();
      surface_map = mapper.Output();
      if (verbose) cout << msg, cout.flush();
    } break;

    default: {
      FatalError("Selected mapping method not implemented");
    } break;
  }

  if (!surface_map->Write(output_name)) {
    if (verbose) cout << " failed" << endl;
    FatalError("Failed to write surface map to " << output_name);
  }
  if (verbose) cout << " done" << endl;

  return 0;
}
