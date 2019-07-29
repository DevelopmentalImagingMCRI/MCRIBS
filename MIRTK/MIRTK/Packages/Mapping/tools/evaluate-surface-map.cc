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

#include "mirtk/EdgeTable.h"
#include "mirtk/Triangle.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PiecewiseLinearMap.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"


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
  cout << "Evaluates quantitative quality measures of a surface map.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input   Surface map.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Write "<name> = <value>" result to output stream
template <class TValue>
void Print(const char *name, const TValue &value)
{
  const streamsize w = cout.width(40);
  cout << left << name << setw(0) << " = " << value << endl;
  cout.width(w);
}

// -----------------------------------------------------------------------------
/// Get copy of surface mesh which discretizes the domain of the piecewise linear map
vtkSmartPointer<vtkPolyData> CopyMapSurface(const PiecewiseLinearMap *map)
{
  if (map->Domain()->GetDataObjectType() != VTK_POLY_DATA) {
    FatalError("Input map is piecewise-linear, but map domain is not a surface?!");
  }
  vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
  surface->ShallowCopy(map->Domain());
  surface->GetPointData()->Initialize();
  surface->GetCellData() ->Initialize();
  surface->BuildLinks();
  return surface;
}

// -----------------------------------------------------------------------------
/// Replace surface mesh points by parametric texture coordinates
vtkSmartPointer<vtkPoints> MapSurfacePoints(const PiecewiseLinearMap *map)
{
  if (map->NumberOfComponents() != 2) {
    FatalError("Surface map codomain dimension must be 2 or 3!");
  }
  vtkPolyData *surface = vtkPolyData::SafeDownCast(map->Domain());
  if (surface == nullptr) {
    FatalError("Surface map domain must be surface mesh (i.e., vtkPolyData)!");
  }
  double u[3] = {0.};
  vtkSmartPointer<vtkPoints> points;
  points.TakeReference(surface->GetPoints()->NewInstance());
  points->SetNumberOfPoints(surface->GetNumberOfPoints());
  for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
    map->GetValue(static_cast<int>(ptId), u);
    points->SetPoint(ptId, u);
  }
  return points;
}

// -----------------------------------------------------------------------------
double MappedArea(const PiecewiseLinearMap *map)
{
  vtkPolyData * const surface = vtkPolyData::SafeDownCast(map->Domain());
  if (surface == nullptr) return 0.;
  vtkSmartPointer<vtkPolyData> mapped;
  mapped.TakeReference(surface->NewInstance());
  mapped->ShallowCopy(surface);
  mapped->SetPoints(MapSurfacePoints(map));
  return Area(mapped);
}

// -----------------------------------------------------------------------------
/// Compute signed area of mapped triangles
vtkSmartPointer<vtkDataArray>
ComputeSignedAreaOfMappedPlanarTriangles(const PiecewiseLinearMap *map)
{
  vtkDataArray * const values  = map->Values();
  vtkPolyData  * const surface = vtkPolyData::SafeDownCast(map->Domain());
  if (surface == nullptr || values->GetNumberOfComponents() != 2) {
    FatalError("Surface map must be piecewise linear 2D parameterization!");
  }

  vtkSmartPointer<vtkDataArray> area;
  area = vtkSmartPointer<vtkFloatArray>::New();
  area->SetName("SignedArea");
  area->SetNumberOfComponents(1);
  area->SetNumberOfTuples(surface->GetNumberOfCells());

  vtkIdType npts, *pts;
  double    a[2], b[2], c[2], A;

  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
    surface->GetCellPoints(cellId, npts, pts);
    if (npts != 3) FatalError("Surface mesh must have triangular faces!");
    values->GetTuple(pts[0], a);
    values->GetTuple(pts[1], b);
    values->GetTuple(pts[2], c);
    A = Triangle::SignedArea2D(a, b, c);
    if (abs(A) < 1e-9) A = 0.;
    area->SetComponent(cellId, 0, A);
  }

  return area;
}

// -----------------------------------------------------------------------------
/// Count number of mapped triangles in x-y plane with inconsistent normal direction
///
/// \param[in] map   Piecewise linear map from surface in 3D to 2D plane.
/// \param[in] areas Previously computed signed areas of mapped triangles.
///
/// \returns Number of flipped triangles.
int NumberOfFlippedPlanarTriangles(const PiecewiseLinearMap     *map,
                                   vtkSmartPointer<vtkDataArray> areas = nullptr)
{
  if (areas == nullptr) {
    areas = ComputeSignedAreaOfMappedPlanarTriangles(map);
  }

  int n_negative = 0;
  int n_positive = 0;

  double area;
  for (vtkIdType cellId = 0; cellId < areas->GetNumberOfTuples(); ++cellId) {
    area = areas->GetComponent(cellId, 0);
    if      (area < 0.) ++n_negative;
    else if (area > 0.) ++n_positive;
  }

  return min(n_positive, n_negative);
}

// -----------------------------------------------------------------------------
/// Count number of mapped triangles with zero area
///
/// \param[in] map   Piecewise linear map from surface in 3D to 2D plane.
/// \param[in] areas Previously computed (signed) areas of mapped triangles.
///
/// \returns Number of degenerated triangles.
int NumberOfDegeneratedPlanarTriangles(const PiecewiseLinearMap     *map,
                                       vtkSmartPointer<vtkDataArray> areas = nullptr)
{
  if (areas == nullptr) {
    areas = ComputeSignedAreaOfMappedPlanarTriangles(map);
  }

  int n = 0;

  double area;
  for (vtkIdType cellId = 0; cellId < areas->GetNumberOfTuples(); ++cellId) {
    area = areas->GetComponent(cellId, 0);
    if (area == 0.) ++n;
  }

  return n;
}

// -----------------------------------------------------------------------------
/// Evaluate edge-length distortion
double EdgeLengthDistortion(const PiecewiseLinearMap *map,
                            vtkPointData *pd = nullptr, vtkCellData *cd = nullptr)
{
  vtkPolyData * const surface = vtkPolyData::SafeDownCast(map->Domain());
  if (surface == nullptr) {
    FatalError("Surface map domain must be surface mesh (i.e., vtkPolyData)");
  }

  vtkDataArray * const u_value = map->Values();
  if (u_value->GetNumberOfComponents() < 2 ||
      u_value->GetNumberOfComponents() > 3) {
    FatalError("Surface map must have codomain dimension 2 or 3!");
  }


  double avg_distortion = 0.;
  vtkSmartPointer<vtkDataArray> pa, ca;

  if (pd) {
    pa = vtkSmartPointer<vtkFloatArray>::New();
    pa->SetName("EdgeLengthDistortion");
    pa->SetNumberOfComponents(1);
    pa->SetNumberOfTuples(surface->GetNumberOfPoints());
    pa->FillComponent(0, 0.);
    pd->AddArray(pa);
  }
  if (cd) {
    ca = vtkSmartPointer<vtkFloatArray>::New();
    ca->SetName("EdgeLengthDistortion");
    ca->SetNumberOfComponents(1);
    ca->SetNumberOfTuples(surface->GetNumberOfCells());
    ca->FillComponent(0, 0.);
    cd->AddArray(ca);
  }

  EdgeTable edgeTable(surface);
  if (edgeTable.NumberOfEdges() > 0) {
    const double eps = 1e-12;

    int       i;
    double    p1[3], p2[3], u1[3] = {0.}, u2[3] = {0.}, scale;
    vtkIdType ptId1, ptId2, cellId;
    vtkSmartPointer<vtkIdList> cells1, cells2;
    if (cd) {
      cells1 = vtkSmartPointer<vtkIdList>::New();
      cells2 = vtkSmartPointer<vtkIdList>::New();
    }

    Vector d1(edgeTable.NumberOfEdges());
    Vector d2(edgeTable.NumberOfEdges());
    EdgeIterator edgeIt(edgeTable);
    for (edgeIt.InitTraversal(); (i = edgeIt.GetNextEdge(ptId1, ptId2)) != -1;) {
      surface->GetPoint(ptId1, p1);
      surface->GetPoint(ptId2, p2);
      u_value->GetTuple(ptId1, u1);
      u_value->GetTuple(ptId2, u2);
      d1(i) = sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) + eps;
      d2(i) = sqrt(vtkMath::Distance2BetweenPoints(u1, u2)) + eps;
    }

    const double norm = d2.Sum() / d1.Sum();
    for (edgeIt.InitTraversal(); (i = edgeIt.GetNextEdge(ptId1, ptId2)) != -1;) {
      scale = norm * d1(i) / d2(i);
      if (pa) {
        pa->SetComponent(ptId1, 0, pa->GetComponent(ptId1, 0) + scale);
        pa->SetComponent(ptId2, 0, pa->GetComponent(ptId2, 0) + scale);
      }
      if (ca) {
        surface->GetPointCells(ptId1, cells1);
        surface->GetPointCells(ptId2, cells2);
        cells1->IntersectWith(cells2);
        for (vtkIdType i = 0; i < cells1->GetNumberOfIds(); ++i) {
          cellId = cells1->GetId(i);
          ca->SetComponent(cellId, 0, ca->GetComponent(cellId, 0) + scale);
        }
      }
      avg_distortion += pow(scale - 1., 2);
    }
    if (pa) {
      int n;
      for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
        n = edgeTable.NumberOfAdjacentPoints(ptId);
        if (n > 0) pa->SetComponent(ptId, 0, pa->GetComponent(ptId, 0) / n);
      }
    }
    if (ca) {
      vtkIdType npts, *pts;
      for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
        surface->GetCellPoints(cellId, npts, pts);
        if (npts > 0) ca->SetComponent(cellId, 0, ca->GetComponent(cellId, 0) / npts);
      }
    }
    avg_distortion /= edgeTable.NumberOfEdges();
  }

  return sqrt(avg_distortion);
}

// -----------------------------------------------------------------------------
/// Evaluate triangle area distortion
double TriangleAreaDistortion(const PiecewiseLinearMap *map,
                              vtkPointData *pd = nullptr, vtkCellData *cd = nullptr)
{
  vtkPolyData * const surface = vtkPolyData::SafeDownCast(map->Domain());
  if (surface == nullptr) {
    FatalError("Surface map domain must be surface mesh (i.e., vtkPolyData)");
  }

  vtkDataArray * const u_value = map->Values();
  if (u_value->GetNumberOfComponents() < 2 ||
      u_value->GetNumberOfComponents() > 3) {
    FatalError("Surface map must have codomain dimension 2 or 3!");
  }

  double avg_distortion = 0.;
  vtkSmartPointer<vtkDataArray> pa, ca;

  if (pd) {
    pa = vtkSmartPointer<vtkFloatArray>::New();
    pa->SetName("AreaDistortion");
    pa->SetNumberOfComponents(1);
    pa->SetNumberOfTuples(surface->GetNumberOfPoints());
    pa->FillComponent(0, 0.);
    pd->AddArray(pa);
  }
  if (cd) {
    ca = vtkSmartPointer<vtkFloatArray>::New();
    ca->SetName("AreaDistortion");
    ca->SetNumberOfComponents(1);
    ca->SetNumberOfTuples(surface->GetNumberOfCells());
    ca->FillComponent(0, 0.);
    cd->AddArray(ca);
  }

  if (surface->GetNumberOfCells()) {
    const double eps = 1e-12;

    double    p1[3], p2[3], p3[3];
    double    u1[3] = {0.}, u2[3] = {0.}, u3[3] = {0.};
    double    scale;
    vtkIdType npts, *pts;

    Vector A1(surface->GetNumberOfCells());
    Vector A2(surface->GetNumberOfCells());
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
      surface->GetCellPoints(cellId, npts, pts);
      if (npts != 3) FatalError("Map domain must be triangulated!");
      surface->GetPoint(pts[0], p1);
      surface->GetPoint(pts[1], p2);
      surface->GetPoint(pts[2], p3);
      u_value->GetTuple(pts[0], u1);
      u_value->GetTuple(pts[1], u2);
      u_value->GetTuple(pts[2], u3);
      A1(static_cast<int>(cellId)) = Triangle::DoubleArea(p1, p2, p3) + eps;
      A2(static_cast<int>(cellId)) = Triangle::DoubleArea(u1, u2, u3) + eps;
    }

    const double norm = A1.Sum() / A2.Sum();
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
      scale = norm * A2(static_cast<int>(cellId)) / A1(static_cast<int>(cellId));
      if (pa) {
        surface->GetCellPoints(cellId, npts, pts);
        for (vtkIdType i = 0; i < npts; ++i) {
          pa->SetComponent(pts[i], 0, pa->GetComponent(pts[i], 0) + scale);
        }
      }
      if (ca) {
        ca->SetComponent(cellId, 0, ca->GetComponent(cellId, 0) + scale);
      }
      avg_distortion += pow(scale - 1., 2);
    }
    if (pa) {
      unsigned short n;
      vtkIdType *cells;
      for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
        surface->GetPointCells(ptId, n, cells);
        if (n > 0) pa->SetComponent(ptId, 0, pa->GetComponent(ptId, 0) / n);
      }
    }
    avg_distortion /= surface->GetNumberOfCells();
  }

  return sqrt(avg_distortion);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
#define ASSERT_IS_LINEAR_MAP() \
  if (!linmap) FatalError("Option " << OPTNAME << " requires a piecewise-linear surface map!")

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(1);

  vtkSmartPointer<vtkPolyData>  surface;
  vtkSmartPointer<vtkDataArray> array;
  vtkPointData                 *pd;
  vtkCellData                  *cd;

  UniquePtr<Mapping> map(Mapping::New(POSARG(1)));
  PiecewiseLinearMap *linmap = dynamic_cast<PiecewiseLinearMap *>(map.get());
  if (linmap) {
    surface = CopyMapSurface(linmap);
    pd = surface->GetPointData();
    cd = surface->GetCellData();
  }

  map->Initialize();

  for (ALL_OPTIONS) {
    if (OPTION("-parametric-area")) {
      ASSERT_IS_LINEAR_MAP();
      int n;
      array = ComputeSignedAreaOfMappedPlanarTriangles(linmap);
      cd->AddArray(array);
      n = NumberOfFlippedPlanarTriangles(linmap, array);
      Print("No. of flipped triangles", n);
      n = NumberOfDegeneratedPlanarTriangles(linmap, array);
      Print("No. of degenerated triangles", n);
    }
    else if (OPTION("-flipped-triangles")) {
      ASSERT_IS_LINEAR_MAP();
      int n = NumberOfFlippedPlanarTriangles(linmap, cd->GetArray("SignedArea"));
      Print("No. of flipped triangles", n);
    }
    else if (OPTION("-degenerated-triangles")) {
      ASSERT_IS_LINEAR_MAP();
      int n = NumberOfDegeneratedPlanarTriangles(linmap, cd->GetArray("SignedArea"));
      Print("No. of degenerated triangles", n);
    }
    else if (OPTION("-edge-distortion") || OPTION("-edge-length-distortion")) {
      ASSERT_IS_LINEAR_MAP();
      double distortion = EdgeLengthDistortion(linmap, nullptr, cd);
      Print("Average edge-length distortion", distortion);
    }
    else if (OPTION("-area-distortion")) {
      ASSERT_IS_LINEAR_MAP();
      double distortion = TriangleAreaDistortion(linmap, nullptr, cd);
      Print("Average area distortion", distortion);
    }
    else if (OPTION("-map-points")) {
      ASSERT_IS_LINEAR_MAP();
      surface->SetPoints(MapSurfacePoints(linmap));
    }
    else if (OPTION("-o") || OPTION("-output")) {
      ASSERT_IS_LINEAR_MAP();
      // Write surface mesh with thus far computed point/cell attributes
      const char *output_name = ARGUMENT;
      if (!WritePolyData(output_name, surface)) {
        FatalError("Failed to write surface mesh to " << output_name);
      }
      // Reset surface mesh
      surface = CopyMapSurface(linmap);
      pd = surface->GetPointData();
      cd = surface->GetCellData();
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  return 0;
}
