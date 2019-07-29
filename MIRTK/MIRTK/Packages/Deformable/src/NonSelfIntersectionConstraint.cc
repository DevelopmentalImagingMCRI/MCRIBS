/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#include "mirtk/NonSelfIntersectionConstraint.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/ObjectFactory.h"
#include "mirtk/VtkMath.h"

#include "vtkUnsignedCharArray.h"
#include "vtkGenericCell.h"
#include "vtkCellData.h"
#include "vtkPoints.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(NonSelfIntersectionConstraint);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace NonSelfIntersectionConstraintUtils {

// -----------------------------------------------------------------------------
typedef SurfaceCollisions::CollisionInfo      CollisionInfo;
typedef SurfaceCollisions::CollisionsSet      CollisionsSet;
typedef SurfaceCollisions::CollisionsIterator CollisionsIterator;
typedef SurfaceCollisions::CollisionsArray    CollisionsArray;

// -----------------------------------------------------------------------------
/// Calculate sum of edge lengths
struct SumEdgeLengths
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _Sum;

  SumEdgeLengths() : _Sum(.0) {}

  SumEdgeLengths(const SumEdgeLengths &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Sum(.0)
  {}

  void join(const SumEdgeLengths &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    vtkIdType ptId1, ptId2;
    double    p1[3], p2[3];

    EdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); it.GetNextEdge(ptId1, ptId2) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      _Sum += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate non-self-intersection penalty
struct Evaluate
{
  const CollisionsArray *_Collisions;
  double                 _MinDistance;
  double                 _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Collisions (other._Collisions),
    _MinDistance(other._MinDistance),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<size_t> &re)
  {
    for (size_t cellId = re.begin(); cellId != re.end(); ++cellId) {
      const CollisionsSet &colls = (*_Collisions)[cellId];
      for (CollisionsIterator it = colls.begin(); it != colls.end(); ++it) {
        const CollisionInfo &coll = *it;
        _Sum += _MinDistance - coll._Distance;
      }
    }
  }
};

// -----------------------------------------------------------------------------
// Evaluate non-self-intersection force
struct EvaluateGradient
{
  typedef NonSelfIntersectionConstraint::GradientType Force;

  vtkSmartPointer<vtkPolyData> _DataSet;
  const CollisionsArray       *_Collisions;
  double                       _MinDistance;
  Force                       *_Gradient;
  int                         *_Count;

  EvaluateGradient() : _Gradient(NULL), _Count(NULL) {}

  EvaluateGradient(const EvaluateGradient &other, split)
  :
    _DataSet    (other._DataSet),
    _Collisions (other._Collisions),
    _MinDistance(other._MinDistance)
  {
    _Gradient = CAllocate<Force>(_DataSet->GetNumberOfPoints());
    _Count    = CAllocate<int>  (_DataSet->GetNumberOfPoints());
  }

  void join(EvaluateGradient &other)
  {
    for (vtkIdType i = 0; i < _DataSet->GetNumberOfPoints(); ++i) {
      _Gradient[i] += other._Gradient[i];
      _Count   [i] += other._Count   [i];
    }
    Deallocate(other._Gradient);
    Deallocate(other._Count);
  }

  void operator ()(const blocked_range<size_t> &re)
  {
    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    double    w, f[3];
    vtkIdType ptId;

#if 0
    int    subId;
    double dist2, pcoords[3], *weights;
    weights = new double[_DataSet->GetMaxCellSize()];
#endif

    for (size_t cellId = re.begin(); cellId != re.end(); ++cellId) {
      const CollisionsSet &colls = (*_Collisions)[cellId];
      for (CollisionsIterator it = colls.begin(); it != colls.end(); ++it) {
        const CollisionInfo &coll = *it;
        if (coll._Distance > .0) {
          _DataSet->GetCell(cellId, cell);
          w = abs(_MinDistance - coll._Distance) / (_MinDistance * coll._Distance);
          f[0] = w * (coll._Point1[0] - coll._Point2[0]);
          f[1] = w * (coll._Point1[1] - coll._Point2[1]);
          f[2] = w * (coll._Point1[2] - coll._Point2[2]);
#if 0 // does not work well for coplanar triangles as long as
      // SurfaceCollisions does not return the center points of the
      // projected overlaps
          if (cell->EvaluatePosition(const_cast<double *>(coll._Point1), NULL,
                                     subId, pcoords, dist2, weights) == 1) {
            for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); ++i) {
              ptId = cell->GetPointId(i);
              _Gradient[ptId]._x -= weights[i] * f[0];
              _Gradient[ptId]._y -= weights[i] * f[1];
              _Gradient[ptId]._z -= weights[i] * f[2];
              ++_Count[ptId];
            }
          } else
#endif
          {
            for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); ++i) {
              ptId = cell->GetPointId(i);
              _Gradient[ptId] -= Force(f[0], f[1], f[2]);
              ++_Count[ptId];
            }
          }
        }
      }
    }

#if 0
    delete[] weights;
#endif
  }
};


} // namespace NonSelfIntersectionConstraintUtils
using namespace NonSelfIntersectionConstraintUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
NonSelfIntersectionConstraint
::NonSelfIntersectionConstraint(const char *name, double weight)
:
  SurfaceConstraint(name, weight),
  _FastCollisionTest(false),
  _MinDistance(-1.0),
  _MaxAngle(60.)
{
  _ParameterPrefix.push_back("Non-self-intersection ");
  _ParameterPrefix.push_back("Non-self-intersection force ");
}

// -----------------------------------------------------------------------------
void NonSelfIntersectionConstraint
::CopyAttributes(const NonSelfIntersectionConstraint &other)
{
  _FastCollisionTest = other._FastCollisionTest;
  _MinDistance       = other._MinDistance;
  _MaxAngle          = other._MaxAngle;
  _Collisions        = other._Collisions;
}

// -----------------------------------------------------------------------------
NonSelfIntersectionConstraint
::NonSelfIntersectionConstraint(const NonSelfIntersectionConstraint &other)
:
  SurfaceConstraint(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
NonSelfIntersectionConstraint &NonSelfIntersectionConstraint
::operator =(const NonSelfIntersectionConstraint &other)
{
  if (this != &other) {
    SurfaceConstraint::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
NonSelfIntersectionConstraint::~NonSelfIntersectionConstraint()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool NonSelfIntersectionConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Distance")                 == 0 ||
      strcmp(param, "Distance threshold")       == 0 ||
      strcmp(param, "Minimum distance")         == 0 ||
      strcmp(param, "Minimum surface distance") == 0) {
    return FromString(value, _MinDistance);
  }
  return SurfaceConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList NonSelfIntersectionConstraint::Parameter() const
{
  ParameterList params = SurfaceConstraint::Parameter();
  InsertWithPrefix(params, "Distance threshold", _MinDistance);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void NonSelfIntersectionConstraint::Init()
{
  AllocateCount(_NumberOfPoints);
}

// -----------------------------------------------------------------------------
void NonSelfIntersectionConstraint::Initialize()
{
  // Reset
  _NumberOfCollisions = 0;
  _Collisions.clear();

  // Initialize base class
  SurfaceConstraint::Initialize();

  // Initialize this class
  NonSelfIntersectionConstraint::Init();

  // Set minimum distance to average edge length if unset
  if (_MinDistance < .0) {
    const int nedges = _PointSet->SurfaceEdges()->NumberOfEdges();
    if (nedges > 0) {
      NonSelfIntersectionConstraintUtils::SumEdgeLengths eval;
      eval._Points    = _PointSet->SurfacePoints();
      eval._EdgeTable = _PointSet->SurfaceEdges();
      parallel_reduce(blocked_range<int>(0, nedges), eval);
      _MinDistance = eval._Sum / nedges;
    }
  }
}

// -----------------------------------------------------------------------------
void NonSelfIntersectionConstraint::Reinitialize()
{
  // Reinitialize base class
  SurfaceConstraint::Reinitialize();

  // Reinitialize this class
  NonSelfIntersectionConstraint::Init();
}

// -----------------------------------------------------------------------------
void NonSelfIntersectionConstraint::Update(bool gradient)
{
  if (!gradient && _NumberOfCollisions == 0) return;

  vtkPolyData *surface = _PointSet->Surface();

  SurfaceCollisions check;
  check.Input(surface);
  check.FastCollisionTest(_FastCollisionTest);
  check.AdjacentIntersectionTestOff();
  check.NonAdjacentIntersectionTestOff();
  check.FrontfaceCollisionTestOn();
  check.BackfaceCollisionTestOn();
  check.StoreCollisionDetailsOn();
  check.MinFrontfaceDistance(_MinDistance);
  check.MinBackfaceDistance(_MinDistance);
  check.MaxSearchRadius(10.0 * _MinDistance);
  check.MaxAngle(_MaxAngle);

  if (!gradient) {
    // For efficiency reasons, perform check only for previously collided cells
    // Update(true) is usually called by the optimizer before each line search
    vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
    mask->SetNumberOfComponents(1);
    mask->SetNumberOfTuples(surface->GetNumberOfCells());
    for (size_t cellId = 0; cellId < _Collisions.size(); ++cellId) {
      mask->SetTuple1(cellId, static_cast<double>(!_Collisions[cellId].empty()));
    }
    check.Mask(mask);
  }

  check.Run();

  _NumberOfCollisions = 0;
  _Collisions.clear();

  if (check.FoundCollisions()) {

    const double maxNormalAngleCos = cos(_MaxAngle);

    bool iscoll;
    double n1[3], n2[3];
    vtkDataArray *normals = _PointSet->SurfaceFaceNormals();
    const SurfaceCollisions::CollisionsArray &colls = check.Collisions();

    for (size_t cellId = 0; cellId < colls.size(); ++cellId) {
      CollisionsIterator coll = colls[cellId].begin();
      while (coll != colls[cellId].end()) {
        iscoll = (coll->_Distance > .0);
        if (iscoll) {
          normals->GetTuple(cellId,        n1);
          normals->GetTuple(coll->_CellId, n2);
          iscoll = (vtkMath::Dot(n1, n2) < maxNormalAngleCos);
        }
        if (iscoll) {
          if (_Collisions.empty()) _Collisions.resize(colls.size());
          _Collisions[cellId].insert(*coll);
          ++_NumberOfCollisions;
        }
        ++coll;
      }
    }
  }
}

// -----------------------------------------------------------------------------
double NonSelfIntersectionConstraint::Evaluate()
{
  if (_NumberOfCollisions == 0) return .0;

  NonSelfIntersectionConstraintUtils::Evaluate eval;
  eval._Collisions  = &_Collisions;
  eval._MinDistance = _MinDistance;
  parallel_reduce(blocked_range<size_t>(0, _Collisions.size()), eval);

  return eval._Sum / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void NonSelfIntersectionConstraint::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));
  memset(_Count,    0, _NumberOfPoints * sizeof(int));

  if (_NumberOfCollisions == 0) return;

  NonSelfIntersectionConstraintUtils::EvaluateGradient eval;
  eval._DataSet     = _PointSet->Surface();
  eval._Collisions  = &_Collisions;
  eval._MinDistance = _MinDistance;
  eval._Gradient    = _Gradient;
  eval._Count       = _Count;
  parallel_reduce(blocked_range<size_t>(0, _Collisions.size()), eval);

  for (int i = 0; i < _NumberOfPoints; ++i) {
    if (_Count[i] > 0) _Gradient[i] /= _Count[i];
  }

  SurfaceConstraint::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void NonSelfIntersectionConstraint::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  if (_NumberOfPoints == 0 && !all) return;

  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
  surface->ShallowCopy(_PointSet->Surface());

  vtkSmartPointer<vtkDataArray> coll;
  coll = vtkSmartPointer<vtkUnsignedCharArray>::New();
  coll->SetName("NumberOfCollisions");
  coll->SetNumberOfComponents(1);
  coll->SetNumberOfTuples(surface->GetNumberOfCells());
  surface->GetCellData()->AddArray(coll);

  if (_Collisions.empty()) {
    coll->FillComponent(0, .0);
  } else {
    coll->FillComponent(0, .0);
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
      coll->SetTuple1(cellId, static_cast<double>(_Collisions[cellId].size()));
    }
  }

  snprintf(fname, sz, "%ssurface%s.vtp", prefix, suffix);
  WritePolyData(fname, surface);
}


} // namespace mirtk
