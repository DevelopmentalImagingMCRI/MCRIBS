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

#include "mirtk/StretchingForce.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/ObjectFactory.h"
#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkDataArray.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(StretchingForce);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace StretchingForceUtils {


// -----------------------------------------------------------------------------
/// Evaluate stretching penalty
struct Evaluate
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _RestLength;
  double           _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _RestLength(other._RestLength),
    _Sum(.0)
  {}

  Evaluate(const Evaluate &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _RestLength(other._RestLength),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    int    ptId1, ptId2, edgeId;
    double p1[3], p2[3], d;

    EdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); (edgeId = it.GetNextEdge(ptId1, ptId2) != -1);) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      d  = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      d -= _RestLength;
      _Sum += d * d;
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of penalty (i.e., negative internal stretching force)
struct EvaluateGradient
{
  typedef StretchingForce::GradientType Force;

  vtkPoints       *_Points;
  vtkDataArray    *_Status;
  const EdgeTable *_EdgeTable;
  double           _RestLength;
  Force           *_Gradient;

  void operator ()(const blocked_range<int> &re) const
  {
    double     p1[3], p2[3], e[3], w, d;
    const int *adjPts;
    int        numAdjPts;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      _Points->GetPoint(ptId, p1);
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
      for (int i = 0; i < numAdjPts; ++i) {
        _Points->GetPoint(adjPts[i], p2);
        vtkMath::Subtract(p2, p1, e);
        d = vtkMath::Norm(e);
        w = 2.0 * (d - _RestLength) / d;
        _Gradient[ptId] -= w * Force(e[0], e[1], e[2]);
      }
    }
  }
};


} // namespace StretchingForceUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
StretchingForce::StretchingForce(const char *name, double weight)
:
  InternalForce(name, weight),
  _RestLength(-1.),
  _AverageLength(0.),
  _UseCurrentAverageLength(true)
{
  _ParameterPrefix.push_back("Stretching ");
  _ParameterPrefix.push_back("Edge stretching ");
  _ParameterPrefix.push_back("Surface stretching ");
  _ParameterPrefix.push_back("Mesh stretching ");
  _ParameterPrefix.push_back("Surface mesh stretching ");
}

// -----------------------------------------------------------------------------
void StretchingForce::CopyAttributes(const StretchingForce &other)
{
  _RestLength              = other._RestLength;
  _AverageLength           = other._AverageLength;
  _UseCurrentAverageLength = other._UseCurrentAverageLength;
}

// -----------------------------------------------------------------------------
StretchingForce::StretchingForce(const StretchingForce &other)
:
  InternalForce(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
StretchingForce &StretchingForce::operator =(const StretchingForce &other)
{
  if (this != &other) {
    InternalForce::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
StretchingForce::~StretchingForce()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool StretchingForce::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Rest length")         == 0 ||
      strcmp(param, "Edge length")         == 0 ||
      strcmp(param, "Rest edge length")    == 0 ||
      strcmp(param, "Average edge length") == 0 ||
      strcmp(param, "Average length")      == 0) {
    string lvalue = ToLower(value);
    if (lvalue.size() > 12 && lvalue.compare(lvalue.size() - 12, 12, " edge length") == 0) {
      lvalue = lvalue.substr(0, lvalue.size() - 12);
    }
    if (lvalue == "average" || lvalue == "avg" || lvalue == "mean") {
      _RestLength = -1;
      return true;
    } else {
      return FromString(value, _RestLength);
    }
  }
  return InternalForce::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList StretchingForce::Parameter() const
{
  ParameterList params = InternalForce::Parameter();
  string rest_length;
  if (_RestLength < 0) rest_length = "Average";
  else                 rest_length = ToString(_RestLength);
  InsertWithPrefix(params, "Rest edge length", rest_length);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void StretchingForce::Initialize()
{
  // Initialize base class
  InternalForce::Initialize();
  if (_NumberOfPoints == 0) return;

  // Initialize this class
  StretchingForce::Init();
}

// -----------------------------------------------------------------------------
void StretchingForce::Reinitialize()
{
  // Reinitialize base class
  InternalForce::Reinitialize();

  // Reinitialize this class
  StretchingForce::Init();
}

// -----------------------------------------------------------------------------
void StretchingForce::Init()
{
  if (_RestLength < .0) {
    if (!_UseCurrentAverageLength) {
      _AverageLength = AverageEdgeLength(_PointSet->Points(), *_PointSet->Edges());
    }
  } else {
    _AverageLength = _RestLength;
  }
}

// -----------------------------------------------------------------------------
void StretchingForce::Update(bool gradient)
{
  // Update base class
  InternalForce::Update(gradient);

  // Update average edge length
  if (_RestLength < .0 && _UseCurrentAverageLength) {
    _AverageLength = AverageEdgeLength(_PointSet->Points(), *_PointSet->Edges());
  }
}

// -----------------------------------------------------------------------------
double StretchingForce::Evaluate()
{
  if (_PointSet->NumberOfEdges() == 0) return .0;
  MIRTK_START_TIMING();
  StretchingForceUtils::Evaluate eval;
  eval._Points     = _PointSet->Points();
  eval._EdgeTable  = _PointSet->Edges();
  eval._RestLength = _AverageLength;
  parallel_reduce(blocked_range<int>(0, _PointSet->NumberOfEdges()), eval);
  MIRTK_DEBUG_TIMING(3, "evaluation of stretching penalty");
  return eval._Sum / _PointSet->NumberOfEdges();
}

// -----------------------------------------------------------------------------
void StretchingForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  MIRTK_START_TIMING();
  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  StretchingForceUtils::EvaluateGradient eval;
  eval._Points     = _PointSet->Points();
  eval._Status     = _PointSet->Status();
  eval._EdgeTable  = _PointSet->Edges();
  eval._RestLength = _AverageLength;
  eval._Gradient   = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  for (int i = 0; i < _NumberOfPoints; ++i) {
    if (_PointSet->Edges()->NumberOfAdjacentPoints(i) > 0) {
      _Gradient[i] /= _PointSet->Edges()->NumberOfAdjacentPoints(i);
    }
  }

  InternalForce::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
  MIRTK_DEBUG_TIMING(3, "evaluation of stretching force");
}


} // namespace mirtk
