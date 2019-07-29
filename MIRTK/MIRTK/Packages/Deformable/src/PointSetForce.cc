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

#include "mirtk/PointSetForce.h"

#include "mirtk/Vtk.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/MultiLevelTransformation.h"

#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkVertex.h"


namespace mirtk {


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace PointSetForceUtils {

// -----------------------------------------------------------------------------
// Typedefs
typedef PointSetForce::GradientType GradientType;

// -----------------------------------------------------------------------------
/// Perform one iteration of gradient averaging
struct AverageGradient
{
  const EdgeTable *_EdgeTable;
  GradientType    *_Input;
  GradientType    *_Output;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int n;
    const int *adjPts;

    GradientType *out = _Output + re.begin();
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, ++out) {
      (*out) = _Input[ptId];
      _EdgeTable->GetAdjacentPoints(ptId, n, adjPts);
      for (int i = 0; i < n; ++i) {
        (*out) += _Input[adjPts[i]];
      }
      (*out) /= (n + 1);
    }
  }
};

// -----------------------------------------------------------------------------
/// Perform one iteration of gradient magnitude averaging
struct AverageGradientMagnitude
{
  const EdgeTable *_EdgeTable;
  GradientType    *_Input;
  GradientType    *_Output;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int        numAdjPts;
    const int *adjPts;
    double     norm, avg_norm;

    GradientType *in  = _Input  + re.begin();
    GradientType *out = _Output + re.begin();

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, ++in, ++out) {
      norm = in->Length();
      if (norm) {
        avg_norm = norm;
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
        for (int i = 0; i < numAdjPts; ++i) {
          avg_norm += _Input[adjPts[i]].Length();
        }
        avg_norm /= (numAdjPts + 1);
        avg_norm /= norm;
        (*out) = (*in) * avg_norm;
      } else {
        (*out) = .0;
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Perform one iteration of signed gradient averaging
struct AverageSignedGradient
{
  const EdgeTable *_EdgeTable;
  GradientType    *_Input;
  GradientType    *_Output;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int        numAdjPts, n;
    const int *adjPts;

    GradientType *in  = _Input  + re.begin();
    GradientType *out = _Output + re.begin();

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, ++in, ++out) {
      (*out) = (*in), n = 1;
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
      for (int i = 0; i < numAdjPts; ++i) {
        const GradientType &adj = _Input[adjPts[i]];
        if (in->DotProduct(adj) > .0) {
          (*out) += adj, ++n;
        }
      }
      (*out) /= n;
    }
  }
};

// -----------------------------------------------------------------------------
/// Perform one iteration of signed gradient magnitude averaging
struct AverageSignedGradientMagnitude
{
  const EdgeTable *_EdgeTable;
  GradientType    *_Input;
  GradientType    *_Output;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int        numAdjPts, n;
    const int *adjPts;
    double     norm, avg_norm;

    GradientType *in  = _Input  + 3 * re.begin();
    GradientType *out = _Output + 3 * re.begin();

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, ++in, ++out) {
      norm = in->Length();
      if (norm) {
        avg_norm = norm, n = 1;
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
        for (int i = 0; i < numAdjPts; ++i) {
          const GradientType &adj = _Input[adjPts[i]];
          if (in->DotProduct(adj) > .0) {
            avg_norm += adj.Length(), ++n;
          }
        }
        avg_norm /= n;
        avg_norm /= norm;
        (*out) = (*in) * avg_norm;
      } else {
        (*out) = .0;
      }
    }
  }
};

// -----------------------------------------------------------------------------
template <class AvgFunc>
void AverageGradientVectors(GradientType *gradient, const EdgeTable *edgeTable, int niter)
{
  const int npoints = edgeTable->NumberOfPoints();
  AvgFunc avg;
  avg._Input     = gradient;
  avg._Output    = Allocate<GradientType>(npoints);
  avg._EdgeTable = edgeTable;
  blocked_range<vtkIdType> ptIds(0, npoints);
  for (int iter = 0; iter < niter; ++iter) {
    parallel_for(ptIds, avg);
    swap(avg._Input, avg._Output);
  }
  if (avg._Output == gradient) {
    delete[] avg._Input;
  } else {
    memcpy(gradient, avg._Output, npoints * sizeof(GradientType));
    delete[] avg._Output;
  }
}


} // namespace PointSetForceUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void PointSetForce::AllocateGradient(int n)
{
  if (_GradientSize < n || n <= 0) {
    Deallocate(_Gradient);
    if (n > 0) {
      _Gradient     = Allocate<GradientType>(n);
      _GradientSize = n;
    } else {
      _GradientSize = 0;
    }
  }
}

// -----------------------------------------------------------------------------
void PointSetForce::AllocateCount(int n)
{
  if (_CountSize < n || n <= 0) {
    Deallocate(_Count);
    if (n > 0) {
      _Count     = Allocate<int>(n);
      _CountSize = n;
    } else {
      _CountSize = 0;
    }
  }
}

// -----------------------------------------------------------------------------
PointSetForce::PointSetForce(const char *name, double weight)
:
  EnergyTerm(name, weight),
  _PointSet(nullptr),
  _GradientAveraging(0),
  _AverageSignedGradients(false),
  _AverageGradientMagnitude(false),
  _SurfaceForce(false),
  _Gradient(nullptr),
  _GradientSize(0),
  _Count(nullptr),
  _CountSize(0),
  _InitialUpdate(false)
{
}

// -----------------------------------------------------------------------------
void PointSetForce::CopyAttributes(const PointSetForce &other)
{
  _PointSet                 = other._PointSet;
  _GradientAveraging        = other._GradientAveraging;
  _AverageSignedGradients   = other._AverageSignedGradients;
  _AverageGradientMagnitude = other._AverageGradientMagnitude;
  _SurfaceForce             = other._SurfaceForce;
  _InitialUpdate            = other._InitialUpdate;
  AllocateGradient(other._GradientSize);
  AllocateCount(other._CountSize);
}

// -----------------------------------------------------------------------------
PointSetForce::PointSetForce(const PointSetForce &other)
:
  EnergyTerm(other),
  _Gradient(nullptr),
  _GradientSize(0),
  _Count(nullptr),
  _CountSize(0)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
PointSetForce &PointSetForce::operator =(const PointSetForce &other)
{
  if (this != &other) {
    EnergyTerm::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
PointSetForce::~PointSetForce()
{
  Deallocate(_Gradient);
  Deallocate(_Count);
}

// =============================================================================
// Point set attributes
// =============================================================================

// -----------------------------------------------------------------------------
vtkDataArray *PointSetForce::PointData(const char *name, bool optional) const
{
  vtkDataArray *data = nullptr;
  vtkPointSet  * const ps = DeformedPointSet();
  vtkPointData * const pd = ps->GetPointData();
  auto it = _PointDataName.find(name);
  if (it == _PointDataName.end()) {
    data = pd->GetArray(name);
  } else {
    data = pd->GetArray(it->second.c_str());
  }
  if (data) {
    if (data->GetNumberOfComponents() <= 0) {
      Throw(ERR_LogicError, __FUNCTION__, "Point data array has no components!");
    }
    if (data->GetNumberOfTuples() != ps->GetNumberOfPoints()) {
      Throw(ERR_LogicError, __FUNCTION__, "Point data array has invalid size!\n"
            "  This indicates that the point data array was not correctly adjusted\n"
            "  during the remeshing of the deformed ", _SurfaceForce ? "surface" : "point set",
            ". Please report\n  this bug or debug the execution in order to fix this issue.");
    }
  } else if (!optional) {
    Throw(ERR_LogicError, __FUNCTION__, _SurfaceForce ? "Surface" : "Point set",
          " has no point data array named: ", name);
  }
  return data;
}

// -----------------------------------------------------------------------------
void PointSetForce::AddPointData(const char *name, vtkSmartPointer<vtkDataArray> &data, bool global)
{
  // Remove previously added array if any
  RemovePointData(name);

  // Remove array from point set attributes to prevent duplicate additions
  vtkPointData * const pd = this->PointData();
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    if (pd->GetArray(i) == data) pd->RemoveArray(i--);
  }

  // Set unique array name
  if (global) {
    data->SetName(name);
  } else {
    string prefix = ParameterNameWithPrefix(name);
    string unique = prefix;
    bool is_unique = true;
    for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
      if (unique == pd->GetArrayName(i)) {
        is_unique = false;
        break;
      }
    }
    if (!is_unique) {
      for (int j = 1; j <= 99; ++j) {
        unique = prefix + ToString(j);
        is_unique = true;
        for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
          if (unique == pd->GetArrayName(i)) {
            is_unique = false;
            break;
          }
        }
        if (is_unique) break;
      }
      if (!is_unique) unique = prefix + "X";
    }
    data->SetName(unique.c_str());
  }

  // Add point data array
  pd->AddArray(data);
  _PointDataName[name] = data->GetName();
}

// -----------------------------------------------------------------------------
vtkDataArray *PointSetForce::AddPointData(const char *name, int c, int type, bool global)
{
  const bool optional = true;
  vtkSmartPointer<vtkDataArray> data = PointData(name, optional);
  if (!data) {
    data = NewVtkDataArray(type);
    data->SetNumberOfComponents(c);
  } else if (data->GetDataType() != type || data->GetNumberOfComponents() != c) {
    if (global) {
      Throw(ERR_LogicError, __FUNCTION__, "Mismatch of global data array type and/or number of components");
    } else {
      data = NewVtkDataArray(type);
      data->SetNumberOfComponents(c);
    }
  }
  if (_NumberOfPoints > 0) data->SetNumberOfTuples(_NumberOfPoints);
  AddPointData(name, data, global);
  return data;
}

// -----------------------------------------------------------------------------
void PointSetForce::RemovePointData(const char *name)
{
  auto it = _PointDataName.find(name);
  if (it == _PointDataName.end()) return;
  PointData()->RemoveArray(it->second.c_str());
  _PointDataName.erase(it);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPoints> PointSetForce::GetInitialPoints() const
{
  vtkSmartPointer<vtkPoints> points;
  double p[3];

  const MultiLevelTransformation *mffd;
  mffd = dynamic_cast<const MultiLevelTransformation *>(_PointSet->Transformation());

  if (_SurfaceForce) {

    if (mffd) {
      points = vtkSmartPointer<vtkPoints>::New();
      points->SetNumberOfPoints(_PointSet->NumberOfSurfacePoints());
      for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
        _PointSet->GetInputSurfacePoint(i, p);
        mffd->GlobalTransform(p[0], p[1], p[2]);
        points->SetPoint(i, p);
      }
    } else {
      points = vtkSmartPointer<vtkPoints>::New();
      points->DeepCopy(_PointSet->InputSurface()->GetPoints());
    }

  } else {

    if (mffd) {
      points = vtkSmartPointer<vtkPoints>::New();
      points->SetNumberOfPoints(_PointSet->NumberOfPoints());
      for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
        _PointSet->GetInputPoint(i, p);
        mffd->GlobalTransform(p[0], p[1], p[2]);
        points->SetPoint(i, p);
      }
    } else {
      points = vtkSmartPointer<vtkPoints>::New();
      points->DeepCopy(_PointSet->InputPointSet()->GetPoints());
    }

  }

  return points;
}

// -----------------------------------------------------------------------------
void PointSetForce::Init()
{
  // Get number of mesh vertices
  if (_SurfaceForce) {
    _NumberOfPoints = _PointSet->NumberOfSurfacePoints();
  } else {
    _NumberOfPoints = _PointSet->NumberOfPoints();
  }

  // Allocate gradient vector
  // Note: _Count must be first allocated by subclass if required
  AllocateGradient(_NumberOfPoints);
  if (_CountSize > 0) AllocateCount(_NumberOfPoints);
}

// -----------------------------------------------------------------------------
void PointSetForce::Initialize()
{
  // Initialize base class
  EnergyTerm::Initialize();

  // Indicates also that force is inactive
  _NumberOfPoints = 0;

  // Free previously allocated memory
  AllocateGradient(0);
  AllocateCount(0);

  // Check input
  if (!_PointSet) {
    cerr << "PointSetForce::Initialize: Input point set not set" << endl;
    exit(1);
  }
  if (!_PointSet->InputPointSet()) {
    cerr << "PointSetForce::Initialize: Undeformed point set not set" << endl;
    exit(1);
  }

  // Initialize this class
  PointSetForce::Init();

  // Delayed initialization upon next Update call
  _InitialUpdate = true;
}

// -----------------------------------------------------------------------------
void PointSetForce::Reinitialize()
{
  PointSetForce::Init();
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void PointSetForce::Update(bool)
{
  if (_InitialUpdate || _PointSet->Transformation()) {
    _PointSet->Update(_InitialUpdate && _PointSet->SelfUpdate());
  }
  _InitialUpdate = false;
}

// -----------------------------------------------------------------------------
void PointSetForce::EvaluateGradient(double *gradient, double, double weight)
{
  using namespace PointSetForceUtils;

  if (_NumberOfPoints > 0) {

    // Smooth node displacements such that adjacent nodes move coherently.
    // Can also be viewed as an averaging of the gradient vectors in a local
    // neighborhood. With decreasing smoothing iterations, a multi-resolution
    // optimization of the deformable surface model can be mimicked.
    if (_GradientAveraging > 0) {
      MIRTK_START_TIMING();
      SharedPtr<const EdgeTable> edgeTable = SharedEdgeTable();
      if (_AverageSignedGradients) {
        if (_AverageGradientMagnitude) {
          typedef struct AverageSignedGradientMagnitude AvgOp;
          AverageGradientVectors<AvgOp>(_Gradient, edgeTable.get(), _GradientAveraging);
        } else {
          typedef struct AverageSignedGradient AvgOp;
          AverageGradientVectors<AvgOp>(_Gradient, edgeTable.get(), _GradientAveraging);
        }
      } else {
        if (_AverageGradientMagnitude) {
          typedef struct AverageGradientMagnitude AvgOp;
          AverageGradientVectors<AvgOp>(_Gradient, edgeTable.get(), _GradientAveraging);
        } else {
          typedef struct AverageGradient AvgOp;
          AverageGradientVectors<AvgOp>(_Gradient, edgeTable.get(), _GradientAveraging);
        }
      }
      MIRTK_DEBUG_TIMING(3, "averaging of"
          << (_AverageSignedGradients ? " signed " : " ")
          << "energy gradient" << (_AverageGradientMagnitude ? " magnitude " : " ")
          << "(#iter=" << _GradientAveraging << ")");
    }

    if (_PointSet->Transformation()) {
      const double t0 = _PointSet->InputTime();
      const double t  = _PointSet->Time();
      const class PointSet &pos = _SurfaceForce ? _PointSet->InputSurfacePoints()
                                                : _PointSet->InputPoints();
      _PointSet->Transformation()->ParametricGradient(pos, _Gradient, gradient, t, t0, weight);
    } else {
      double         *g;
      vtkIdTypeArray *origPtIds = _PointSet->OriginalSurfacePointIds();
      if (_SurfaceForce && origPtIds) {
        for (int i = 0, j; i < _NumberOfPoints; ++i) {
          j = static_cast<int>(origPtIds->GetComponent(i, 0));
          g = gradient + 3 * j;
          g[0] += weight * _Gradient[i]._x;
          g[1] += weight * _Gradient[i]._y;
          g[2] += weight * _Gradient[i]._z;
        }
      } else {
        double *g = gradient;
        for (int i = 0; i < _NumberOfPoints; ++i, g += 3) {
          g[0] += weight * _Gradient[i]._x;
          g[1] += weight * _Gradient[i]._y;
          g[2] += weight * _Gradient[i]._z;
        }
      }
    }
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void PointSetForce::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  if (_NumberOfPoints == 0 && !all) return;

  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  if (_SurfaceForce) {
    snprintf(fname, sz, "%ssurface%s.vtp", prefix, suffix);
    WritePolyData(fname, _PointSet->Surface());
  } else {
    snprintf(fname, sz, "%spointset%s%s", prefix, suffix, _PointSet->DefaultExtension());
    _PointSet->Write(fname);
  }
}

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkFloatArray> ToFloatArray(const PointSetForce::GradientType *v, int n)
{
  vtkSmartPointer<vtkFloatArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetNumberOfComponents(3);
  array->SetNumberOfTuples(n);
  for (int i = 0; i < n; ++i) array->SetTuple3(i, v[i]._x, v[i]._y, v[i]._z);
  return array;
}

// -----------------------------------------------------------------------------
void PointSetForce::WriteGradient(const char *p, const char *suffix) const
{
  if (_NumberOfPoints == 0) return;

  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  vtkSmartPointer<vtkDataArray> gradient = ToFloatArray(_Gradient, _NumberOfPoints);
  gradient->SetName("gradient");

  if (_SurfaceForce) {

    vtkSmartPointer<vtkPoints> points;
    points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(_NumberOfPoints);

    vtkSmartPointer<vtkCellArray> vertices;
    vertices = vtkSmartPointer<vtkCellArray>::New();
    vertices->Allocate(_NumberOfPoints);

    double pt[3];
    for (vtkIdType i = 0; i < _NumberOfPoints; ++i) {
      _PointSet->GetInputSurfacePoint(i, pt);
      points->SetPoint(i, pt);
      vertices->InsertNextCell(1, &i);
    }

    vtkSmartPointer<vtkPolyData> output;
    output = vtkSmartPointer<vtkPolyData>::New();
    output->SetPoints(points);
    output->SetVerts(vertices);
    output->GetPointData()->AddArray(gradient);

    snprintf(fname, sz, "%sgradient%s.vtp", prefix, suffix);
    WritePolyData(fname, output);

  } else {

    snprintf(fname, sz, "%sgradient%s%s", prefix, suffix, _PointSet->DefaultExtension());
    _PointSet->Write(fname, gradient);

  }
}


} // namespace mirtk
