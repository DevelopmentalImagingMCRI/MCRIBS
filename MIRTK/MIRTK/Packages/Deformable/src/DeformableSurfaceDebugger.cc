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

#include "mirtk/DeformableSurfaceDebugger.h"

#include "mirtk/DeformableSurfaceModel.h"
#include "mirtk/EulerMethod.h"
#include "mirtk/PointSetIO.h"

#include "mirtk/CommonExport.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"


namespace mirtk {


// global "debug" flag (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int debug;


// -----------------------------------------------------------------------------
DeformableSurfaceDebugger
::DeformableSurfaceDebugger(const DeformableSurfaceModel *model, const char *prefix)
:
  _Prefix(prefix),
  _Iteration(0),
  _Model(model),
  _Interval(1)
{
}

// -----------------------------------------------------------------------------
DeformableSurfaceDebugger::~DeformableSurfaceDebugger()
{
}

// -----------------------------------------------------------------------------
void DeformableSurfaceDebugger::HandleEvent(Observable *obj, Event event, const void *data)
{
  if (_Model == NULL) return;

  const EulerMethod *euler = dynamic_cast<EulerMethod *>(obj);

  const int sz = 8;
  char suffix[sz];
  switch (event) {

    // -------------------------------------------------------------------------
    // Write intermediate results after each gradient step
    case IterationEvent:
    case IterationStartEvent: {
      ++_Iteration;
      if (_Iteration == 1) {
        _Model->WriteDataSets(_Prefix.c_str(), "_000", debug >= 3);
      }
    } break;

    case IterationEndEvent: {
      if (_Iteration % _Interval == 0) {
        snprintf(suffix, sz, "_%03d", _Iteration);
        _Model->WriteDataSets(_Prefix.c_str(), suffix, debug >= 3);
      }
      if (debug >= 2 && (_Iteration == 1 || (_Iteration % _Interval) == 0)) {
        snprintf(suffix, sz, "_%03d", _Iteration);
        _Model->WriteGradient(_Prefix.c_str(), suffix);
        if (euler) {
          const int sz = 1024;
          char fname[sz];
          snprintf(fname, sz, "%sgradient%s.vtp", _Prefix.c_str(), suffix);
          WritePointSet(fname, _Model->Output());
        }
      }
    } break;

    // -------------------------------------------------------------------------
    // Unhandled event
    default: break;
  }
}


} // namespace mirtk
