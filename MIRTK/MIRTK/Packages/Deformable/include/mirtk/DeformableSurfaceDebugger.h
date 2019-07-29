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

#ifndef MIRTK_DeformableSurfaceDebugger_H
#define MIRTK_DeformableSurfaceDebugger_H

#include "mirtk/Observer.h"


namespace mirtk {


// Forward declaration of observed object
class DeformableSurfaceModel;


/**
 * Writes intermediate surfaces to the current working directory
 *
 * This object observes the optimizer of the deformable surface model and must
 * therefore be attached to the respective LocalOptimizer instance
 * (typically EulerMethod or a subclass of it).
 */
class DeformableSurfaceDebugger : public Observer
{
  mirtkObjectMacro(DeformableSurfaceDebugger);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Prefix for output file names
  mirtkPublicAttributeMacro(string, Prefix);

  /// Iteration counter
  mirtkPublicAttributeMacro(int, Iteration);

  /// Reference to the deformable surface model
  mirtkPublicAggregateMacro(const DeformableSurfaceModel, Model);

  /// Write intermediate results only every n gradient steps
  mirtkPublicAttributeMacro(int, Interval);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy construction
  /// \note Intentionally not implemented.
  DeformableSurfaceDebugger(const DeformableSurfaceDebugger &);

  /// Assignment operator
  /// \note Intentionally not implemented.
  DeformableSurfaceDebugger &operator =(const DeformableSurfaceDebugger &);

public:

  /// Constructor
  DeformableSurfaceDebugger(const DeformableSurfaceModel * = NULL, const char * = "");

  /// Destructor
  ~DeformableSurfaceDebugger();

  /// Handle event and print message to output stream
  void HandleEvent(Observable *, Event, const void *);

};


} // namespace mirtk

#endif // MIRTK_DeformableSurfaceDebugger_H
