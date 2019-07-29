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

#ifndef MIRTK_MinActiveStoppingCriterion_H
#define MIRTK_MinActiveStoppingCriterion_H

#include "mirtk/StoppingCriterion.h"


namespace mirtk {


/**
 * Stops mesh deformation when the ratio of active nodes drops below a threshold
 */
class MinActiveStoppingCriterion : public StoppingCriterion
{
  mirtkObjectMacro(MinActiveStoppingCriterion);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Minimum ratio of active nodes
  mirtkPublicAttributeMacro(double, Threshold);

  /// Number of iterations a node must be inactive to be labeled passive for
  /// the remainder of the iterative optimization process
  mirtkReadOnlyAttributeMacro(int, StreakOfPassiveIterations);

  /// Current ratio of active nodes
  mirtkReadOnlyAttributeMacro(double, Active);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const MinActiveStoppingCriterion &other);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  MinActiveStoppingCriterion(const ObjectiveFunction * = NULL);

  /// Copy constructor
  MinActiveStoppingCriterion(const MinActiveStoppingCriterion &);

  /// Assignment operator
  MinActiveStoppingCriterion &operator =(const MinActiveStoppingCriterion &);

  /// Create new copy of this instance
  virtual StoppingCriterion *New() const;

  /// Destructor
  virtual ~MinActiveStoppingCriterion();

  /// Initialize stopping criterion after input and parameters are set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Test stopping criterion
  ///
  /// \param[in] iter  Current number of iterations.
  /// \param[in] value Objective function value at current iteration (unused).
  /// \param[in] delta Last change of objective function parameters.
  ///
  /// \returns Whether stopping criterion is fulfilled.
  virtual bool Fulfilled(int iter, double value, const double *delta);

  // ---------------------------------------------------------------------------
  // Logging

  /// Print current stopping criterion status / value
  ///
  /// This function must be called after Fulfilled, which should update any
  /// cached values that are needed by this function to avoid a costly
  /// reevaluation of the stopping criterion.
  ///
  /// \note The printed string is expected to be considerably short and must
  ///       not end with a newline or space characters.
  virtual void Print(ostream &) const;

};


} // namespace mirtk

#endif // MIRTK_MinActiveStoppingCriterion_H
