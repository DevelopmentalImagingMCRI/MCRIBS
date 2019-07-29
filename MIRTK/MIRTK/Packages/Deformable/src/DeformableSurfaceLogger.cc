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

#include "mirtk/DeformableSurfaceLogger.h"

#include "mirtk/Math.h"
#include "mirtk/DeformableSurfaceModel.h"
#include "mirtk/LocalOptimizer.h"
#include "mirtk/EulerMethod.h"
#include "mirtk/StoppingCriterion.h"
#include "mirtk/Terminal.h"

#include "mirtk/CommonExport.h"


namespace mirtk {


// Global "debug_time" flag (cf. mirtkProfiling.cc)
MIRTK_Common_EXPORT extern int debug_time;


// =============================================================================
// Auxiliaries
// =============================================================================

namespace DeformableSurfaceLoggerUtils {


// -----------------------------------------------------------------------------
ostream &PrintNumber(ostream &os, double value)
{
  const ios::fmtflags fmt = os.flags();
  if (value != .0 && (abs(value) < 1.e-5 || abs(value) >= 1.e5)) {
    os << scientific << setprecision(5) << setw(16) << value; // e-0x
  } else os << fixed << setprecision(5) << setw(12) << value << "    ";
  os.flags(fmt);
  return os;
}

// -----------------------------------------------------------------------------
ostream &PrintNormalizedNumber(ostream &os, double value)
{
  const ios::fmtflags fmt = os.flags();
  os << fixed << setprecision(5) << setw(8) << value;
  os.flags(fmt);
  return os;
}

// -----------------------------------------------------------------------------
ostream &PrintWeight(ostream &os, double weight, int nterms)
{
  const ios::fmtflags fmt = os.flags();
  const int w = iround(abs(weight) * 100.0);
  if      (w ==   0)               os << "< 1";
  else if (w == 100 && nterms > 1) os << ">99";
  else                             os << fixed << setw(3) << w;
  os.flags(fmt);
  return os;
}

// -----------------------------------------------------------------------------
ostream &PrintDelta(ostream &os, double delta, double dt)
{
  const ios::fmtflags fmt = os.flags();
  static const int    p = 4;
  static const double t = pow(0.1, p);
  const int w = (ifloor(dt) + 10) / 10 + p + 1;
  if (delta < t) {
    os << scientific << setprecision(0) << setw(w) << delta;
  } else {
    os << fixed << setprecision(p) << setw(w) << delta;
  }
  os.flags(fmt);
  return os;
}


} // namespace DeformableSurfaceLoggerUtils
using namespace DeformableSurfaceLoggerUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
DeformableSurfaceLogger::DeformableSurfaceLogger(ostream *stream)
:
  _Verbosity  (0),
  _Stream     (stream),
  _Color      (stdout_color),
  _FlushBuffer(true)
{
}

// -----------------------------------------------------------------------------
DeformableSurfaceLogger::~DeformableSurfaceLogger()
{
}

// =============================================================================
// Logging
// =============================================================================

// -----------------------------------------------------------------------------
void DeformableSurfaceLogger::HandleEvent(Observable *obj, Event event, const void *data)
{
  if (!_Stream) return;
  ostream &os = *_Stream;

  // Change/Remember stream output format
  const streamsize w = os.width(0);
  const streamsize p = os.precision(6);

  // Deformable surface model (i.e., energy function)
  LocalOptimizer *optimizer   = dynamic_cast<LocalOptimizer *>(obj);
  EulerMethod    *eulermethod = dynamic_cast<EulerMethod *>(optimizer);
  if (optimizer == NULL) return;
  DeformableSurfaceModel *model = dynamic_cast<DeformableSurfaceModel *>(optimizer->Function());
  if (model == NULL) return;

  // Typed pointers to event data (which of these is valid depends on the type of event)
  const Iteration      *iter = reinterpret_cast<const Iteration      *>(data);
  const LineSearchStep *step = reinterpret_cast<const LineSearchStep *>(data);
  const char           *msg  = reinterpret_cast<const char           *>(data);

  // Note: endl forces flush of stream buffer! Use \n instead and only flush buffer
  //       at end of this function if _FlushBuffer flag is set.

  switch (event) {

    // Start of optimization
    case StartEvent: {
      if (_Verbosity == 0) {
        if (_Color) os << xboldblack;
        os << "          Energy     ";
        string name;
        int    i;
        size_t w;
        for (i = 0; i < model->NumberOfTerms(); ++i) {
          const EnergyTerm *term = model->Term(i);
          if (term->Weight() != .0) {
            w = 16;
            if (model->NumberOfTerms() < 5 &&
                term->DivideByInitialValue()) w += 10;
            name                   = term->Name();
            if (name.empty()) name = term->NameOfClass();
            if (name.length() > w-3) {
              name = name.substr(0, w-3), name += "...";
            }
            os << string(3 + (w - name.length())/2, ' ');
            os << name;
            os << string((w - name.length() + 1)/2, ' ');
          }
        }
        if (optimizer->NumberOfStoppingCriteria() > 0 || eulermethod) {
          os << "     ";
          if (optimizer->NumberOfStoppingCriteria() == 1) {
            os << "Stopping criterion";
          } else {
            os << "Stopping criteria";
          }
        }
        if (_Color) os << xreset;
      }
      _NumberOfGradientSteps = 0;
    } break;

    // Next gradient step
    case IterationEvent:
    case IterationStartEvent: {
      ++_NumberOfGradientSteps;
      os << "\n";
      if (debug_time || _Verbosity > 0) {
        os << "Iteration " << setw(3) << left << _NumberOfGradientSteps << right << " ";
      } else {
        os << " " << setw(3) << _NumberOfGradientSteps << " ";
      }
      if (debug_time) os << "\n";
      // If this event is followed by an EndEvent without an IterationEndEvent
      // or LineSearchIterationStartEvent+LineSearchIterationEndEvent before,
      // one of the stopping criteria must have been fullfilled
      _NumberOfIterations = 0;
      _NumberOfSteps      = 0;
    } break;

    // Start of line search given current gradient direction
    case LineSearchStartEvent: {
      if (_Verbosity > 0) {

        model->Value(); // update cached value of individual terms if necessary

        if (debug_time) os << "\n              ";
        os << "Line search with ";
        if (step->_Info) os << step->_Info << " ";
        os << "step length in [";
        os.unsetf(ios::floatfield);
        os << setprecision(2) << step->_MinLength << ", " << step->_MaxLength << "]";

        string name;
        int i, nterms = 0;
        bool divide_by_init = false;
        for (i = 0; i < model->NumberOfTerms(); ++i) {
          const EnergyTerm *term = model->Term(i);
          if (term->Weight() != .0) {
            ++nterms;
            if (term->DivideByInitialValue()) divide_by_init = true;
          }
        }

        os << "\n\n  Energy = ";
        if (_Color) os << xboldblue;
        PrintNumber(os, step->_Current);
        if (_Color) os << xreset;
        for (i = 0; i < model->NumberOfTerms(); ++i) {
          const EnergyTerm *term = model->Term(i);
          if (term->Weight() != .0) {
            if (nterms > 1) {
              os << " = ";
              if (term->DivideByInitialValue()) {
                PrintNormalizedNumber(os, model->Value(i));
                os << " (";
                PrintNumber(os, model->RawValue(i));
                os << ")";
              } else {
                if (divide_by_init) os << "          ";
                PrintNumber(os, model->Value(i));
                if (divide_by_init) os << " ";
              }
            }
            os << "    ";
            PrintWeight(os, term->Weight(), nterms);
            name                   = term->Name();
            if (name.empty()) name = term->NameOfClass();
            os << "% " << name << "\n";
            break;
          }
        }
        for (i = i + 1; i < model->NumberOfTerms(); ++i) {
          const EnergyTerm *term = model->Term(i);
          if (term->Weight() != .0) {
            os << "                            + ";
            if (term->DivideByInitialValue()) {
              PrintNormalizedNumber(os, model->Value(i));
              os << " (";
              PrintNumber(os, model->RawValue(i));
              os << ")";
            } else {
              if (divide_by_init) os << "          ";
              PrintNumber(os, model->Value(i));
              if (divide_by_init) os << " ";
            }
            os << "    ";
            PrintWeight(os, term->Weight(), nterms);
            name                   = term->Name();
            if (name.empty()) name = term->DefaultName();
            os << "% " << name << "\n";
          }
        }
        if (_Color) os << xboldblack;
        os << "\n                 Energy        Step Length        Max. Delta\n\n";
        if (_Color) os << xreset;
      }
      _NumberOfIterations = 0;
      _NumberOfSteps      = 0;
      break;
    }

    case LineSearchIterationStartEvent:
      if (_Verbosity > 0) {
        if (debug_time) os << "\nStep " << left << setw(3) << iter->Count() << "\n";
        else            os <<   "     " << setw(3) << iter->Count() << "   ";
      }
      break;

    case AcceptedStepEvent:
    case RejectedStepEvent: {
      if (_Verbosity > 0) {
        if (debug_time) os << "\n           ";
        if (_Color) os << (event == AcceptedStepEvent ? xgreen : xbrightred);
        PrintNumber(os, step->_Value ) << "  ";
        PrintNumber(os, step->_Length) << "  ";
        PrintNumber(os, step->_Delta );
        if (_Color) os << xreset;
        else os << "    " << ((event == AcceptedStepEvent) ? "Accepted" : "Rejected");
        os << "\n";
      }
      if (event == AcceptedStepEvent) ++_NumberOfSteps;
      // Increment counter of actual line search iterations performed
      // Note that the Delta convergence criterium on the minimum maximum
      // DoF change can cause the line search to stop immediately
      ++_NumberOfIterations;
      break;
    }

    case LineSearchIterationEndEvent: {
      if (_Verbosity > 0 && iter->Count() == iter->Total()) {
        if (_Color) os << xboldblack;
        os << "\n              Maximum number of iterations exceeded\n";
        if (_Color) os << xreset;
      }
      break;
    }

    // End of line search
    case LineSearchEndEvent: {
      if (_Verbosity > 0) {
        // The minimum maximum DoF change convergence criterium kicked in immediately...
        if (_NumberOfIterations == 0) {
          os << "                              ";
          if (_Color) os << xboldred;
          PrintNumber(os, step->_Delta) << "\n";
        }
        // Report if line search was successful or no improvement
        if (_Color) os << xboldblack;
        if (_NumberOfSteps > 0) {
          if (_Verbosity > 0) {
            os << "\n               Step length = ";
            PrintNumber(os, step->_TotalLength) << " / ";
            PrintNumber(os, step->_Unit)        << "\n";
          }
        } else {
          os << "\n         No further improvement within search range\n";
        }
        if (_Color) os << xreset;
        os << "\n";
      }
      break;
    }

    // End of iteration
    case IterationEndEvent: {
      if (_Verbosity <= 0) {
        if (_Color) os << xboldblue;
        PrintNumber(os, model->Value());
        if (_Color) os << xreset;
        int nterms = 0;
        for (int i = 0; i < model->NumberOfTerms(); ++i) {
          const EnergyTerm *term = model->Term(i);
          if (term->Weight() != .0) {
            if (model->NumberOfTerms() < 5) {
              if (nterms++ == 0) os << " = ";
              else               os << " + ";
              if (term->DivideByInitialValue()) {
                PrintNormalizedNumber(os, model->Value(i));
                os << " (";
                PrintNumber(os, model->RawValue(i));
                os << ")";
              } else {
                os << " ";
                PrintNumber(os, model->Value(i));
              }
            } else {
              os << "    ";
              PrintNumber(os, model->RawValue(i));
            }
          }
        }
        if (optimizer->NumberOfStoppingCriteria() > .0 || eulermethod) {
          os << "  [";
          for (int i = 0; i < optimizer->NumberOfStoppingCriteria(); ++i) {
            if (i > 0) os << ", ";
            optimizer->StoppingCriterion(i)->Print(os);
          }
          if (eulermethod) {
            if (optimizer->NumberOfStoppingCriteria() > .0) os << ", ";
            os << "delta = ";
            PrintDelta(os, eulermethod->LastDelta(), eulermethod->StepLength());
            os << "mm";
          }
          os << "]";
        }
      }
      if (_NumberOfSteps == 0) _NumberOfSteps = 1; // no line search (cf. EulerMethod)
    } break;

    // Energy function modified after convergence and optimization restarted
    case RestartEvent: {
      os << "\n";
      if (_Verbosity > 0) {
        os << "\nContinue with modified energy function\n";
      }
    } break;

    // End of optimization
    case EndEvent: {
      os << "\n";
      if (optimizer->Converged()) {
        os << "          ";
        if (_Color) os << xboldblack;
        if (eulermethod && eulermethod->LastDelta() <= eulermethod->Delta()) {
          os << "Minimum delta stopping criterion fulfilled";
        } else {
          os << "No further improvement";
          if (optimizer->NumberOfStoppingCriteria() > 0) {
            os << " or stopping criterion fulfilled";
          }
        }
        os << "\n";
        if (_Color) os << xreset;
      }
      break;
    }

    // Status message broadcasted by registration filter
    case StatusEvent: {
      if (_Color) os << xboldblack;
      os << msg;
      if (_Color) os << xreset;
      break;
    }

    // Log message broadcasted by registration filter
    case LogEvent: {
      if (_Verbosity > 0) os << msg;
      break;
    }

    default: break;
  }

  // Flush output buffer
  if (_FlushBuffer) os.flush();

  // Reset stream output format
  os.width(w);
  os.precision(p);
}


} // namespace mirtk
