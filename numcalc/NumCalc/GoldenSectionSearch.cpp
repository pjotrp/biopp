//
// File: GoldenSectionSearch.cpp
// Created by: Julien Dutheil 
// Created on: Mon Nov 10 10:42:17 2003
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "GoldenSectionSearch.h"
#include "NumTools.h"
#include "NumConstants.h"
#include "OneDimensionOptimizationTools.h"

// From Utils:
#include <Utils/TextTools.h>

using namespace bpp;

/******************************************************************************/

bool GoldenSectionSearch::GSSStopCondition::isToleranceReached() const
{
  // NRC Test for done:
  const GoldenSectionSearch * gss = dynamic_cast<const GoldenSectionSearch *>(_optimizer);
  double x0 = gss -> x0;
  double x1 = gss -> x1;
  double x2 = gss -> x2;
  double x3 = gss -> x3;
  return NumTools::abs(x3 - x0) <= _tolerance * (NumTools::abs(x1) + NumTools::abs(x2));
}
    
/******************************************************************************/

GoldenSectionSearch::GoldenSectionSearch(Function * function) :
AbstractOptimizer(function)
{
  _nbEvalMax = 10000;
  _defaultStopCondition = new GSSStopCondition(this);
  _stopCondition = dynamic_cast<OptimizationStopCondition *>(_defaultStopCondition->clone());
  _isInitialIntervalSet = false;
}

/******************************************************************************/

void GoldenSectionSearch::doInit(const ParameterList & params) throw (Exception)
{
  // Set the initial value (no use here! Use setInitialValues() instead).
  if(params.size() != 1) throw Exception("GoldenSectionSearch::init(). This optimizer only deals with one parameter.");

  // Bracket the minimum.
  Bracket bracket = OneDimensionOptimizationTools::bracketMinimum(_xinf, _xsup, _function, _parameters);
  if(_verbose > 0)
  {
    printMessage("Initial bracketing:");
    printMessage("A: x = " + TextTools::toString(bracket.a.x) + ", f = " + TextTools::toString(bracket.a.f));
    printMessage("B: x = " + TextTools::toString(bracket.b.x) + ", f = " + TextTools::toString(bracket.b.f));
    printMessage("C: x = " + TextTools::toString(bracket.c.x) + ", f = " + TextTools::toString(bracket.c.f));
  }
  
  // At any given time we will keep track of four points, x0, x1, x2 and x3.
  x0 = bracket.a.x;
  x3 = bracket.c.x;
  if (NumTools::abs(bracket.c.x - bracket.b.x)
      > NumTools::abs(bracket.b.x - bracket.a.x))
  {
    // Make x0 to x1 the smaller segment,
    x1 = bracket.b.x;
    // and fill in the new point to be tried.
    x2 = bracket.b.x + NumConstants::GOLDEN_RATIO_C * (bracket.c.x - bracket.b.x);
  }
  else
  {
    x2 = bracket.b.x;
    x1 = bracket.b.x - NumConstants::GOLDEN_RATIO_C * (bracket.b.x - bracket.a.x);
  }
  // The initial function evaluations.
  // Note that we never need to evaluate the function at the original endpoints.
  _parameters[0]->setValue(x1); f1 = _function->f(_parameters);
  _parameters[0]->setValue(x2); f2 = _function->f(_parameters);
}

/******************************************************************************/

void GoldenSectionSearch::setInitialInterval(double inf, double sup)
{
  if(sup > inf)
  {
    _xinf = inf; _xsup = sup;
  }
  else
  {
    _xinf = sup; _xsup = inf;
  }
  _isInitialIntervalSet = true;
}

/******************************************************************************/

double GoldenSectionSearch::doStep() throw (Exception)
{
  if(!_isInitialIntervalSet) throw Exception("GoldenSectionSearch::step. Initial interval not set: call the 'setInitialInterval' method first!");
  
  _nbEval++;

  if (f2 < f1)
  {
    // One possible outcome, its housekeeping,
    NumTools::shift<double>(x0, x1, x2);
    x2 = NumConstants::GOLDEN_RATIO_R * x1 + NumConstants::GOLDEN_RATIO_C * x3;
    // and a new function evaluation.
    _parameters[0]->setValue(x2);
    _tolIsReached = _nbEval > 2 && _stopCondition->isToleranceReached();
    NumTools::shift<double>(f1, f2, _function->f(_parameters));
    return f2;
  }
  else
  {
    // The other outcome,
    NumTools::shift<double>(x3, x2, x1);
    x1 = NumConstants::GOLDEN_RATIO_R * x2 + NumConstants::GOLDEN_RATIO_C * x0;
    // and its new function evaluation.
    _parameters[0]->setValue(x1);
    _tolIsReached = _nbEval > 2 && _stopCondition->isToleranceReached();
    NumTools::shift<double>(f2, f1, _function->f(_parameters));
    return f1;
  }
}

/******************************************************************************/

double GoldenSectionSearch::getFunctionValue() const throw (NullPointerException)
{
  if(_function == NULL) throw NullPointerException("GoldenSectionSearch::getFunctionValue. No function associated to this optimizer.");
  //return NumTools::min(f1, f2); 
  return _currentValue; 
}

/******************************************************************************/

