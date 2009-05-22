//
// File: SimpleNewtonMultiDimensions.cpp
// Created by: Julien Dutheil
// Created on: Thu Apr 26 15:29 2007
//

/*
Copyright or © or Copr. CNRS, (November 19, 2004)

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

/******************************************************************************/

#include "SimpleNewtonMultiDimensions.h"

using namespace bpp;

/******************************************************************************/

SimpleNewtonMultiDimensions::SimpleNewtonMultiDimensions(DerivableSecondOrder * function):
  AbstractOptimizer(function)
{
  _defaultStopCondition = new FunctionStopCondition(this);
  _stopCondition = dynamic_cast<OptimizationStopCondition *>(_defaultStopCondition->clone());
  _nbParams = 0;
  _optimizer = new NewtonOneDimension(function);
  _stepChar = "";
}

/******************************************************************************/

SimpleNewtonMultiDimensions::SimpleNewtonMultiDimensions(const SimpleNewtonMultiDimensions & opt):
  AbstractOptimizer(opt)
{
  _nbParams = opt._nbParams;
  if(opt._optimizer) _optimizer = dynamic_cast<NewtonOneDimension *>(opt._optimizer->clone());
  else               _optimizer = NULL;
}

/******************************************************************************/

SimpleNewtonMultiDimensions & SimpleNewtonMultiDimensions::operator=(const SimpleNewtonMultiDimensions & opt)
{
  AbstractOptimizer::operator=(opt);
  _nbParams = opt._nbParams;
  if(opt._optimizer) _optimizer = dynamic_cast<NewtonOneDimension *>(opt._optimizer->clone());
  else               _optimizer = NULL;
  return *this;
}

/******************************************************************************/

SimpleNewtonMultiDimensions::~SimpleNewtonMultiDimensions()
{
  delete _optimizer;
}

/******************************************************************************/

void SimpleNewtonMultiDimensions::setFunction(Function * function)
{
  AbstractOptimizer::setFunction(function);
  _optimizer->setFunction(function);
}

/******************************************************************************/

void SimpleNewtonMultiDimensions::doInit(const ParameterList & params) throw (Exception)
{
  _nbParams = params.size();
  if(_nbParams == 0) return;

  // Initialize optimizers:
  unsigned int nbEvalMax = _nbEvalMax / _nbParams;
  _optimizer->setMaximumNumberOfEvaluations(nbEvalMax);
  _optimizer->setProfiler(_profiler);
  _optimizer->setMessageHandler(_messageHandler);
  _optimizer->getStopCondition()->setTolerance(getStopCondition()->getTolerance());
  _optimizer->setConstraintPolicy(_constraintPolicy);
  _optimizer->setVerbose( _verbose > 0 ? _verbose - 1 : 0);
  _optimizer->setMaximumNumberOfCorrections(10);
  _function->setParameters(_parameters);
}

/******************************************************************************/

double SimpleNewtonMultiDimensions::doStep() throw (Exception)
{
  double f = _function->getValue();
  for(unsigned int i = 0; i < _nbParams; i++)
  {
    if(_verbose > 0)
    {
      cout << _parameters[i]->getName() << ":";
      cout.flush();
    }
    // Re-init optimizer according to new values:
    _optimizer->init(_parameters.subList(i));

    // Optimize through this dimension:
    f = _optimizer->optimize();
    if(_verbose > 0) cout << endl;
    _parameters.matchParametersValues(_function->getParameters());
    _nbEval += _optimizer->getNumberOfEvaluations(); 
  }
  _tolIsReached = _nbParams <= 1;
  return f;
}

/******************************************************************************/

