//
// File: MetaOptimizer.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 12 16:05 2007
// From file: NewtonBrentMetaOptimizer.cpp
// Created on: Tue Nov 17 17:22 2004
// 
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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

/**************************************************************************/

#include "MetaOptimizer.h"

// From Utils:
#include <Utils/ApplicationTools.h>

using namespace bpp;

/**************************************************************************/

string MetaOptimizerInfos::IT_TYPE_STEP = "step";
string MetaOptimizerInfos::IT_TYPE_FULL = "full";

/**************************************************************************/

MetaOptimizer::MetaOptimizer(
    Function * function,
    MetaOptimizerInfos * desc,
    unsigned int n):
  AbstractOptimizer(function),
  _optDesc(desc), _n(n)
{
  _defaultStopCondition = new FunctionStopCondition(this);
  _stopCondition = dynamic_cast<OptimizationStopCondition *>(_defaultStopCondition->clone());
  _precisionStep = log10(_stopCondition->getTolerance()) / _n;
  _stepCount = 0;
  _stepChar = "";
  _optParameters.resize(desc->getNumberOfOptimizers());
  _nbParameters.resize(desc->getNumberOfOptimizers());
}

/**************************************************************************/

MetaOptimizer::MetaOptimizer(
    const MetaOptimizer & opt):
  AbstractOptimizer(opt)
{
  _optDesc       = dynamic_cast<MetaOptimizerInfos *>(opt._optDesc->clone());
  _optParameters = opt._optParameters;
  _nbParameters  = opt._nbParameters;
  _n             = opt._n;
  _precisionStep = opt._precisionStep;
  _stepCount     = opt._stepCount;
}

/**************************************************************************/

MetaOptimizer & MetaOptimizer::operator=(
    const MetaOptimizer & opt)
{
  AbstractOptimizer::operator=(opt);
  _optDesc       = dynamic_cast<MetaOptimizerInfos *>(opt._optDesc->clone());
  _optParameters = opt._optParameters;
  _nbParameters  = opt._nbParameters;
  _n             = opt._n;
  _precisionStep = opt._precisionStep;
  _stepCount     = opt._stepCount;
  return *this;
}

/**************************************************************************/

MetaOptimizer::~MetaOptimizer()
{
  // Delete all optimizers:
  delete _optDesc;
}

/**************************************************************************/

void MetaOptimizer::doInit(const ParameterList & parameters)
  throw (Exception)
{
  _optParameters.resize(_optDesc->getNumberOfOptimizers());
  for (unsigned int i = 0; i < _optDesc->getNumberOfOptimizers(); i++)
  {
    _optParameters[i].reset();
    for (unsigned int j = 0; j < _optDesc->getParameterNames(i).size(); j++)
    {
      string pname = _optDesc->getParameterNames(i)[j];
      if (parameters.hasParameter(pname))
      {
        _optParameters[i].addParameter(parameters.getParameter(pname));
      }
    }
    _nbParameters[i] = _optParameters[i].size();
  }

  // Initialize optimizers:
  for(unsigned int i = 0; i < _optDesc->getNumberOfOptimizers(); i++)
  {
    if(_nbParameters[i] > 0)
    {
      Optimizer * opt = _optDesc->getOptimizer(i);
      dynamic_cast<AbstractOptimizer *>(opt)->updateParameters(_updateParameters);
      opt->setProfiler(_profiler);
      opt->setMessageHandler(_messageHandler);
      opt->setConstraintPolicy(_constraintPolicy);
      opt->setVerbose(_verbose > 0 ? _verbose - 1 : 0);
    }
  }
  
  // Actualize parameters:
  _parameters.matchParametersValues(_function->getParameters());
  
  _function->setParameters(_parameters);
  _initialValue = _function->getValue();
  // Reset counter:
  _stepCount = 1;
  // Recompute step if precision has changed:
  _precisionStep = (log10(_stopCondition->getTolerance()) - log10(_initialValue)) / _n;
}

/**************************************************************************/

double MetaOptimizer::doStep() throw (Exception)
{
  _stepCount++;
  
  int tolTest = 0;
  double tol = _stopCondition->getTolerance();
  if(_stepCount <= _n)
  {
    tol = _initialValue * pow(10, _stepCount * _precisionStep);
  }
  
  for(unsigned int i = 0; i < _optDesc->getNumberOfOptimizers(); i++)
  {
    if(_nbParameters[i] > 0)
    {
      if(_verbose > 1 && ApplicationTools::message)
      {
        *ApplicationTools::message << endl << _optDesc->getName(i) << endl;
        ApplicationTools::message->flush();
      }
      if(_optDesc->requiresFirstOrderDerivatives(i))
        dynamic_cast<DerivableFirstOrder *>(_function)->enableFirstOrderDerivatives(true);
      if(_optDesc->requiresSecondOrderDerivatives(i))  
        dynamic_cast<DerivableSecondOrder *>(_function)->enableSecondOrderDerivatives(true);

      _optParameters[i].matchParametersValues(_parameters);
      Optimizer * opt = _optDesc->getOptimizer(i);
      opt->getStopCondition()->setTolerance(tol);
      opt->init(_optParameters[i]);
      if(_optDesc->getIterationType(i) == MetaOptimizerInfos::IT_TYPE_STEP)
        opt->step();
      else if(_optDesc->getIterationType(i) == MetaOptimizerInfos::IT_TYPE_FULL)
        opt->optimize();
      else throw Exception("MetaOptimizer::step. Unknown iteration type specified.");
      _nbEval += opt->getNumberOfEvaluations();
      if(_optDesc->requiresFirstOrderDerivatives(i))
        dynamic_cast<DerivableFirstOrder *>(_function)->enableFirstOrderDerivatives(false);
      if(_optDesc->requiresSecondOrderDerivatives(i))  
        dynamic_cast<DerivableSecondOrder *>(_function)->enableSecondOrderDerivatives(false);
      if(_verbose > 1) cout << endl;
      _parameters.matchParametersValues(opt->getParameters());
    }
    tolTest += _nbParameters[i] > 0 ? 1 : 0;
  }
  _tolIsReached = (tolTest == 1);
   
  return _function->getValue();
}

/**************************************************************************/

