//
// File: ConjugateGradientMultiDimensions.cpp
// Created by: Julien Dutheil
// Created on: Wed Apr 11 16:51 2007
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

#include "ConjugateGradientMultiDimensions.h"
#include "OneDimensionOptimizationTools.h"

using namespace bpp;

/******************************************************************************/

ConjugateGradientMultiDimensions::ConjugateGradientMultiDimensions(DerivableFirstOrder * function):
  AbstractOptimizer(function), _f1dim(function)
{
  _defaultStopCondition = new FunctionStopCondition(this);
  _stopCondition = dynamic_cast<OptimizationStopCondition *>(_defaultStopCondition->clone());
  _optimizer = new BrentOneDimension(function);
}

/******************************************************************************/

void ConjugateGradientMultiDimensions::doInit(const ParameterList & params) throw (Exception)
{
  unsigned int nbParams = params.size();
  _g.resize(nbParams);
  _h.resize(nbParams);
  _xi.resize(nbParams);
  dynamic_cast<DerivableFirstOrder *>(_function)->enableFirstOrderDerivatives(true);
  _function->setParameters(params);
  getGradient(_xi);
  for(unsigned int i = 0; i < nbParams; i++)
  {
    _g[i]  = -_xi[i];
    _xi[i] = _h[i] = _g[i];
  }
}

/******************************************************************************/

double ConjugateGradientMultiDimensions::doStep() throw (Exception)
{
  double gg, gam, f, dgg;
  unsigned int n = _parameters.size();
  //Loop over iterations.
  dynamic_cast<DerivableFirstOrder *>(_function)->enableFirstOrderDerivatives(false);
  _nbEval += OneDimensionOptimizationTools::lineMinimization(_f1dim, _parameters, _xi, _stopCondition->getTolerance(), NULL, NULL, _verbose > 0 ? _verbose - 1 : 0);
  dynamic_cast<DerivableFirstOrder *>(_function)->enableFirstOrderDerivatives(true);
  f = _function->f(_parameters);
  if(_tolIsReached)
  {
    return f;
  }
  getGradient(_xi);
  dgg = gg = 0.0;
  for(unsigned j = 0; j < n; j++)
  {
    gg += _g[j] * _g[j];
    /* dgg += xi[j] * xi[j]; */ //This statement for Fletcher-Reeves.
    dgg += (_xi[j] + _g[j]) * _xi[j]; //This statement for Polak-Ribiere.
  }
  if (gg == 0.0)
  { 
    //Unlikely. If gradient is exactly zero then
    return f;
  }
  gam = dgg / gg;
  for(unsigned int j = 0; j < n; j++)
  {
    _g[j] = -_xi[j];
    _xi[j] = _h[j] = _g[j] + gam * _h[j];
  }
  
  return f;
}

/******************************************************************************/

void ConjugateGradientMultiDimensions::getGradient(vector<double> & gradient) const
{
  for(unsigned int i = 0; i < gradient.size(); i++)
  {
    gradient[i] = dynamic_cast<DerivableFirstOrder *>(_function)->getFirstOrderDerivative(_parameters[i]->getName());
  }
}

/******************************************************************************/

