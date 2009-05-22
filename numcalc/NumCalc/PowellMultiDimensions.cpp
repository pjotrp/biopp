//
// File: PowellMultiDimensions.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov 17 15:16:45 2003
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

#include "PowellMultiDimensions.h"

#include "NumTools.h"
#include "BrentOneDimension.h"
#include "OneDimensionOptimizationTools.h"

using namespace bpp;

/******************************************************************************/

bool PowellMultiDimensions::PMDStopCondition::isToleranceReached() const
{
  // NRC Test for done:
  const PowellMultiDimensions * pmd = dynamic_cast<const PowellMultiDimensions *>(_optimizer);
  double fp   = pmd->_fp;
  double fret = pmd->_fret;
  return 2.0 * NumTools::abs(fp - fret) <= _tolerance * (NumTools::abs(fp) + NumTools::abs(fret));
}
    
/******************************************************************************/

PowellMultiDimensions::PowellMultiDimensions(Function * function) :
AbstractOptimizer(function), _f1dim(function)
{
  _defaultStopCondition = new PMDStopCondition(this);
  _stopCondition = dynamic_cast<OptimizationStopCondition *>(_defaultStopCondition->clone());
}

/******************************************************************************/

void PowellMultiDimensions::doInit(const ParameterList & params) throw (Exception)
{
  // Build the initial matrix:
  unsigned int n = params.size();
  _xi.resize(n);
  for(unsigned int i = 0; i < n; i++)
  {
    // Copy the parameter list:
    _xi[i].resize(n);
    for(unsigned int j = 0; j < n; j++)
    {
      // Set the directions to unit vectors:
      _xi[i][j] = (j == i) ? 1 : 0;
    }
  }
  
  // Starting point:
  _fret = _function->f(_parameters);
  _pt   = _parameters;
}
  
/******************************************************************************/
  
double PowellMultiDimensions::doStep() throw (Exception)
{
  unsigned int n = _parameters.size();
  _fp = _fret;
  unsigned int ibig = 0;
  double del = 0.0; // Will be the biggest function decrease
  Vdouble xit(n);
  
  // In each iteration, loop over all directions in the set.
  double fptt;
  for(unsigned int i = 0; i < n; i++)
  {
    // Copy the direction:
    for(unsigned int j = 0; j < n; j++)
    {
      xit[j] = _xi[j][i];
    }
    fptt = _fret;
    _nbEval += OneDimensionOptimizationTools::lineMinimization(_f1dim, _parameters, xit, _stopCondition->getTolerance(), NULL, _messageHandler, _verbose > 0 ? _verbose - 1 : 0);
    _fret = _function->f(_parameters);
    if(_verbose > 2) printPoint(_parameters, _fret);
    if(_fret > _fp) throw Exception("DEBUG: PowellMultiDimensions::doStep(). Line minimization failed!");
    if (fptt - _fret > del)
    {
      del = fptt - _fret;
      ibig = i;
    }
  }

  ParameterList ptt = _parameters;
  for(unsigned int j = 0; j < n; j++)
  {
    ptt[j]->setValue(2.0 * _parameters[j]->getValue() - _pt[j]->getValue());
    xit[j] = _parameters[j]->getValue() - _pt[j]->getValue();
    _pt[j]->setValue(_parameters[j]->getValue());
  }
  fptt = _function->f(ptt);
  if (fptt < _fp)
  {
    double t = 2.0 * (_fp - 2.0 * _fret + fptt) * NumTools::sqr(_fp - _fret - del) - del * NumTools::sqr(_fp - fptt);
    if (t < 0.0)
    {
      //cout << endl << "New direction: drection " << ibig << " removed." << endl;
      _nbEval += OneDimensionOptimizationTools::lineMinimization(_f1dim, _parameters, xit, _stopCondition->getTolerance(), NULL, _messageHandler, _verbose > 0 ? _verbose - 1 : 0);
      _fret = _function->f(_parameters);
      if(_fret > _fp) throw Exception("DEBUG: PowellMultiDimensions::doStep(). Line minimization failed!");
      for(unsigned int j = 0; j < n; j++)
      {
        _xi[j][ibig]  = _xi[j][n - 1];
        _xi[j][n - 1] = xit[j];
      }
    }
  }
  else _function->setParameters(_parameters);

  return _fret;
}

/******************************************************************************/

double PowellMultiDimensions::optimize() throw (Exception)
{
  AbstractOptimizer::optimize();
  // Apply best parameter:
  return _function->f(_parameters);
}

/******************************************************************************/

