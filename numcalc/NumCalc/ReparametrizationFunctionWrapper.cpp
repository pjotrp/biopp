//
// File: ReparametrizationFunctionWrapper.h
// Created by: Julien Dutheil
// Created on: Fri Jan  30 09:30 2009
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

#include "ReparametrizationFunctionWrapper.h"

// From Utils:
#include <Utils/ApplicationTools.h>

using namespace bpp;

ReparametrizationFunctionWrapper::ReparametrizationFunctionWrapper(Function * function, bool verbose) :
  AbstractParametrizable(""),
  _function(function)
{
  _functionParameters = function->getParameters();
  for(unsigned int i = 0; i < _functionParameters.size(); i++)
  {
    Parameter * p = _functionParameters[i];
    Constraint * constraint = p->getConstraint();
    const string name = p->getName();
    double value = p->getValue();
    if(! constraint)
    {
      if(verbose) ApplicationTools::displayMessage("Parameter " + name + " does not need to be transformed.");
      PlaceboTransformedParameter p(name, value);
      addParameter_(p);
    }
    else
    {
      try
      {
        Interval* interval = dynamic_cast<Interval *>(constraint);
        if(interval)
        {
          if(verbose) ApplicationTools::displayMessage("Parameter " + name + " was tanh transformed.");
          IntervalTransformedParameter p(name, value, interval->getLowerBound(), interval->getUpperBound());
          addParameter_(p);
        }
        else throw Exception("");
      }
      catch(exception & e)
      {
        try
        {
          ExcludingPositiveReal* epr = dynamic_cast<ExcludingPositiveReal *>(constraint);
          if(epr)
          {
            if(verbose) ApplicationTools::displayMessage("Parameter " + name + " was log transformed.");
            RTransformedParameter p(name, value, epr->getLowerBound(), true);
            addParameter_(p);
          }
          else throw Exception("");
        }
        catch(exception & e)
        {
          try
          {
            IncludingPositiveReal* ipr = dynamic_cast<IncludingPositiveReal *>(constraint);
            if(ipr)
            {
              if(verbose) ApplicationTools::displayMessage("Parameter " + name + " was log transformed.");
              RTransformedParameter p(name, value, ipr->getLowerBound(), true);
              addParameter_(p);
            }
            else throw Exception("");
          }
          catch(exception & e)
          {
            try
            {
              ExcludingNegativeReal* enr = dynamic_cast<ExcludingNegativeReal *>(constraint);
              if(enr)
              {
                if(verbose) ApplicationTools::displayMessage("Parameter " + name + " was log transformed.");
                RTransformedParameter p(name, value, enr->getUpperBound(), false);
                addParameter_(p);
              }
              else throw Exception("");
            }
            catch(exception & e)
            {
              try
              {
                IncludingNegativeReal* inr = dynamic_cast<IncludingNegativeReal *>(constraint);
                if(inr)
                {
                  if(verbose) ApplicationTools::displayMessage("Parameter " + name + " was log transformed.");
                  RTransformedParameter p(name, value, inr->getUpperBound(), false);
                  addParameter_(p);
                }
                else throw Exception("");
              }
              catch(exception & e)
              {
                if(verbose) ApplicationTools::displayWarning("No transformation found for this constraint! Parameter " + p->getName());
                PlaceboTransformedParameter p(name, value);
                addParameter_(p);
              }
            }
          }
        }
      }
    }
  }
}

void ReparametrizationFunctionWrapper::fireParameterChanged (const ParameterList &parameters)
{
  //Recompute function parameters:
  //We could update only the parameter that actually changed,
  //but that would implied a quick sort on parameter names (nlog(n))
  //whereas using a loop over the set is in o(n). It should hence be 
  //more efficient in most cases.
  for(unsigned int i = 0; i < getNumberOfParameters(); i++)
  {
    double x = dynamic_cast<TransformedParameter *>(&getParameter_(i))->getOriginalValue();
    _functionParameters[i]->setValue(x);
  }
}

