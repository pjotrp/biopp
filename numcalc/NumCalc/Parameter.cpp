//
// File: Parameter.cpp
// Created by: Julien Dutheil
// Created on: Wed Oct 15 15:40:47 2003
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

#include "Parameter.h"

//From Utils:
#include <Utils/TextTools.h>

using namespace bpp;

#include <iostream>
using namespace std;

/******************************************************************************/

ParameterEvent::ParameterEvent(Parameter * parameter): _parameter(parameter) {}

/** Constructors: *************************************************************/

Parameter::Parameter(const string & name, double value, Constraint * constraint, bool attachConstraint) 
throw (ConstraintException)
{
	_name = name;
	_constraint = constraint;
	// This may throw a ConstraintException:
	setValue(value);
  _attach = attachConstraint;
}

Parameter::Parameter(const Parameter & p)
{
	_name           = p._name;
	_value          = p._value;
  _attach         = p._attach;
	if(p._attach && p._constraint)
    _constraint   = p._constraint->clone();
  else
    _constraint   = p._constraint;
  _listeners      = p._listeners;
  _listenerAttach = p._listenerAttach;
  for(unsigned int i = 0; i < _listeners.size(); i++)
    if(_listenerAttach[i])
      _listeners[i] = dynamic_cast<ParameterListener *>(p._listeners[i]->clone());
}

Parameter & Parameter::operator=(const Parameter & p)
{
	_name           = p._name;
	_value          = p._value;
  _attach         = p._attach;
	if(p._attach && p._constraint)
    _constraint   = p._constraint->clone();
  else
    _constraint   = p._constraint;
  _listeners      = p._listeners;
  _listenerAttach = p._listenerAttach;
  for(unsigned int i = 0; i < _listeners.size(); i++)
    if(_listenerAttach[i])
      _listeners[i] = dynamic_cast<ParameterListener *>(p._listeners[i]->clone());
	return *this;	
}

/** Destructor: ***************************************************************/

Parameter::~Parameter()
{
  if(_attach && _constraint) delete _constraint;
  for(unsigned int i = 0; i < _listeners.size(); i++)
    if(_listenerAttach[i])
      delete _listeners[i];
} 

/** Value: ********************************************************************/

void Parameter::setValue(double value) throw (ConstraintException)
{
	if(_constraint != NULL && !_constraint->isCorrect(value)) 
		throw ConstraintException("Parameter::setValue", this, value);
	_value = value;
  ParameterEvent event(this);
  fireParameterValueChanged(event);
}

/** Constraint: ***************************************************************/

const Constraint * Parameter::removeConstraint()
{
	const Constraint * c = _constraint;
	_constraint = NULL;
	return c;
}

/******************************************************************************/

void Parameter::removeParameterListener(const string & listenerId)
{
  for(unsigned int i = 0; i < _listeners.size(); i++)
  {
    if(_listeners[i]->getId() == listenerId)
    {
      if(_listenerAttach[i]) delete _listeners[i];
      _listeners.erase(_listeners.begin() + i);
      _listenerAttach.erase(_listenerAttach.begin() + i);
    }
  }
}

/******************************************************************************/

IncludingPositiveReal Parameter::R_PLUS(0);

ExcludingPositiveReal Parameter::R_PLUS_STAR(0);

IncludingNegativeReal Parameter::R_MINUS(0);

ExcludingNegativeReal Parameter::R_MINUS_STAR(0);

IncludingInterval Parameter::PROP_CONSTRAINT_IN(0, 1);

ExcludingInterval Parameter::PROP_CONSTRAINT_EX(0, 1);

/******************************************************************************/

