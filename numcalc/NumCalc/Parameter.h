//
// File: Parameter.h
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

#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include "ParameterExceptions.h"
#include "Constraints.h"

// From the STL:
#include <string>
#include <iostream>
using namespace std;

//From Utils:
#include <Utils/Clonable.h>

namespace bpp
{

class Parameter;

class ParameterEvent:
  public Clonable
{
  protected:
    Parameter * _parameter;

  public:
    ParameterEvent(Parameter * parameter);

#ifndef NO_VIRTUALCOV
    Clonable *
#else
    ParameterEvent *
#endif
    clone() const { return new ParameterEvent(*this); }

  public:
    const Parameter * getParameter() const { return _parameter; }
    Parameter * getParameter() { return _parameter; }
};

/**
 * @The parameter listener interface.
 *
 * Imlementing this interface allows to catch events associated to parameters modifications.
 * Listeners must have an identifier that will be used to pinpoint it when attached to a list.
 * This identifier needs not be unique though, but listeners with identical id will be undistinguishable.
 */
class ParameterListener:
  public virtual Clonable
{
  public:
#ifndef NO_VIRTUALCOV
    Clonable *
#else
    ParameterListener *
#endif
    clone() const = 0;

  public:

    /**
     * @retun The identifier of this listener.
     */
    virtual const string & getId() const = 0;

    /**
     * @brief Notify a renaming action.
     *
     * @param event Event associated to the acion.
     */
    virtual void parameterNameChanged(ParameterEvent & event) = 0;
    
    /**
     * @brief Notify a value change.
     *
     * @param event Event associated to the acion.
     */
    virtual void parameterValueChanged(ParameterEvent & event) = 0;
};

/**
 * @brief This class is designed to facilitate the manipulation of parameters.
 *
 * A parameter object contains a <i>value</i> stored as a double.
 * It also contains a <i>name</i> and optionaly a constraint.
 * Constraint objects allows to apply restriction on the value of the parameter,
 * for instance positive number, or a particular interval and so on.
 *
 * @see ParameterList, Parametrizable, Constraint.
 */
class Parameter:
  public virtual Clonable
{
  protected:
    string _name;             //Parameter name
    double _value;            //Parameter value
    Constraint * _constraint; //A constraint on the value
    bool _attach;
    vector<ParameterListener *> _listeners;
    vector<bool> _listenerAttach;
  
  public: // Class constructors and destructors:
    
    /**
     * @brief Build a new parameter.
     *
     * @param name       The parameter name.
     * @param value      The parameter value.
     * @param constraint An optional pointer toward a constraint Object.
     * @param attachConstraint Tell if the constraint must be attached to this parameter, or shared
     * between different objects (the default behavior, for backward compatibility).
     * If the first case, the constraint object will be destroyed when the parameter is destroyed,
     * and duplicated whe the parameter is copied.
     * @throw ConstraintException If the parameter value does not match the contraint.
     */
    Parameter(const string & name = "", double value = 0, Constraint * constraint = NULL, bool attachConstraint = false)
    throw (ConstraintException);

    /**
     * @brief Copy constructor.
     */
    Parameter(const Parameter & param);
    
    /**
     * @brief Assignment operator.
     */
    Parameter & operator=(const Parameter & param);
  
    virtual ~Parameter();
    
#ifndef NO_VIRTUALCOV
    Clonable *
#else
    Parameter *
#endif
    clone() const { return new Parameter(*this); }
    
  public:

    /**
     * @brief Set the name of this parameter.
     *
     * @param name the new parameter name.
     */
    virtual void setName(const string & name)
    {
      _name = name;
      ParameterEvent event(this);
      fireParameterNameChanged(event);
    }
  
    /**
     * @brief Set the value of this parameter.
     *
     * @param value the new parameter value.
     */
    virtual void setValue(double value) throw (ConstraintException);
  
    /**
     * @brief Get the name of this parameter.
     *
     * @return The parameter name.
     */
    virtual string getName() const { return _name; }
  
    /**
     * @brief Get the value of this parameter.
     *
     * @return The parameter value.
     */
    virtual double getValue() const { return _value; }
    
    /**
     * @brief Return the constraint associated to this parameter if there is one.
     *
     * @return A pointer toward the constraint, or NULL if there is no constraint.
     */
    virtual const Constraint * getConstraint() const { return _constraint; }
    
    /**
     * @brief Return the constraint associated to this parameter if there is one.
     *
     * @return A pointer toward the constraint, or NULL if there is no constraint.
     */
    virtual Constraint * getConstraint() { return _constraint; }

    /**
     * @brief Tells if this parameter has a constraint.
     *
     * @return True if this parameter has a contraint.
     */
    virtual bool hasConstraint() const { return _constraint != NULL; }
    
    /**
     * @brief Remove the constraint associated to this parameter.
     *
     * Warning! The contraint objet is not deleted.
     *
     * @return A pointer toward the formerly used contraint.
     */
    virtual const Constraint * removeConstraint();

    virtual void setConstraint(Constraint * constraint) throw (ConstraintException)
    {
      if(_constraint && _constraint->isCorrect(_value))
        _constraint = constraint;
      else throw ConstraintException("Parameter::setConstraint", this, _value);
    }

    /**
     * @brief Add a new listener to this parameter.
     *
     * @param listener The listener to add.
     * @param attachListener Tell if the parameter will own this listener.
     * If so, deep copies will be made when cloning the parameter, and the listener will be destroyed upon
     * destruction of the parameter or upon removal. Alternatively, only superficial copies will be made,
     * and the listener will persist if the parameter is destroyed.
     */
    virtual void addParameterListener(ParameterListener * listener, bool attachListener = true)
    {
      _listeners.push_back(listener);
      _listenerAttach.push_back(attachListener);
    }

    /**
     * @brief Remove all listeners with a given id from this parameter.
     *
     * @param listenerId The id of listener to remove.
     */
    virtual void removeParameterListener(const string & listenerId);

  protected:
    void fireParameterNameChanged(ParameterEvent & event)
    {
      for(vector<ParameterListener *>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
        (*it)->parameterNameChanged(event);
    }
    void fireParameterValueChanged(ParameterEvent & event)
    {
      for(vector<ParameterListener *>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
        (*it)->parameterValueChanged(event);
    }
  
  public:
    static IncludingPositiveReal R_PLUS;
    static ExcludingPositiveReal R_PLUS_STAR;
    static IncludingNegativeReal R_MINUS;
    static ExcludingNegativeReal R_MINUS_STAR;
    static IncludingInterval PROP_CONSTRAINT_IN;
    static ExcludingInterval PROP_CONSTRAINT_EX;
};

} //end of namespace bpp.

#endif  //_PARAMETER_H_

