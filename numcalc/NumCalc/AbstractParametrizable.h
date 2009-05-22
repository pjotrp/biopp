//
// File: AbstractParametrizable.h
// Created by: Julien Dutheil
// Created on: Sun Mar 29 09:10 2009
// Created from file Parametrizable.h
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

#ifndef _ABSTRACTPARAMETRIZABLE_H_
#define _ABSTRACTPARAMETRIZABLE_H_

#include "Parametrizable.h"

//From the STL:
#include <map>
using namespace std;

namespace bpp
{

/**
 * @brief A partial implementation of the Parametrizable interface.
 *
 * Parameters are stored in a protected ParameterList object.
 *
 * The abstract fireParameterChanged() method is provided so that the derived class
 * know when a parameter has changed, and can be updated.
 * All methods call the corresponding method in ParameterList and then call the
 * fireParameterChanged() method.
 */
class AbstractParametrizable:
  public virtual Parametrizable
{
  private:
    ParameterList parameters_;
    string prefix_;

  public:
    AbstractParametrizable(const string& prefix) : prefix_(prefix) {}

    virtual ~AbstractParametrizable() {}

  public:
    bool hasParameter(const string & name) const { return parameters_.hasParameter(prefix_ + name); }

    const ParameterList & getParameters() const { return parameters_; }
    
    const Parameter & getParameter(const string & name) const throw (ParameterNotFoundException)
    {
      return parameters_.getParameter(prefix_ + name);
    }
  
    double getParameterValue(const string & name) const
      throw (ParameterNotFoundException)
    { 
      return getParameter(name).getValue();
    }

    void setAllParametersValues(const ParameterList & parameters) 
      throw (ParameterNotFoundException, ConstraintException)
    {
      parameters_.setAllParametersValues(parameters);
      fireParameterChanged(parameters);
    }

    void setParameterValue(const string & name, double value) 
      throw (ParameterNotFoundException, ConstraintException)
    {
      parameters_.setParameterValue(prefix_ + name, value);
      fireParameterChanged(parameters_.subList(prefix_ + name));
    }

    void setParametersValues(const ParameterList & parameters)
      throw (ParameterNotFoundException, ConstraintException)
    { 
      parameters_.setParametersValues(parameters);
      fireParameterChanged(parameters);
    }

    void matchParametersValues(const ParameterList & parameters)
      throw (ConstraintException)
    { 
      parameters_.matchParametersValues(parameters);
      fireParameterChanged(parameters);
    }

    unsigned int getNumberOfParameters() const { return parameters_.size(); }
     
    void setNamespace(const string& prefix);
    
    string getNamespace() const { return prefix_; }
    
    string getParameterNameWithoutNamespace(const string& name) const;

    /**
     * @brief Notify the class when one or several parameters have changed.
     *
     * @param parameters A ParameterList object with parameters that changed.
     */
    virtual void fireParameterChanged(const ParameterList & parameters) = 0;

  protected:
    void addParameter_(const Parameter& parameter)
    {
      parameters_.addParameter(parameter);
    }

    void addParameters_(const ParameterList& parameters)
    {
      parameters_.addParameters(parameters);
    }

    void deleteParameter_(unsigned int index) throw (IndexOutOfBoundsException)
    {
      if(index >= parameters_.size())
        throw IndexOutOfBoundsException("AbstractParametrizable::deleteParameter_.", index, 0, parameters_.size() - 1);
      string name = parameters_[index]->getName();
      parameters_.deleteParameter(index);
    }

    void resetParameters_()
    {
      parameters_.reset();
    }

    /**
     * @param The name of the parameter.
     * @return A reference toward the corresponding parameter.
     * @throw ParameterNotFoundException If no parameter with that name is found in the list.
     */
    Parameter& getParameter_(const string & name) throw (ParameterNotFoundException)
    {
      return parameters_.getParameter(prefix_ + name);
    }
  
    /**
     * @param The name of the parameter, including its namespace.
     * @return A reference toward the corresponding parameter.
     * @throw ParameterNotFoundException If no parameter with that name is found in the list.
     */
    Parameter& getParameterWithNamespace_(const string& name) throw (ParameterNotFoundException)
    {
      return getParameter_(name);
    }
    /**
     * @param The name of the parameter, including its namespace.
     * @return A reference toward the corresponding parameter.
     * @throw ParameterNotFoundException If no parameter with that name is found in the list.
     */
    const Parameter& getParameterWithNamespace_(const string& name) const throw (ParameterNotFoundException)
    {
      return getParameter(name);
    }

    Parameter& getParameter_(unsigned int index) throw (IndexOutOfBoundsException)
    {
      if(index >= parameters_.size())
        throw IndexOutOfBoundsException("AbstractParametrizable::getParameter_.", index, 0, parameters_.size() - 1);
      return *parameters_[index];
    }
    const Parameter& getParameter_(unsigned int index) const throw (IndexOutOfBoundsException)
    {
      if(index >= parameters_.size())
        throw IndexOutOfBoundsException("AbstractParametrizable::getParameter_.", index, 0, parameters_.size() - 1);
      return *parameters_[index];
    }
    
    ParameterList& getParameters_() { return parameters_; }
};

} //end of namespace bpp.

#endif //_ABSTRACTPARAMETRIZABLE_H_

