//
// File: AbstractParameterAliasable.h
// Created by: Julien Dutheil
// Created on: Thu May 14 17:08 2009
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

#ifndef _ABSTRACTPARAMETERALIASABLE_H_
#define _ABSTRACTPARAMETERALIASABLE_H_

#include "AbstractParametrizable.h"
#include "ParameterAliasable.h"

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
class AbstractParameterAliasable:
  public AbstractParametrizable,
  public virtual ParameterAliasable
{
  private:

    mutable ParameterList independentParameters_;

    class AliasParameterListener:
      public ParameterListener
    {
      protected:
        string id_;
        unsigned int alias_;
        ParameterList *pl_;
        string name_;

      public:
        AliasParameterListener(const string& id, unsigned int alias, ParameterList* pl): id_(id), alias_(alias), pl_(pl)
        {
          //This allow us to check if the parameter position have changed at some point...
          name_ = (*pl_)[alias]->getName();
        }
        AliasParameterListener * clone() const { return new AliasParameterListener(*this); }

      public:
        const string & getId() const { return id_; }

        void setParameterList(ParameterList *pl) { pl_ = pl; }

        void parameterNameChanged(ParameterEvent & event) throw (Exception)
        {
          Parameter * p = (*pl_)[alias_];
          if(p->getName() != name_)
            throw Exception("AbstractParametrizable::AliasParameterListener::parameterNameChanged. Error, aliased parameter have change, maybe because a parameter was removed?");
          p->setName(event.getParameter()->getName());
          name_ = p->getName();
        }
    
        void parameterValueChanged(ParameterEvent & event) throw (Exception)
        {
          Parameter * p = (*pl_)[alias_];
          if(p->getName() != name_)
            throw Exception("AbstractParametrizable::AliasParameterListener::parameterValueChanged. Error, aliased parameter have change, maybe because a parameter was removed?");
          p->setValue(event.getParameter()->getValue());
        }

        const string& getName() const { return name_; }

        void rename(const string& name) { name_ = name; }
      
    };

    /**
     * Contains all parameter listeners for maintening alias relationships.
     * The registry will be updated appropriately upon cloning and deleting.
     */
    map<string, AliasParameterListener *> aliasListenersRegister_;
  
  public:
    AbstractParameterAliasable(const string& prefix) : AbstractParametrizable(prefix) {}

    AbstractParameterAliasable(const AbstractParameterAliasable & ap);
    
    AbstractParameterAliasable & operator=(const AbstractParameterAliasable & ap);

    virtual ~AbstractParameterAliasable();

  public:
    void setNamespace(const string& prefix);
 
    const ParameterList & getIndependentParameters() const { return independentParameters_; }
    
    unsigned int getNumberOfIndependentParameters() const { return independentParameters_.size(); }

    void aliasParameters(const string & p1, const string & p2) throw (ParameterNotFoundException, Exception);

    void unaliasParameters(const string & p1, const string & p2) throw (ParameterNotFoundException, Exception);

    void fireParameterChanged(const ParameterList& parameters)
    {
      independentParameters_.matchParametersValues(getParameters());
    }

  protected:
    void addParameter_(const Parameter& parameter)
    {
      AbstractParametrizable::addParameter_(parameter);
      independentParameters_.addParameter(parameter);
    }

    void addParameters_(const ParameterList& parameters)
    {
      AbstractParametrizable::addParameters_(parameters);
      independentParameters_.addParameters(parameters);
    }

    void deleteParameter_(unsigned int index) throw (IndexOutOfBoundsException)
    {
      string name = getParameter_(index).getName();
      AbstractParametrizable::deleteParameter_(index);
      if(independentParameters_.hasParameter(name))
        independentParameters_.deleteParameter(name);
    }

    void resetParameters_()
    {
      AbstractParametrizable::resetParameters_();
      independentParameters_.reset();
    }

};

} //end of namespace bpp.

#endif //_ABSTRACTPARAMETERALIASABLE_H_

