//
// File: ParameterAliasable.h
// Created by: Julien Dutheil
// Created on: Thu May 14 16:53 2009
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

#ifndef _PARAMETERALIASABLE_H_
#define _PARAMETERALIASABLE_H_

#include "Parametrizable.h"
#include "ParameterExceptions.h"
#include "ParameterList.h"

//From the STL:
#include <string>

/**
 * @brief Extend the Parametrizable interface with support for parameter aliases.
 *
 * Parameter aliases allows several parameter to be constrained together, and 
 * for instance, be jointly estimated.
 */
namespace bpp
{

class ParameterAliasable :
  public virtual Parametrizable
{
  public:
    ParameterAliasable() {}
    virtual ~ParameterAliasable() {}

  public:

    /**
     * @brief Get the number of independent parameters.
     *
     * @return The number of independent parameters.
     * If no parameters are aliased, this is equivalent to the getNumberOfParameters() method.
     */
    virtual unsigned int getNumberOfIndependentParameters() const = 0;

    /**
     * @brief Set two parameters as 'aliased'.
     *
     * The values of the two parameters will be synchronized, so that setting the value of one parameter will automatically set the value of the other one accordingly.
     * @param p1 Original parameter.
     * @param p2 Aliased parameter.
     * @throw ParameterNotFoundException if p1 or p2 do not correspond to existing parameters.
     * @throw Exception when trying to perform non-valid association.
     */
    virtual void aliasParameters(const string & p1, const string & p2) throw (ParameterNotFoundException, Exception) = 0; 

    /**
     * @brief Detach two parameters previously set as 'aliased'.
     *
     * The values of the two parameters will now be independent.
     * @param p1 Original parameter.
     * @param p2 Aliased parameter.
     * @throw ParameterNotFoundException if p1 or p2 do not correspond to existing parameters.
     * @throw Exception when trying to perform non-valid dissociation.
      */
    virtual void unaliasParameters(const string & p1, const string & p2) throw (ParameterNotFoundException, Exception)  = 0;

    /**
     * @brief Get the minimal list of parameters to set the model.
     *
     * If no parameters are aliased, this is the same a getParameters().
     *
     * @return A minimal set of parameters.
     */
    virtual const ParameterList & getIndependentParameters() const = 0;

};




/**
 * @brief A low-level implementation of the ParameterAliasable interface with void functions.
 *
 * @see Parameter, ParameterList, ParameterAliasable
 */
class ParameterAliasableAdapter:
  public ParametrizableAdapter
{
	public:
		ParameterAliasableAdapter() {}
		virtual ~ParameterAliasableAdapter() {}

	public:

		/**
		 * @name The ParameterAliasable interface.
		 *
		 * @{
		 */
		const ParameterList & getIndependentParameters() const { return getParameters(); }
    void aliasParameters(const string & p1, const string & p2) throw (ParameterNotFoundException, Exception) {}
    void unaliasParameters(const string & p1, const string & p2) throw (ParameterNotFoundException, Exception) {}
    unsigned int getNumberOfIndependentParameters() const{ return 0; }
		/** @} */

};

} // end of namespace bpp.

#endif // _PARAMETERALIASABLE_H_

