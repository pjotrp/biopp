//
// File: SimpleDiscreteDistribution.h
// Created by: Julien Dutheil
// Created on: ?
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

#ifndef _SIMPLEDISCRETEDISTRIBUTION_H_
#define _SIMPLEDISCRETEDISTRIBUTION_H_

#include "AbstractDiscreteDistribution.h"

// From the STL:
#include <map>
using namespace std;

namespace bpp
{

/**
 * @brief A DiscreteDistribution object builds from a map<double,double> object.
 *
 * Keys are taken to be interval values, and map values to be the corresponding probabilities.
 */
class SimpleDiscreteDistribution:
  public AbstractDiscreteDistribution
{
  public:
    /**
     * @brief Build a new SimpleDiscreteDistribution object.
     *
     * @param distribution The map object to use.
     */
    SimpleDiscreteDistribution(const map<double, double> & distribution);

    virtual ~SimpleDiscreteDistribution() {}

#if defined(NO_VIRTUAL_COV)
    Clonable * clone() const { return new SimpleDiscreteDistribution(*this); }
#else
    SimpleDiscreteDistribution * clone() const { return new SimpleDiscreteDistribution(*this); }
#endif

  public:
    void fireParameterChanged(const ParameterList & parameters) {}
    Domain getDomain() const;
};

} //end of namespace bpp.

#endif  //_SIMPLEDISCRETEDISTRIBUTION_H_

