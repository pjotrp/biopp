//
// File: Constants.h
// Created by: Julien Dutheil
// Created on: Tue Feb 03 14:21 2009
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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

#ifndef _NUMCONSTANTS_H_
#define _NUMCONSTANTS_H_

namespace bpp {

/**
 * @brief this static class contains several useful constant values.
 */
class NumConstants
{
  public:
    /**
     * @name Golden ratio.
     *
     * The golden ratio, @f$\phi@f$ is equal to @f$\frac{1+\sqrt{5}}{2} = 1.6180339887498948482\ldots@f$.
     * We also define @f$R=\phi-1@f$ and @f$C = 1 - R@f$.
     * @{
     */
    static const double GOLDEN_RATIO_PHI;
    static const double GOLDEN_RATIO_R;
    static const double GOLDEN_RATIO_C;
    /** @} */

    static const double TINY;
    static const double VERY_TINY;
    static const double VERY_BIG;

    static const double PI;
};

}//end of namespace bpp.

#endif	//_NUMCONSTANTS_H_

