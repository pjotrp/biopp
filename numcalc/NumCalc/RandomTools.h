//
// File RandomTools.h
// Author : Julien Dutheil
//          Sylvain Gaillard
// Last modification : Thu November 6 2008
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

#ifndef _RANDOMTOOLS_H_
#define _RANDOMTOOLS_H_

// From the STL:
#include <cmath>
#include <cassert>
#include <ctime>
#include <vector>
using namespace std;

#include "RandomFactory.h"

// From Utils:
#include <Utils/Exceptions.h>

namespace bpp
{

/**
 * @brief Utilitary function dealing with random numbers.
 *
 * This class uses Uniform01K generator by default.
 * It is possible to change this by setting the DEFAULT_GENERATOR variable.
 *
 * This class is adapted from Pupko's SEMPHY library.
 * It also borrow some code from Yang's PAML package.
 *
 * @see RandomFactory
 */
class RandomTools
{
  public:
    RandomTools() {}
    virtual ~RandomTools() {}

  public:
    static RandomFactory * DEFAULT_GENERATOR;
    
    /**
     * @brief Get a double random value (between 0 and specified range).
     *
     * Note : the number you get is between 0 and entry not including entry !
     * @param entry Max number to reach.
     * @param generator Random number generator to use.
     */
    static double giveRandomNumberBetweenZeroAndEntry(double entry, const RandomFactory * generator = DEFAULT_GENERATOR);

    /**
     * @brief Get a boolean random value.
     *
     * @param generator Random number generator to use.
     */
    static bool flipCoin(const RandomFactory * generator = DEFAULT_GENERATOR);

    /**
     * @brief Get an integer random value (between 0 and specified range).
     *
     * Note : the number you get is between 0 and entry not including entry !
     * @param entry Max number to reach.
     * @param generator Random number generator to use.
     */
    static int giveIntRandomNumberBetweenZeroAndEntry(int entry, const RandomFactory * generator = DEFAULT_GENERATOR);

    /**
     * @brief Set the default generator seed.
     *
     * @param seed New seed.
     */
    static void setSeed(long seed);

    /**
     * @return A random number drawn from a normal distribution.
     * @param mean The mean of the law.
     * @param variance The variance of the law.
     * @param generator The uniform generator to use.
     */
    static double randGaussian(double mean, double variance, const RandomFactory * generator = DEFAULT_GENERATOR);
    
    /**
     * @return A random number drawn from a gamma distribution with unit scale (beta=1).
     * @param dblAlpha The alpha parameter.
     * @param generator The uniform generator to use.
     */
    static double randGamma(double dblAlpha, const RandomFactory * generator = DEFAULT_GENERATOR);

    /**
     * @return A random number drawn from a gamma distribution.
     * @param alpha The alpha parameter.
     * @param beta The beta parameter.
     * @param generator The uniform generator to use.
     */
    static double randGamma(double alpha, double beta, const RandomFactory * generator = DEFAULT_GENERATOR);
  
    /**
     * @return A random number drawn from an exponential distribution.
     * @param mean The mean of the distribution.
     * @param generator The uniform generator to use.
     */
    static double randExponential(double mean, const RandomFactory * generator = DEFAULT_GENERATOR);

    /**
     * @brief Sample a vector.
     *
     * The sample is a new vector of the specified size.
     * If the size of the sample is identical to the original vector,
     * the result is a shuffle of the original vector.
     *
     * @param v The vector to sample.
     * @param size The size of the sample.
     * @param replace Should sampling be with replacement?
     * @return A vector which is a sample of v.
     * @throw IndexOutOfBoundException if the sample size exceeds the original size
     * when sampling without replacement.
     */
    template<class T> 
    static vector<T> getSample(const vector<T> & v, unsigned int size, bool replace = false) throw (IndexOutOfBoundsException) {
      if (size > v.size() && !replace)
        throw IndexOutOfBoundsException("RandomTools::getSample: size exceeded v.size.", size, 0, v.size());
      vector<unsigned int> hat;
      for (unsigned int i = 0 ; i < v.size() ; i++)
        hat.push_back(i);
      vector<T> sample;
      for (unsigned int i = 0 ; i < size ; i++) {
        unsigned int pos = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(hat.size());
        sample.push_back(v[hat[pos]]);
        if (!replace) {
          hat[pos] = hat[hat.size() - 1];
          hat.pop_back();
        }
      }
      return sample;
    }

    /**
     * @brief Get a random state from a set of probabilities/scores.
     *
     * The input probabilities are scaled so that they sum to one.
     * If 'x' probabilities are provided as input, the output vector will contain values between 0 and 'x-1'.
     *
     * @param n The sample size.
     * @param probs The set of intput probabilities.
     * @return A vector of int values corresponding to the output states. States are supposed to be in the same order as the input probabilities, the first state being '0'.
     */ 
    static vector<unsigned int> randMultinomial(unsigned int n, const vector<double>& probs);

    /**
     * @name Probability functions.
     *
     * Adapted from Yang's PAML package.
     *
     * @{
     */
    /**
     * @brief Normal quantile function.
     *
     * Returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
     * returns (-9999) if in error
     * Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
     * Applied Statistics 22: 96-97 (AS70)
     *
     * Newer methods:
     *  Wichura MJ (1988) Algorithm AS 241: the percentage points of the
     *    normal distribution.  37: 477-484.
     *  Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
     *    points of the normal distribution.  26: 118-121.
     *
     * @param prob The probability.
     * @return The quantile corresponding to prob.
     */
    static double qNorm(double prob);
    
    /**
     * @brief Computes \f$ln\left(\Gamma\left(\alpha\right)\right)\f$ given \f$\alpha\f$.
     * 
     * Returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
     * Stirling's formula is used for the central polynomial part of the procedure.
     * Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
     * Communications of the Association for Computing Machinery, 9:684
     *
     * @param alpha Alpha parameter.
     * @return \f$ln\left(\Gamma\left(\alpha\right)\right)\f$
     */
    static double lnGamma (double alpha);

    /**
     * @brief Returns the incomplete gamma ratio I(x,alpha).
     *
     * X is the upper limit of the integration and alpha is the shape parameter.
     * returns (-1) if in error
     * ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
     * (1) series expansion     if (alpha>x || x<=1)
     * (2) continued fraction   otherwise
     * RATNEST FORTRAN by
     * Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
     * 19: 285-287 (AS32)
     *
     * @param x the upper limit of the integration.
     * @param alpha the shape parameter.
     * @param ln_gamma_alpha ln(Gamma(alpha)).
     */
    static double incompleteGamma(double x, double alpha, double ln_gamma_alpha);

    /**
     * @brief \f$\chi^2\f$ quantile function.
     * 
     * returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
     * returns -1 if in error.   0.000002<prob<0.999998
     * RATNEST FORTRAN by
     * Best DJ & Roberts DE (1975) The percentage points of the 
     * Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
     * Converted into C by Ziheng Yang, Oct. 1993.
     *
     * @param prob The probability.
     * @param v number of degree of freedom.
     * @return The quantile corresponding to prob.
     */
    static double qChisq(double prob, double v);

    /**
     * @brief \f$\chi^2\f$ cumulative probability function.
     *
     * @param x The quantile for which the probability should be computed.
     * @param v number of degree of freedom.
     * @return The corresponding probability of the quantile.
     */
    static double pChisq(double x, double v)
    {
      if(x < 0) return 0;
      return pGamma(x, v / 2, 0.5);
    }

    /**
     * @brief The Gamma quantile function.
     *
     * @param prob The probability.
     * @param alpha Alpha parameter.
     * @param beta  Beta parameter.
     * @return The quantile corresponding to prob.
     */
    static double qGamma(double prob, double alpha, double beta)
    {
      return qChisq(prob,2.0*(alpha))/(2.0*(beta));
    }

    /**
     * @brief \f$\Gamma\f$ cumulative probability function.
     *
     * @param x The quantile for which the probability should be computed.
     * @param alpha Alpha parameter.
     * @param beta  Beta parameter.
     * @return The corresponding probability of the quantile.
     */
    static double pGamma(double x, double alpha, double beta)
    {
      return incompleteGamma(beta*x, alpha, lnGamma(alpha));
    }

    /** @} */
    

  private:
    static double DblGammaGreaterThanOne(double dblAlpha, const RandomFactory * generator);
    static double DblGammaLessThanOne(double dblAlpha, const RandomFactory * generator);
};

} //end of namespace bpp.

#endif  //_RANDOMTOOLS_H_

