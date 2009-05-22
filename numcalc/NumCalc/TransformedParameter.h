//
// File: TransformedParameter.h
// Created by: Julien Dutheil
// Created on: Fri Jan 30 09:42 2009
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

#ifndef _TRANSFORMEDPARAMETER_H_
#define _TRANSFORMEDPARAMETER_H_

#include "Parameter.h"
#include "NumConstants.h"

#include <cmath>
using namespace std;

namespace bpp
{

/**
 * @brief The TransformedParameter interface.
 *
 * This interface extends the Parameter interface.
 * A transformed parameter does not have a constraint attached to it, and is supposed to range from -inf to +inf.
 * It uses a transformation in order to do this, typically using a bijection.
 * The exact function used to achieve the transformation depends on the implementation of the interface.
 */ 
class TransformedParameter:
  public virtual Parameter
{
	public:
		
#ifndef NO_VIRTUAL_COV
		TransformedParameter * clone() const = 0;
#endif
	
	public:
    /**
     * @brief Set the value of the parameter using the orignal coordinate system.
     *
     * @param value Parameter value in original coordinates.
     * @throw ConstraintException if the value is not correct.
     */
    virtual void setOriginalValue(double value) throw (ConstraintException) = 0;

    /**
     * @return The current value of the parameter in orignal coordinates.
     */
    virtual double getOriginalValue() const = 0;
	
};

/**
 * @brief Parameter transformation from ] b, +inf [ or ] -inf, b [ to ]-inf, + inf [.
 *
 * The equation of the tranformation is 
 * @f[
 * x' = \begin{cases}
 *   \log(a\cdot(x-b)) & \text{if $x < b+1$},\\
 *   a(x-1-b)          & \text{if $a \geq b+1$}.
 * \end{cases}
 * @f]
 * for a transformation from ] b, +inf [ to ]-inf, + inf [.
 * The 'b' parameter is the lower bound and 'a' is a scaling factor set to 1 by default.
 * For a transformation from  ] -inf, b [, the transformation is then
 * @f[
 * x' = \begin{cases}
 *   -\log(-a\cdot(x-b)) & \text{if $x < b-1$},\\
 *   -a(x-1-b)           & \text{if $a \geq b-1$}.
 * \end{cases}
 * @f]
 */
class RTransformedParameter:
  public virtual TransformedParameter
{
  protected:
    double _scale;
    double _bound;
    bool _positive;

  public:
    /**
     * @brief Build a new RTransformedParameter, with given bound and scale.
     *
     * @param name the name of the parameter.
     * @param value the value of th eparameter, in orginal coordinates.
     * @param bound the inerval bound to use.
     * @param positive tell if the original interval is positive or negative.
     * @param scale the scaling factor.
     */
    RTransformedParameter(const string & name, double value, double bound = 0, bool positive = true, double scale = 1):
      Parameter(name, 1.),
      _scale(scale),
      _bound(bound),
      _positive(positive)
    {
      setOriginalValue(value);
    }

#ifndef NO_VIRTUAL_COV
    RTransformedParameter * clone() const { return new RTransformedParameter(*this); }
#endif

  public:
    void setOriginalValue(double value) throw (ConstraintException) 
    {
      if(_positive ? value <= _bound : value >= _bound) throw ConstraintException("RTransformedParameter::setValue", this, value);
      if(_positive  & (value < 1 + _bound)) setValue(log(_scale * (value - _bound)));
      if(_positive  & (value >= 1 + _bound)) setValue(_scale * (value - 1. - _bound));
      if(!_positive & (value > -1 + _bound)) setValue(log(-_scale * (value - _bound)));
      if(!_positive & (value <= -1 + _bound)) setValue(-_scale * (value - 1. - _bound));
    }

    double getOriginalValue() const
    {
      double x = getValue();
      if(_positive)
        if(x < 0) return exp(x)/_scale + _bound;
        else      return x/_scale + 1. + _bound;
      else 
        if(x < 0) return - exp(-x)/_scale + _bound;
        else      return -x/_scale - 1. + _bound;
    }
};

/**
 * @brief Parameter transformation from ] a, b [ to ]-inf, + inf [.
 *
 * The equation of the tranformation is 
 * @f[
 * x' = s\tan\left(\pi\frac{x-a}{b-a} - \frac{\pi}{2}\right)
 * @f]
 * The 'a' and 'b' parameters are the lower and upper bounds and 's' is a sclaing factor set to 1 by default.
 * If the hyperbolic option is set to true (the default), then the following transformation is used instead:
 * @f[
 * x' = s\,\text{atanh}\left(2\frac{x-a}{b-a} - 1\right)
 * @f]
 *
 */
class IntervalTransformedParameter:
  public virtual TransformedParameter
{
  protected:
    double _scale;
    double _lowerBound;
    double _upperBound;
    bool _hyper;
    double _tiny;

  public:
    /**
     * @brief Build a new IntervalTransformedParameter, with given bounds and scale.
     *
     * @param name the name of the parameter.
     * @param value the value of th eparameter, in orginal coordinates.
     * @param lowerBound the inerval lower bound to use.
     * @param upperBound the inerval lower bound to use.
     * @param scale the scaling factor.
     * @param hyper tell if the hyberboic function should be used (true by default).
     */
    IntervalTransformedParameter(const string & name, double value, double lowerBound = 0, double upperBound = 1, double scale = 1, bool hyper = true):
      Parameter(name, hyper ? scale * atanh(2. * (value - lowerBound) / (upperBound - lowerBound) - 1.) : scale * tan(3.141593 * (value - lowerBound)/(upperBound - lowerBound) - 1.570796)),
      _scale(scale),
      _lowerBound(lowerBound),
      _upperBound(upperBound),
      _hyper(hyper),
      _tiny(NumConstants::TINY)
    {}

#ifndef NO_VIRTUAL_COV
    IntervalTransformedParameter * clone() const { return new IntervalTransformedParameter(*this); }
#endif

  public:
    void setOriginalValue(double value) throw (ConstraintException) 
    {
      if(value <= _lowerBound || value >= _upperBound) throw ConstraintException("IntervalTransformedParameter::setValue", this, value);
      setValue(_hyper ?_scale * atanh(2. * (value - _lowerBound) / (_upperBound - _lowerBound) - 1.) : _scale * tan(NumConstants::PI * (value - _lowerBound)/(_upperBound - _lowerBound) - NumConstants::PI / 2.));
    }

    double getOriginalValue() const
    {
      double x = getValue();
      double x2 = _hyper ?  (tanh(x / _scale) + 1.) / 2. : (atan(x / _scale) + NumConstants::PI / 2.) / NumConstants::PI;
      //Prevent rounding error:
      if(x2 < _lowerBound + _tiny) x2 = _lowerBound + _tiny;
      if(x2 > _upperBound - _tiny) x2 = _upperBound - _tiny;
      return x2;
    }
};

/**
 * @brief 'Placebo' parameter transformation from ] b, +inf [ or ] -inf, b [ to ]-inf, + inf [.
 *
 * The class create a Transformed parameter which is exactly the same as a standard parameter.
 * It only implements the setOriginalValue and getOriginalValue methods, and remove the constraint.
 */
class PlaceboTransformedParameter:
  public virtual TransformedParameter
{
  public:
    PlaceboTransformedParameter(const string & name, double value):
      Parameter(name, value)
    {}

#ifndef NO_VIRTUAL_COV
    PlaceboTransformedParameter * clone() const { return new PlaceboTransformedParameter(*this); }
#endif

  public:
    void setOriginalValue(double value) throw (ConstraintException) 
    {
      setValue(value);
    }

    double getOriginalValue() const
    {
      return getValue();
    }
};


} //end of namespace bpp.

#endif	//_TRANSFORMEDPARAMETER_H_

