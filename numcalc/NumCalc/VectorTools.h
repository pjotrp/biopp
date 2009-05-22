//
// File: VectorTools.h
// Created by: Julien Dutheil
// Created on: Fri Mar 14 14:16:32 2003
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

#ifndef _VECTORTOOLS_H_
#define _VECTORTOOLS_H_

#include "VectorExceptions.h"
#include "NumTools.h"

#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

namespace bpp
{

typedef vector<double> Vdouble;
typedef vector<Vdouble> VVdouble;
typedef vector<VVdouble> VVVdouble;
typedef vector<VVVdouble> VVVVdouble;

typedef vector<int> Vint;
typedef vector<Vint> VVint;
typedef vector<VVint> VVVint;
typedef vector<VVVint> VVVVint;

/**
 * @name Element-wise operations.
 * @{
 */

template<class T>
vector<T>  operator+ (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
  unsigned int size;
  if(v1.size() != v2.size()) {
    throw DimensionException("VectorOperators::operator+", v1.size(), v2.size());
  } else {
    size = v1.size();
  }
  vector<T> result(size);
  for(unsigned int i = 0; i < size; i++) {
    result[i] = v1[i] + v2[i];
  }
  return result;
}

template<class T>
vector<T> operator- (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
  unsigned int size;
  if(v1.size() != v2.size()) {
    throw DimensionException("VectorOperators::operator-", v1.size(), v2.size());
  } else {
    size = v1.size();
  }
  vector<T> result(size);
  for(unsigned int i = 0; i < size; i++) {
    result[i] = v1[i] - v2[i];
  }
  return result;
}

template<class T>
vector<T> operator* (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
  unsigned int size;
  if(v1.size() != v2.size()) {
    throw DimensionException("VectorOperators::operator*", v1.size(), v2.size());
  } else {
    size = v1.size();
  }
  vector<T> result(size);
  for(unsigned int i = 0; i < size; i++) {
    result[i] = v1[i] * v2[i];
  }
  return result;
}

template<class T>
vector<T> operator/ (const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
  unsigned int size;
  if(v1.size() != v2.size()) {
    throw DimensionException("VectorOperators::operator/", v1.size(), v2.size());
  } else {
    size = v1.size();
  }
  vector<T> result(size);
  for(unsigned int i = 0; i < size; i++) {
    result[i] = v1[i] / v2[i];
  }
  return result;
}



template<class T, class C>
vector<T> operator+ (const vector<T> & v1, const C & c)
{
  vector<T> result(v1.size());
  for(unsigned int i = 0; i < result.size(); i++) {
    result[i] = v1[i] + c;
  }
  return result;
}
template<class T, class C>
vector<T> operator+ (const C & c, const vector<T> & v1)
{
  vector<T> result(v1.size());
  for(unsigned int i = 0; i < result.size(); i++) {
    result[i] = c + v1[i];
  }
  return result;
}

template<class T, class C>
vector<T> operator- (const vector<T> & v1, const C & c)
{
  vector<T> result(v1.size());
  for(unsigned int i = 0; i < result.size(); i++) {
    result[i] = v1[i] - c;
  }
  return result;
}
template<class T, class C>
vector<T> operator- (const C & c, const vector<T> & v1)
{
  vector<T> result(v1.size());
  for(unsigned int i = 0; i < result.size(); i++) {
    result[i] = c - v1[i];
  }
  return result;
}

template<class T, class C>
vector<T> operator* (const vector<T> & v1, const C & c)
{
  vector<T> result(v1.size());
  for(unsigned int i = 0; i < result.size(); i++) {
    result[i] = v1[i] * c;
  }
  return result;
}
template<class T, class C>
vector<T> operator* (const C & c, const vector<T> & v1)
{
  vector<T> result(v1.size());
  for(unsigned int i = 0; i < result.size(); i++) {
    result[i] = c * v1[i];
  }
  return result;
}

template<class T, class C>
vector<T> operator/ (const vector<T> & v1, const C & c)
{
  vector<T> result(v1.size());
  for(unsigned int i = 0; i < result.size(); i++) {
    result[i] = v1[i] / c;
  }
  return result;
}
template<class T, class C>
vector<T> operator/ (const C & c, const vector<T> & v1)
{
  vector<T> result(v1.size());
  for(unsigned int i = 0; i < result.size(); i++) {
    result[i] = c / v1[i];
  }
  return result;
}



template<class T>
void operator+= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
  for(unsigned int i = 0; i < v1.size(); i++) {
    v1[i] += v2[i];
  }
}

template<class T>
void operator-= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
  for(unsigned int i = 0; i < v1.size(); i++) {
    v1[i] -= v2[i];
  }
}

template<class T>
void operator*= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
  for(unsigned int i = 0; i < v1.size(); i++) {
    v1[i] *= v2[i];
  }
}

template<class T>
void operator/= (vector<T> & v1, const vector<T> & v2) throw (DimensionException)
{
  for(unsigned int i = 0; i < v1.size(); i++) {
    v1[i] /= v2[i];
  }
}



template<class T, class C>
void operator+= (vector<T> & v1, const C & c)
{
  for(unsigned int i = 0; i < v1.size(); i++) {
    v1[i] += c;
  }
}

template<class T, class C>
void operator-= (vector<T> & v1, const C & c)
{
  for(unsigned int i = 0; i < v1.size(); i++) {
    v1[i] -= c;
  }
}

template<class T, class C>
void operator*= (vector<T> & v1, const C & c)
{
  for(unsigned int i = 0; i < v1.size(); i++) {
    v1[i] *= c;
  }
}

template<class T, class C> 
void operator/= (vector<T> & v1, const C & c)
{
  for(unsigned int i = 0; i < v1.size(); i++) {
    v1[i] /= c;
  }
}
/** @} */

/******************************************************************************/

class VectorTools
{
  public:
    VectorTools() {}
    virtual ~VectorTools() {}

  public:
    /**
     * @brief Build a sequence vector.
     *
     * Build a vector from a value to another with a specified step.
     * This works for numerical values for which additions, subtractions and division
     * makes sens.
     *
     * @param from The begining.
     * @param to The end.
     * @param by The step.
     * @return A vector containing the sequence.
     */
    template <class T>
    static vector<T> seq(T from, T to, T by)
    {
      vector<T> v;
      if (from < to) {
        // for (T i = from ; i <= to ; i += by) {           // Not good for double, 'to'
        for (T i = from ; i <= to + (by / 100) ; i += by) { // must be a little bit larger
          v.push_back(i);
        }
      } else {
        for (T i = from ; i >= to - (by / 100) ; i -= by) {
          v.push_back(i);
        }
      }
      return v;
    }

    /**
     * @brief Send the position of the first occurence of 'which'.
     *
     * Comparisons are performed using the == operator.
     * Maximum complexity: O(v.size()).
     *
     * @param v The vector to search.
     * @param which The element to search.
     * @return The position of which in v.
     */
    template<class T>
    static unsigned int which(const vector<T> & v, const T & which) throw (ElementNotFoundException<T>)
    {
      for(unsigned int i = 0; i < v.size(); i++)
        if(v[i] == which) return i;
      throw ElementNotFoundException<T>("VectorTools::which.", &v, &which);
    }

    /**
     * @brief Send the positions of all occurences of 'which'.
     *
     * Comparisons are performed using the == operator.
     * Complexity: O(v.size()).
     *
     * @param v The vector to search.
     * @param which The element to search.
     * @return A vector containing the positions of which in v.
     */
    template<class T>
    static vector<unsigned int> whichAll(const vector<T> & v, const T & which) throw (ElementNotFoundException<T>)
    {
      vector<unsigned int> w;
      for(unsigned int i = 0; i < v.size(); i++)
        if(v[i] == which) w.push_back(i);
      if (w.size())
        return w;
      throw ElementNotFoundException<T>("VectorTools::whichAll.", &v, &which);
    }

    /**
     * @brief Send a new vector with unique elements.
     *
     * The input vector is copied, and the copy is sorted using QuickSort algorithm.
     * A one-pass loop then look for duplicates and copy unique element to a result vector.
     * The output vector is hence sorted.
     * 
     * If v is empty, it is passed 'as is' in return (after being copied).
     *
     * @param v the vector to parse.
     */
    template<class T>
    static vector<T> unique(const vector<T> & v)
    {
      if(v.size() == 0) return v;
      vector<T> sortedV(v.begin(), v.end());
      sort(sortedV.begin(), sortedV.end());
      vector<T> uniq;
      uniq.push_back(sortedV[0]);
      for(unsigned int i = 1; i < sortedV.size(); i++)
      {
        if(sortedV[i] != sortedV[i-1]) uniq.push_back(sortedV[i]);
      }
      return uniq;
    }

    /**
     * @brief Tell if the vector as unique elements.
     *
     * The input vector is copied, and the copy is sorted using QuickSort algorithm.
     * A one-pass loop then look for duplicates.
     * 
     * If v is empty, the method returns 'true'.
     *
     * @param v the vector to parse.
     */
    template<class T>
    static bool isUnique(const vector<T> & v)
    {
      if(v.size() == 0) return true;
      vector<T> sortedV(v.begin(), v.end());
      sort(sortedV.begin(), sortedV.end());
      for(unsigned int i = 1; i < sortedV.size(); i++)
      {
        if(sortedV[i] == sortedV[i-1]) return false;
      }
      return true;
    }

  /**
   * @author Laurent Gueguen
   * @return the vector of the selected elements, in the order of the
   *  required positions
   * @param v1 the vector of elements, v2 the vector of the selected positions
   */
  template<class T>
  static vector<T> extract(const vector<T> & v1, const vector<int> & v2)
  {
    vector<T> v(v2.size());
    for (unsigned int i = 0; i < v2.size(); i++)
      v[i] = v1[v2[i]];
    return(v);
  }

    /**
     * @brief Count each element of a vector.
     *
     * @return A map with keys = unique vector values and values = count for each vector value.
     * @param v the vector to parse.
     */
    template<class T>
    static map<T, unsigned int> countValues(const vector<T> & v)
    {
      map<T, unsigned int> c;
      for (unsigned int i = 0 ; i < v.size() ; i++)
      {
        c[v[i]]++;
      }
      return c;
    }

    /**
     * @brief Get the break points for a given number of classes.
     *
     * Given a vector of values, return the values that cut the range of values
     * in a given number of classes.
     *
     * @param v The vector to parse.
     * @param n The expected number of clasees.
     * @return a vector of size = n + 1 containing the breaking points.
     */
    static vector<double> breaks(const vector<double> & v, unsigned int n);

    /**
     * @brief Get the optimal class number following Scott's method.
     *
     * Use Scott's (1979) method to compute the optimal class number for histogram.
     *
     * Scott, D.W. (1979) On optimal and data-based histograms. Biometrika, 66, 605¿610.
     *
     * @param v The vector to parse.
     * @return The number of classes.
     */
    template<class T>
    static unsigned int nclassScott(const vector<T> & v) {
      vector<T> r1 = VectorTools::range(v);
      T r = r1[1] - r1[0];
      double n = v.size();
      double h = 3.5 * VectorTools::sd<T, double>(v) * NumTools::pow(n, -1. / 3);
      return (unsigned int) ceil(r / h);
    }

    /**
     * @return The product of all elements in a vector.
     * @param v1 A vector.
     */
    template<class T>
    static T prod(const vector<T> & v1)
    {
      T p = 1;
      for(unsigned int i = 0; i < v1.size(); i++) p *= v1[i];
      return p;
    }

    /**
     * @return The sum of all elements in a vector.
     * @param v1 A vector.
     */
    template<class T>
    static T sum(const vector<T> & v1)
    {
      T p = 0;
      for(unsigned int i = 0; i < v1.size(); i++) p += v1[i];
      return p;
    }

    /**
     * @author Laurent Gueguen
     * @Log-normalize vector v1, ie add a constant to the elements of v
     *  such that @f$\sum_i(\exp(v_i)) = 1@f$.
     * @param v vector.
     */
    template<class T>
    static void lognorm(vector<T> & v)
    {
      T M = std::max(v);
      T x = std::exp(v[0] - M);
      for (unsigned int i = 1; i < v.size(); i++)
        x += std::exp(v[i] - M);
      v -= M + log(x);
    }
  
    /**
     * @author Laurent Gueguen
     * @return From vectors v1 and v2, return @f$\log(\sum_i(v2_i*\exp(v1_i)))@f$.
     * @param v1 and v2 two vectors.
     */
    template<class T>
    static T logsumexp(const vector<T> & v1, const vector<T> & v2)
    {
      unsigned int size;
      if(v1.size() != v2.size())
        throw DimensionException("VectorOperators::logsumexp", v1.size(), v2.size());
      else
        size = v1.size();
    
      T M = std::max(v1);
      T x = v2[0] * std::exp(v1[0] - M);
      for (unsigned int i = 1; i < size; i++)
        x += v2[i] * std::exp(v1[i] - M);
      return(std::log(x) + M);
    }

    
    /**
     * @name These methods apply the corresponding function to each element
     * and return the result in a new vector.
     *
     * @{
     */
    template<class T>
    static vector<double> log(const vector<T> & v1)
    {
      vector<double> v2(v1.size());
      for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::log(v1[i]);
      return v2;
    }
    template<class T>
    static vector<double> log(const vector<T> & v1, double base)
    {
      vector<double> v2(v1.size());
      for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::log(v1[i]) / std::log(base);
      return v2;
    }

    template<class T>
    static vector<double> exp(const vector<T> & v1)
    {
      vector<double> v2(v1.size());
      for(unsigned int i = 0; i < v2.size(); i++) v2[i] = std::exp(v1[i]);
      return v2;
    }

    template<class T>
    static vector<double> log10(const vector<T> & v1)
    {
      vector<double> v2(v1.size());
      for(unsigned int i = 0; i < v1.size(); i++) v2[i] = std::log10(v1[i]);
      return v2;
    }

    template<class T>
    static vector<T> fact(const vector<T> & v1)
    {
      vector<T> v2(v1.size());
       for(unsigned int i = 0; i < v1.size(); i++) v2[i] = NumTools::fact<T>(v1[i]);
      return v2;
    }

    template<class T>
    static vector<T> sqr(const vector<T> & v1)
    {
      vector<T> v2(v1.size());
       for(unsigned int i = 0; i < v1.size(); i++) v2[i] = NumTools::sqr<T>(v1[i]);
      return v2;
    }

    template<class T>
    static vector<T> pow(const vector<T> & v1, T & b)
    {
      vector<T> v2(v1.size());
       for(unsigned int i = 0; i < v1.size(); i++) v2[i] = NumTools::pow<T>(v1[i], b);
      return v2;
    }
    /** @} */

    /**
     * @brief Concatenate a vector after converting to string.
     *
     * @param v The vector to concatenate.
     * @param delim A string which is used to separate the values (default is " ").
     */
    template<class T>
    static string paste(const vector<T> & v, const string & delim = " ")
    {
      ostringstream out;
      for(unsigned int i = 0; i < v.size(); i++)
      {
        out << v[i];
        if (i < v.size() - 1)
          out << delim;
      }
      return out.str();
    }

    /**
     * @brief Print a vector to a stream.
     * @param v1 A vector.
     * @param out A stream.
     * @param delim A string which is used to separate the values (default is " ").
     */
    template<class T>
    static void print(const vector<T> & v1, ostream & out = cout, const string & delim = " ")
    {
      for(unsigned int i = 0; i < v1.size(); i++)
      {
        out << v1[i];
        if (i < v1.size() - 1)
          out << delim;
      }
      out << endl;
    }

    /**
     * @return The scalar product of two vectors.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @throw DimensionException If the two vector do not have the same length.
     */
    template<class InputType, class OutputType>
    static OutputType scalar(const vector<InputType> & v1, const vector<InputType> & v2) throw (DimensionException)
    {
      if(v1.size() != v2.size())
      {
        throw DimensionException("VectorFunctions::scalar", v1.size(), v2.size());
      }
      OutputType result = 0;  
      for(unsigned int i = 0; i < v1.size(); i++)
      {
        result += v1[i] * v2[i];
      }
      return result;
    }
    /**
     * This dt product correspond to the dot product <v1,v2> in the space defined by
     * @f[
     * M =
     * \begin{pmatrix}
     * w_1 & \ldots & \\
     * \vdots & w_2  & \ldots\\
     *        & \vdots & \ddots\\
     * \end{pmatrix}
     * @f]
     * @return The "weighted" scalar product of two vectors.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @param w A vector of weights.
     * @throw DimensionException If the two vector do not have the same length or do not match the length of the weights.
     */
    template<class InputType, class OutputType>
    static OutputType scalar(const vector<InputType> & v1, const vector<InputType> & v2, const vector<InputType> & w) throw (DimensionException)
    {
      if(v1.size() != w.size())
      {
        throw DimensionException("VectorFunctions::scalar", v1.size(), w.size());
      }
      if(v2.size() != w.size())
      {
        throw DimensionException("VectorFunctions::scalar", v2.size(), w.size());
      }
      OutputType result = 0;  
      for(unsigned int i = 0; i < v1.size(); i++)
      {
        result += v1[i] * v2[i] * w[i];
      }
      return result;
    }

    /**
     * @return The scalar Kronecker product of two vectors.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @throw DimensionException If the two vector do not have the same length.
     */
    template<class T>
    static vector<T> kroneckerMult(const vector<T> & v1, const vector<T> & v2) throw (DimensionException)
    {
      unsigned int n1 = v1.size();
      unsigned int n2 = v2.size();
      vector<T> v3(n1*n2);
      for(unsigned int i = 0; i < n1; i++)
      {
        T v1i = v1[i];
        for(unsigned int j = 0; j < n2; j++)
        {
          v3[i * n2 + j] = v1i * v2[j];
        }
      }
      return v3;
    }

    /**
     * @return The norm of a vector (@f$\sqrt{\sum_i^n x_i^2}@f$).
     * @param v1 A vector.
     */
    template<class InputType, class OutputType>
    static OutputType norm(const vector<InputType> & v1)
    {
      OutputType result = 0;
      for(unsigned int i = 0; i < v1.size(); i++)
        result += v1[i] * v1[i];
      return sqrt(result);
    }
    
    /**
     * @return The "weighted" norm of a vector (@f$\sqrt{\sum_i^n x_i^2}@f$).
     * @param v1 A vector.
     * @param w A vector of weights.
     * @throw DimensionException If v1 and w do not have the same length.
     * @see scalar.
     */
    template<class InputType, class OutputType>
    static OutputType norm(const vector<InputType> & v1, const vector<InputType> & w) throw (DimensionException)
    { 
      if(v1.size() != w.size())
      {
        throw DimensionException("VectorFunctions::norm", v1.size(), w.size());
      }
      OutputType result = 0;
      for(unsigned int i = 0; i < v1.size(); i++)
        result += v1[i] * v1[i] * w[i];
      return sqrt(result);
    }
    
    /**
     * @return The cosinus of the angle of two vectors.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @throw DimensionException If the two vector do not have the same length.
     */
    template<class InputType, class OutputType>
    static OutputType cos(const vector<InputType> & v1, const vector<InputType> & v2) throw (DimensionException)
    {
      return scalar<InputType, OutputType>(v1, v2)
        / (norm<InputType, OutputType>(v1) * norm<InputType, OutputType>(v2));
    }

    /**
     * @return The weighted cosinus of the angle of two vectors.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @param w A vector of weights.
     * @throw DimensionException If the two vector do not have the same length.
     */
    template<class InputType, class OutputType>
    static OutputType cos(const vector<InputType> & v1, const vector<InputType> & v2, const vector<InputType> & w) throw (DimensionException)
    {
      return scalar<InputType, OutputType>(v1, v2, w)
        / (norm<InputType, OutputType>(v1, w) * norm<InputType, OutputType>(v2, w));
    }

    /**
     * @name Extrema.
     *
     * @{
     */
 
    /**
     * @brief Template function to get the minimum value of a vector.
     *
     * The < operator must be defined for the specified class.
     *
     * @param v The input vector.
     * @return The minimum value in the vector.
     * @throw EmptyVectorException If the input vector is empty.
     */
    template<class T>
    static T min(const vector<T> & v) throw (EmptyVectorException<T>)
    {
      if(v.size() == 0) throw EmptyVectorException<T>("VectorFunctions::min()", & v);
      T mini = v[0];
      for(unsigned int i = 1; i < v.size(); i++)
        if(v[i] < mini) mini = v[i];
      return mini;
    }

    /**
     * @brief Template function to get the maximum value of a vector.
     *
     * The > operator must be defined for the specified class.
     *
     * @param v The input vector.
     * @return The maximum value in the vector.
     * @throw EmptyVectorException If the input vector is empty.
     */
    template<class T>
    static T max(const vector<T> & v) throw (EmptyVectorException<T>)
    {
      if(v.size() == 0) throw EmptyVectorException<T>("VectorFuntions::max()", & v);
      T maxi = v[0];
      for(unsigned int i = 1; i < v.size(); i++)
        if(v[i] > maxi) maxi = v[i];
      return maxi;
    }

    /**
     * @brief Template function to get the index of the minimum value of a vector.
     *
     * The < operator must be defined for the specified class.
     * The position sent is the first one matching the minimum value.
     *
     * @param v The input vector.
     * @return The position of the minimum value in the vector.
     * @throw EmptyVectorException If the input vector is empty.
     */
    template<class T>
    static unsigned int posmin(const vector<T> & v) throw (EmptyVectorException<T>)
    {
      T mini = min(v);
      for(unsigned int i = 0; i < v.size(); i++)
        if(v[i] == mini) return i;
      // This is never reached but must be here, otherwise a warning is issued:
      return 0;
    }

    /**
     * @brief Template function to get the index of the maximum value of a vector.
     *
     * The > operator must be defined for the specified class.
     * The position sent is the first one matching the maximum value.
     *
     * @param v The input vector.
     * @return The position of the maximum value in the vector.
     * @throw EmptyVectorException If the input vector is empty.
     */
    template<class T>
    static unsigned int whichmax(const vector<T> & v) throw (EmptyVectorException<T>)
    {
      T maxi = max(v);
      for(unsigned int i = 0; i < v.size(); i++)
        if(v[i] == maxi) return i;
      // This is never reached but must be here, otherwise a warning is issued:
      return 0;
    }

    /**
     * @brief Template function to get both extrema of a vector.
     *
     * Both < and > operators must be defined for the specified class.
     *
     * @param v The input vector.
     * @return A vector of size 2 which values are min(v) and max(v).
     * throw EmptyVectorException If the input vector is empty.
     */
    template<class T>
    static vector<T> range(const vector<T> & v) throw (EmptyVectorException<T>)
    {
      if (v.size() == 0) throw EmptyVectorException<T>("VectorFuntions::range()", & v);
      vector<T> r(2);
      r[0] = r[1] = v[0];
      for (unsigned int i = 1; i < v.size(); i++) {
        if(v[i] < r[0]) r[0] = v[i];
        if(v[i] > r[1]) r[1] = v[i];
      }
      return r;
    }

    /** @} */

    /**
     * @return The mean value of the vector.
     * @param v1 A vector.
     */
    template<class InputType, class OutputType>
    static OutputType mean(const vector<InputType> & v1)
    { 
      return (OutputType)sum<InputType>(v1) / (OutputType)v1.size();
    }
    /**
     * @return The weighted mean value of the vector.
     * @param v1 A vector.
     * @param w A vector of weights.
     * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
     */
    template<class InputType, class OutputType>
    static OutputType mean(const vector<InputType> & v1, const vector<InputType> & w, bool normalizeWeights = true)
    { 
      if(normalizeWeights) 
      {
        vector<InputType> wn = w / sum(w);
        return scalar<InputType, OutputType>(v1, wn);
      }
      else
      {
        return scalar<InputType, OutputType>(v1, w);
      }
    }

    /**
     * @return The median value of the vector.
     * @param v1 A vector.
     */
    template<class InputType>
    static InputType median(vector<InputType> & v1)
    {
      InputType med = 0;
      if(v1.size() == 0) return med;
      if(v1.size() == 1) return v1[0];
      sort(v1.begin(), v1.end());
      unsigned int i = v1.size() / 2;
      if(v1.size() % 2 == 0)
      {
        //Vector size is pair
        med = double((v1[i-1] + v1[i]) / 2);
      }
      else
      {
        //Vector size is impair
        med = v1[i];
      }
      return med;
    }

    /**
     * @brief Set the mean of a vector to be 0.
     * 
     * @return A vector with mean 0.
     * @param v1 A vector.
     */
    template<class InputType, class OutputType>
    static vector<OutputType> center(const vector<InputType> & v1)
    { 
      OutputType m = mean<InputType,OutputType>(v1);
      vector<OutputType> v(v1.size());
      for(unsigned int i = 0; i < v1.size(); i++)
      {
        v[i] = (OutputType)v1[i] - m;
      } 
      return v;
    }
    /**
     * @brief Set the weighted mean of a vector to be 0.
     * 
     * @return A vector with mean 0.
     * @param v1 A vector.
     * @param w A vector of weights.
     * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
     */
    template<class InputType, class OutputType>
    static vector<OutputType> center(const vector<InputType> & v1, const vector<InputType> & w, bool normalizeWeights = true)
    { 
      OutputType m = mean<InputType, OutputType>(v1, w, normalizeWeights);
      vector<OutputType> v(v1.size());
      for(unsigned int i = 0; i < v1.size(); i++)
      {
        v[i] = (OutputType)v1[i] - m;
      } 
      return v;
    }
 
    /**
     * @return The covariance of two vectors.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @param unbiased Tell if an unbiased estimate must be computed.
     * @throw DimensionException If the two vector do not have the same length.
     */
    template<class InputType, class OutputType>
    static OutputType cov(const vector<InputType> & v1, const vector<InputType> & v2, bool unbiased = true) throw (DimensionException)
    {
      OutputType n = (OutputType)v1.size();
      OutputType x =  scalar<InputType,OutputType>(
          center<InputType, OutputType>(v1),
          center<InputType, OutputType>(v2)
          ) / n;
      if(unbiased) x = x * n / (n - 1);
      return x;
    }

    /**
     * @return The weighted covariance of two vectors.
     * To have a population estimate you have to multiply by \f$\frac{n}{n-1}\f$.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @param w A vector of weights.
     * @param unbiased Tell if an unbiased estimate must be computed.
     * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
     * @throw DimensionException If the two vector do not have the same length.
     */
    template<class InputType, class OutputType>
    static OutputType cov(const vector<InputType> & v1, const vector<InputType> & v2, const vector<InputType> & w, bool unbiased = true, bool normalizeWeights = true) throw (DimensionException)
    { 
      if(normalizeWeights) 
      {
        vector<InputType> wn = w / sum(w);
        OutputType x = scalar<InputType, OutputType>(
            center<InputType, OutputType>(v1, wn, false),
            center<InputType, OutputType>(v2, wn, false),
            wn
          );
        if(unbiased)
        {
          x = x / (1 - sum(sqr<double>(wn)));
        }
        return x;
      }
      else
      {
         OutputType x = scalar<InputType, OutputType>(
            center<InputType, OutputType>(v1, w, false),
            center<InputType, OutputType>(v2, w, false),
            w
          );
        if(unbiased)
        {
          x = x / (1 - sum(sqr(w)));
        }
        return x;
      }
    }
    /**
     * @return The variance of the vector.
     * @param v1 The sample vector.
     * @param unbiased Tell if an unbiased estimate must be computed.
     */
    template<class InputType, class OutputType>
    static OutputType var(const vector<InputType> & v1, bool unbiased = true)
    {
      return cov<InputType, OutputType>(v1, v1, unbiased);
    }
    /**
     * @return The weighted variance of the vector.
     * @param v1 The sample vector.
     * @param w A vector of weights.
     * @param unbiased Tell if an unbiased estimate must be computed.
     * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
     * @throw DimensionException If v1 and w do not have the same length.
     */
    template<class InputType, class OutputType>
    static OutputType var(const vector<InputType> & v1, const vector<InputType> & w, bool unbiased = true, bool normalizeWeights = true) throw (DimensionException)
    {
      return cov<InputType, OutputType>(v1, v1, w, unbiased, normalizeWeights);
    }

    /**
     * @return The standard deviation of the vector.
     * @param v1 The sample vector.
     * @param unbiased Tell if an unbiased estimate must be computed.
     */
    template<class InputType, class OutputType>
    static OutputType sd(const vector<InputType> & v1, bool unbiased = true)
    {
      return sqrt(var<InputType, OutputType>(v1, unbiased));
    }

    /**
     * @return The weighted standard deviation of the vector.
     * @param v1 The sample vector.
     * @param w A vector of weights.
     * @param unbiased Tell if an unbiased estimate must be computed.
     * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
     * @throw DimensionException If v1 and w do not have the same length.
     */
    template<class InputType, class OutputType>
    static OutputType sd(const vector<InputType> & v1, const vector<InputType> & w, bool unbiased = true, bool normalizeWeights = true) throw (DimensionException)
    {
      return sqrt(var<InputType, OutputType>(v1, w, unbiased, normalizeWeights));
    }

    /**
     * @return The Pearson correlation coefficient of two vectors.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @throw DimensionException If the two vector do not have the same length.
     */
    template<class InputType, class OutputType>
    static OutputType cor(const vector<InputType> & v1, const vector<InputType> & v2) throw (DimensionException)
    {
      return cov<InputType, OutputType>(v1, v2)
        / ( sd<InputType, OutputType>(v1) * sd<InputType, OutputType>(v2) );
    }

    /**
     * @return The weighted Pearson correlation coefficient of two vectors.
     * @param v1 First vector.
     * @param v2 Second vector.
     * @param w A vector of weights.
     * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
     * @throw DimensionException If the two vector do not have the same length.
     */
    template<class InputType, class OutputType>
    static OutputType cor(const vector<InputType> & v1, const vector<InputType> & v2, const vector<InputType> & w, bool normalizeWeights = true) throw (DimensionException)
    {
      if(normalizeWeights) 
      {
        vector<InputType> wn = w / sum(w);
        return cov<InputType, OutputType>(v1, v2, wn, false, false)
          / ( sd<InputType, OutputType>(v1, wn, false, false) * sd<InputType, OutputType>(v2, wn, false, false) );
      }
      else
      {
        return cov<InputType, OutputType>(v1, v2, w, false, false)
          / ( sd<InputType, OutputType>(v1, w, false, false) * sd<InputType, OutputType>(v2, w, false, false) );
      }
    }

    /**
     * @return The Shannon entropy indice of the vector.
     * @param v1 The vector.
     * @param base The base of the logarithm to use.
     */
    template<class InputType, class OutputType>
    static double shannon(const vector<InputType> & v, double base = 2.7182818)
    {
      OutputType s = 0;
      for(unsigned int i = 0; i < v.size(); i++)
        if(v[i] > 0) s += v[i] * std::log(v[i]) / std::log(base);
      return -s;
    }

    /**
     * @return 'true' if the two vectors contains the same elements, whatever their order in the container.
     * @param v1 First vector.
     * @param v2 Second vector.
     */
    template<class T>
    static bool haveSameElements(const vector<T> & v1, const vector<T> & v2)
    {
      vector<T> u1(v1);
      vector<T> u2(v2);
      if(u1.size() != u2.size()) return false;
      std::sort(u1.begin(), u1.end());
      std::sort(u2.begin(), u2.end());
      return (u1 == u2);
    }

    /**
     * @return 'true' if the two vectors contains the same elements, <b>in the same frequency</b>, whatever their order in the container.
     *
     * @warning The two input vectors will be sorted.
     *
     * @param v1 First vector.
     * @param v2 Second vector.
     */
    template<class T>
    static bool haveSameElements(vector<T> & v1, vector<T> & v2)
    {
      if(v1.size() != v2.size()) return false;
      std::sort(v1.begin(), v1.end());
      std::sort(v2.begin(), v2.end());
      return (v1 == v2);
    }

    /**
     * @return 'true' if a the input vector contains the given element.
     * @param vec The vector to check.
     * @param el The element to look for.
     */
    template<class T>
    static bool contains(const vector<T> & vec, T el)
    {
      for(unsigned int i = 0; i < vec.size(); i++)
        if(vec[i] == el) return true;
      return false;
    }

    /**
     * @return 'true' if a the first vector contains all elements of the second vector.
     *
     * @warning The two input vectors will be sorted.
     *
     * @param v1 The first vector to check.
     * @param v2 The second vector to check.
     */
    template<class T>
    static bool containsAll(vector<T> & v1, vector<T> & v2)
    {
      std::sort(v1.begin(), v1.end());
      std::sort(v2.begin(), v2.end());
      unsigned int j = 0;
      for(unsigned int i = 0; i < v2.size(); i++)
      {
        if(i > 0 && v2[i] == v2[i-1]) continue;
        while(j < v1.size() - 1 && v1[j] < v2[i]) j++;
        if(v1[j] != v2[i]) return false;
      }
      return true;
    }

    /**
     * @return A vector which is the union of two vectors passed as input.
     * Duplicate element will be removed.
     * @param vec1 Vector 1.
     * @param vec2 Vector 2.
     */
    template<class T>
    static vector<T> vectorUnion(const vector<T> & vec1, const vector<T> & vec2)
    {
      vector<T> unionEl = vec1;
      for(unsigned int j = 0; j < vec2.size(); j++)
      {
        if(!contains(unionEl, vec2[j]))
        unionEl.push_back(vec2[j]);
      }
      return unionEl;
    }

    /**
     * @return A vector which is the union of all vectors passed as input.
     * Duplicate element will be removed.
     * @param vecElementL A vector of vectors.
     */
    template<class T>
    static vector<T> vectorUnion(const vector< vector<T> > & vecElementL)
    {
      vector<T> unionEl;
      for(unsigned int i = 0; i < vecElementL.size(); i++)
      {
        for(unsigned int j = 0; j < vecElementL[i].size(); j++)
        {
          if(!contains(unionEl, vecElementL[i][j]))
          unionEl.push_back(vecElementL[i][j]);
        }
      }
      return unionEl;
    }

    /**
     * @return A vector which is the intersection of two vectors passed as input.
     * @param vec1 Vector 1.
     * @param vec2 Vector 2.
     */
    template<class T>
    static vector<T> vectorIntersection(const vector<T> & vec1, const vector<T> & vec2)
    {
      vector<T> interEl;
      for(unsigned int i = 0; i < vec1.size(); i++)
      {
        if(contains(vec2, vec1[i])) interEl.push_back(vec1[i]);
      }
      return interEl;
    }

    /**
     * @return A vector which is the intersection of all vectors passed as input.
     * @param vecElementL A vector of vectors.
     */
    template<class T>
    static vector<T> vectorIntersection(const vector< vector<T> > & vecElementL)
    {
      if(vecElementL.size() == 1) return vecElementL[0];
      vector<T> interEl;
      if(vecElementL.size() == 0) return interEl;
      for(unsigned int i = 0; i < vecElementL[0].size(); i++)
      {
        bool test = true;
        for(unsigned int j = 1; test && j < vecElementL.size(); j++)
        {
          if(!contains(vecElementL[j], vecElementL[0][i])) test = false;
        }
        if(test) interEl.push_back(vecElementL[0][i]);
      }
      return interEl;
    }

    /**
     * @brief Append the content of a vector to another one.
     * @param vec1 Vector 1.
     * @param vec2 Vector 2.
     */
    template<class T>
    static void append(vector<T> & vec1, const vector<T> & vec2)
    {
      for(unsigned int i = 0; i < vec2.size(); i++)
      {
        vec1.push_back(vec2[i]);
      }
    }

    /**
     * @return A single vector made of the concatenation of the vectors passed as input.
     * @param vecElementL A vector of vectors.
     */
    template<class T>
    static vector<T> append(const vector< vector<T> > & vecElementL)
    {
      if(vecElementL.size() == 1) return vecElementL[0];
      vector<T> v;
      if(vecElementL.size() == 0) return v;
      for(unsigned int i = 0; i < vecElementL[0].size(); i++)
      {
        v.push_back(vecElementL[0][i]);
      }
      return v;
    }

    /**
     * @brief This function returns the difference of two vectors.
     *
     * @warning The two input vectors will be sorted. As a consequence, the output vector will be also sorted.
     *
     * @param v1 First vector.
     * @param v2 Second vector.
     * @param v3 A vector to be populated with all elements in v1 that are not found in v2.
     */  
    template<class T>
    static void diff(vector<T> & v1, vector<T> & v2, vector<T> & v3)
    {
      if(v2.size() == 0) append(v3, v1);
      std::sort(v1.begin(), v1.end());
      std::sort(v2.begin(), v2.end());
      unsigned int j = 0;
      for(unsigned int i = 0; i < v1.size(); i++)
      {
        if(i > 0 && v1[i] == v1[i-1]) continue;
        while(j < v2.size() - 1 && v2[j] < v1[i]) j++;
        if(v2[j] != v1[i]) v3.push_back(v1[i]);
      }
    };

    /**
     * @brief Test function.
     * @return true if all tests are passed.
     */
    static bool test();

};

} //end of namespace bpp.

#endif  //_VECTORTOOLS_H_

