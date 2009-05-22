//
// File: VectorExceptions.h
// Created by: Julien Dutheil
// Created on: 2003
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

#ifndef _VECTOREXCEPTIONS_H_
#define _VECTOREXCEPTIONS_H_

//From utils:
#include <Utils/Exceptions.h>
#include <Utils/TextTools.h>

//From the STL/
#include <vector>
using namespace std;

namespace bpp
{

/**
 * @brief General Exception dealing with vectors.
 */
template<class T>
class VectorException : public Exception
{

	protected:
		const vector<T> * _vect;
			
	public:
		VectorException(const char *   text, const vector<T> * vect = NULL) :
			Exception("VectorException: " + string(text)),
			_vect(vect) {};

		VectorException(const string & text, const vector<T> * vect = NULL) :
			Exception("VectorException: " + text),
			_vect(vect) {};
	
		virtual ~VectorException() throw () {};
		
	public:
		virtual const vector<T> * getVector() const { return _vect; }
};

/**
 * @brief Exception thrown when an empty vector was found.
 */
template<class T>
class EmptyVectorException : public VectorException<T>
{

	public:
		EmptyVectorException(const char *   text, const vector<T> * vect = NULL) :
			VectorException<T>("EmptyVectorException: " + string(text), vect) {};
	
		EmptyVectorException(const string & text, const vector<T> * vect = NULL) :
			VectorException<T>("EmptyVectorException: " + text, vect) {};
	
		virtual ~EmptyVectorException() throw () {}
};

/**
 * @brief Exception thrown when a dimension problem occured.
 */
class DimensionException : public Exception
{
	protected:
		unsigned int _dimension;
		unsigned int _correctDimension; 
			
	public:
		DimensionException(const char *   text, unsigned int dimension, unsigned int correctDimension) :
			  Exception("DimensionException (found " + TextTools::toString(dimension) +", should be " + TextTools::toString(correctDimension) + string(text)),
        _dimension(dimension),
        _correctDimension(correctDimension) {}
				
		DimensionException(const string & text, unsigned int dimension, unsigned int correctDimension) :
		    Exception("DimensionException (found " + TextTools::toString(dimension) + ", should be " + TextTools::toString(correctDimension) + text),
        _dimension(dimension),
        _correctDimension(correctDimension) {};
	
		virtual ~DimensionException() throw () {}

	public:
		virtual unsigned int getDimension() const { return _dimension; }
		virtual unsigned int getCorrectDimension() const { return _correctDimension; }
};

/**
 * @brief Exception thrown when a given element was not found in the vector.
 */
template<class T> class ElementNotFoundException : public VectorException<T>
{

	protected:
		const T * _element;
			
	public:
		ElementNotFoundException(const char *   text, const vector<T> * vect = NULL, const T * element = NULL) :
			VectorException<T>("ElementNotFoundException: " + string(text), vect),
			_element(element) {};

		ElementNotFoundException(const string & text, const vector<T> * vect = NULL, const T * element = NULL) :
			VectorException<T>("ElementNotFoundException: " + text, vect),
			_element(element) {};

		virtual ~ElementNotFoundException() throw () {};
		
	public:
		virtual const T * getElement() const { return _element; }
};

} //end of namespace bpp.

#endif //_VECTOREXCEPTIONS_H_

