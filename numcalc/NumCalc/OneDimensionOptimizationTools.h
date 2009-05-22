//
// File: OneDimensionOptimizationTools.h
// Created by: Julien Dutheil
// Created on: Mon Nov 17 11:15:22 2003
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

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

#ifndef _ONEDIMENSIONOPTIMIZATIONTOOLS_H_
#define _ONEDIMENSIONOPTIMIZATIONTOOLS_H_

#include "Functions.h"
#include "DirectionFunction.h"

namespace bpp
{

class Point
{		
	public: // Constructor and destructor:
		Point() {}
		Point(double x, double f): x(x), f(f) {}
		virtual ~Point() {}
		
	public:
		void set(double x, double f);
		
	public:
		double x;
		double f;			
};
		
class Bracket
{	
	public: // Constructor and destructor::
		Bracket() {}
		virtual ~Bracket() {}
			
	public: // Methods:
		void setA(double xa, double fa);
		void setB(double xb, double fb);
		void setC(double xc, double fc);
			
	public:
		Point a, b, c;
};

/**
 * @brief Tools of one parameter-functions optimizations.
 *
 * For now, contains only one method to bracket a minimum.
 */
class OneDimensionOptimizationTools
{
	public:
		OneDimensionOptimizationTools() {}
		virtual ~OneDimensionOptimizationTools() {}

	public:
		
		/**
		 * @brief Bracket a minimum.
		 *
		 * Given a function func, and given distinct initial points x1 and x2,
		 * this routine searches in the downhill direction (defined by the function as
		 * evaluated at the initial points) and returns a Bracket object with new points
		 * a.x, b.x and c.x that bracket a minimum of the function. Also returned are the
		 * function values at the three points, a.f, b.f and c.f.
		 *
		 * @param a, b       Two initial values for the parameter.
		 * @param function   The function to bracket.
		 * @param parameters The parameter to use as a variable.
		 * @return           A bracket object.
		 */
		static Bracket bracketMinimum(double a, double b, Function * function, ParameterList parameters);
	
    static unsigned int lineMinimization(DirectionFunction & f1dim, ParameterList & parameters, vector<double> & xi, double tolerance, ostream * profiler = NULL, ostream * messenger = NULL, int verbose = 2);

	public:
		
		/**
		 * @brief Maximum magnification allowed for a parabolic-fit step.
		 */
		static double GLIMIT;
	
};

} //end of namespace bpp.

#endif	//_ONEDIMENSIONOPTIMIZATIONTOOLS_H_

