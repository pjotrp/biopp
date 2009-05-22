//
// File: MatrixTools.h
// Created by: Julien Dutheil
// Created on: Mon Jan 19 16:42:25 2004
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

#ifndef _MATRIXTOOLS_H_
#define _MATRIXTOOLS_H_

#include "VectorTools.h"
#include "Matrix.h"
#include "LUDecomposition.h"
#include "EigenValue.h"

#include <cstdio>
#include <iostream>
using namespace std;

namespace bpp
{

/**
 * @brief Functions dealing with matrices.
 */
class MatrixTools
{
	public:
		MatrixTools() {}
		~MatrixTools() {}

	public:

		/**
     * @brief Copy operation. This function supplies the lack of inheritence of the assigment operator :D .
     *
		 * @param A [in] Original matrix.
		 * @param O [out] A copy of the given matrix.
		 */
		template<class MatrixA, class MatrixO>
		static void copy(const MatrixA & A, MatrixO & O)
    { 
      O.resize(A.nRows(), A.nCols());
      for(unsigned int i = 0; i < A.nRows(); i++)
        for(unsigned int j = 0; j < A.nCols(); j++)
          O(i, j) = A(i, j);
    }
	
		/**
     * @brief Get a identity matrix of a given size.
     *
		 * @param n the size of the matrix.
		 * @return O [out] A identity matrix of size n.
		 */
		template<class Matrix>
		static void getId(unsigned int n, Matrix & O)
		{
			O.resize(n, n);
			for(unsigned int i = 0; i < n; i++)
      {
				for(unsigned int j = 0; j < n; j++) O(i, j) = (i == j) ? 1 : 0;
			}
		}

		/**
		 * @param D [in] A vector of diagonal elements.
		 * @param O [out] A diagonal matrix with diagonal elements taken from a vector.
		 */
		template<class Scalar>
		static void diag(const vector<Scalar> & D, Matrix<Scalar> & O)
		{
			unsigned int n = D.size();
			O.resize(n, n);
			for(unsigned int i = 0; i < n; i++)
      {
				for(unsigned int j = 0; j < n; j++) O(i, j) = (i == j) ? D[i] : 0;
			}
		}

		/**
		 * @param M [in] The matrix.
		 * @param O [out] The diagonal elements of a square matrix as a vector.
		 * @throw DimensionException If M is not a square matrix.
		 */
		template<class Scalar>
		static void diag(const Matrix<Scalar> & M, vector<Scalar> & O) throw (DimensionException)
		{
			unsigned int nc = M.nCols();
			unsigned int nr = M.nRows();
			if(nc != nr) throw DimensionException("MatrixTools::diag(). M must be a square matrix.", nr, nc); 
			O.resize(nc);
			for(unsigned int i = 0; i < nc; i++) O[i] = M(i, i);
		}

		/**
		 * @brief Set all elements in M to value x.
		 * @param M A matrix.
		 * @param x The value to use.
		 */
		template<class Matrix, class Scalar>
		static void fill(Matrix & M, Scalar x)
		{
			for(unsigned int i = 0; i < M.nRows(); i++)
      {
				for(unsigned int j = 0; j < M.nCols(); j++)
        {
					M(i, j) = x;
				}
			}
		}

		/**
		 * @brief Multiply all elements of a matrix by a given value, and add a constant.
		 *
		 * Performs \f$\forall i \forall j m_{i,j} = a.m_{i,j}+b\f$.
		 * 
		 * @param A A matrix.
		 * @param a Multiplicator.
		 * @param b Constant.
		 */
		template<class Matrix, class Scalar>
		static void scale(Matrix & A, Scalar a, Scalar b = 0)
    {
			for(unsigned int i = 0; i < A.nRows(); i++)
      {
				for(unsigned int j = 0; j < A.nCols(); j++)
        {
					A(i,j) = a * A(i, j) + b;
				}
			}
		}

		/**
		 * @param A [in] First matrix.
		 * @param B [in] Second matrix.
		 * @param O [out] The dot product of two matrices.
		 */
		template<class Scalar>
		static void mult(const Matrix<Scalar> & A, const Matrix<Scalar> & B, Matrix<Scalar> & O) throw (DimensionException)
		{
			unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
			if(ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA); 
			O.resize(nrA, ncB);
			for(unsigned int i = 0; i < nrA; i++)
      {
				for(unsigned int j = 0; j < ncB; j++)
        {
					O(i, j) = 0;
					for(unsigned int k = 0; k < ncA; k++)
          {
						O(i, j) += A(i, k) * B(k, j);
					}
				}
			}
		}

		/**
		 * @brief Compute A . D . B where D is a diagonal matrix in O(n^3).
		 *
		 * Since D is a diagonal matrix, this function is more efficient than doing
		 * mult(mult(A, diag(D)), B), which involves two 0(n^3) operations.
		 *
		 * @param A [in] The first matrix.
		 * @param D [in] The diagonal matrix (only diagonal elements in a vector)
		 * @param B [in] The second matrix.
		 * @param O [out] The result matrix.
		 * @throw DimensionException If matrices have not the appropriate size.
		 */
		template<class Scalar>
		static void mult(const Matrix<Scalar> & A, const vector<Scalar> & D, const Matrix<Scalar> & B, Matrix<Scalar> & O) throw (DimensionException)
		{
			unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
			if(ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA); 
			if(ncA != D.size()) throw DimensionException("MatrixTools::mult(). Vector size is not eual to matrix size.", D.size(), ncA); 
			O.resize(nrA, ncB);
			for(unsigned int i = 0; i < nrA; i++)
      {
				for(unsigned int j = 0; j < ncB; j++)
        {
					O(i, j) = 0;
					for(unsigned int k = 0; k < ncA; k++)
          {
						O(i, j) += A(i, k) * B(k, j) * D[k];
					}
				}
			}
		}

		/**
		 * @brief Add matrix B to matrix A.
		 *
		 * @param A [in] Matrix A
		 * @param B [in] Matrix B
		 * @throw DimensionException If A and B have note the same size.
		 */
		template<class MatrixA, class MatrixB>
		static void add(MatrixA & A, const MatrixB & B) throw (DimensionException)
		{
	 		unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
			if(ncA != ncB) throw DimensionException("MatrixTools::operator+(). A and B must have the same number of colums.", ncB, ncA); 
			if(nrA != nrB) throw DimensionException("MatrixTools::operator+(). A and B must have the same number of rows.", nrB, nrA); 
			for(unsigned int i = 0; i < A.nRows(); i++)
      {
				for(unsigned int j = 0; j < A.nCols(); j++)
        {
					A(i, j) += B(i, j);
				}
			}
		}		
	
		/**
     * @brief Compute the power of a given matrix.
     *
		 * @param A [in] The matrix.
		 * @param p The number of multiplications.
		 * @param O [out]\f$\prod_{i=1}^p m\f$.
		 * If p = 0, sends the identity matrix.
		 * @throw DimensionException If m is not a square matrix.
		 */
		template<class Matrix>
		static void pow(const Matrix & A, unsigned int p, Matrix & O) throw (DimensionException)
		{
			unsigned int n = A.nRows();
			if(n != A.nCols()) throw DimensionException("MatrixTools::pow(). nrows != ncols.", A.nCols(), A.nRows()); 
			getId<Matrix>(n, O);
      RowMatrix<double> tmp;
      for(unsigned int i = 0; i < p; i++)
      {
        tmp = O;
				mult(tmp, A, O);
			}
		}
	
		/**
     * @brief Perform matrix exponentiation.
     *
     * @warning This method currently relies only on diagonalization, so it won't work if your matrix is not diagonalizable.
     * The function may be extended later to deal with other cases.
     *
		 * @param A [in] The matrix.
		 * @param O [out]\f$\prod_{i=1}^p m\f$.
		 * @throw DimensionException If m is not a square matrix.
		 */
		template<class Scalar>
		static void exp(const Matrix<Scalar> & A, Matrix<Scalar> & O) throw (DimensionException)
		{
			unsigned int n = A.nRows();
			if(n != A.nCols()) throw DimensionException("MatrixTools::exp(). nrows != ncols.", A.nCols(), A.nRows()); 
      EigenValue<Scalar> eigen(A);
      RowMatrix<Scalar> rightEV, leftEV;
      rightEV = eigen.getV();
      inv(rightEV, leftEV);
      mult(rightEV, VectorTools::exp(eigen.getRealEigenValues()), leftEV, O);
		}
	
		/**
		 * @return The position of the maximum value in the matrix.
		 * @param m The matrix.
		 */
		template<class Matrix>
		static vector<unsigned int> whichmax(const Matrix & m)
		{
			unsigned int nrows = m.nRows();
			unsigned int ncols = m.nCols();
			vector<unsigned int> pos(2);
			unsigned int imax = 0;
			unsigned int jmax = 0;
			double currentMax = log(0.);
			for(unsigned int i = 0; i < nrows; i++)
      {
				for(unsigned int j = 0; j < ncols; j++)
        {
					double currentValue = m(i, j);
					//cout << currentValue << "\t" << (currentValue > currentMax) << endl;
					if(currentValue > currentMax)
          {
						imax = i;
						jmax = j;
						currentMax = currentValue;
					}
				}
			}
			pos[0] = imax;
			pos[1] = jmax;
			return pos;
		}

		/**
		 * @return The position of the minimum value in the matrix.
		 * @param m The matrix.
		 */
		template<class Matrix>
		static vector<unsigned int> whichmin(const Matrix & m)
		{
			unsigned int nrows = m.nRows();
			unsigned int ncols = m.nCols();
			vector<unsigned int> pos(2);
			unsigned int imin = 0;
			unsigned int jmin = 0;
			double currentMin = -log(0.);
			for(unsigned int i = 0; i < nrows; i++)
      {
				for(unsigned int j = 0; j < ncols; j++)
        {
					double currentValue = m(i, j);
					if(currentValue < currentMin)
          {
						imin = i;
						jmin = j;
						currentMin = currentValue;
					}
				}
			}
			pos[0] = imin;
			pos[1] = jmin;
			return pos;
		}

		/**
		 * @brief Print a matrix to a stream.
		 * 
		 * @param m The matrix to print.
		 * @param out The stream to use.
		 */
		template<class Matrix>
		static void print(const Matrix & m, ostream & out = cout)
		{
			out << m.nRows() << "x" << m.nCols() << endl;
			out << "[" << endl;
			for(unsigned int i = 0; i < m.nRows(); i++)
      {
				out << "[";
				for(unsigned int j = 0; j < m.nCols() - 1; j++)
        {
					out << m(i, j) << ", ";
				}
				if(m.nCols() > 0) out << m(i, m.nCols() - 1) << "]" << endl;
			}
			out << "]" << endl;
		}
		
		/**
		 * @brief Print a vector to a stream.
		 * 
		 * @param v The vector to print.
		 * @param out The stream to use.
		 */
		template<class Real>
		static void print(const vector<Real> & v, ostream & out = cout)
		{
			out << v.size() << endl;
			out << "[";
			for(unsigned int i = 0; i < v.size() - 1; i++)
      {
				out << v[i] << ", ";
			}
			if(v.size() > 0) out << v[v.size() - 1];
			out << "]" << endl;
		}

		/**
		 * @return True if the matrix is a square matrix.
		 * @param A A matrix.
		 */
		template<class Matrix>
		static bool isSquare(const Matrix & A) { return A.nRows() == A.nCols(); }

		/**
		 * @param A [in] The matrix to inverse.
		 * @param O [out] The inverse matrix of A.
		 * @throw DimensionException If A is not a square matrix.
		 */
		template<class Scalar>
		static void inv(const Matrix<Scalar> & A, Matrix<Scalar> & O) throw (DimensionException)
		{
			if(! isSquare(A)) throw DimensionException("MatrixTools::inv(). Matrix A is not a square matrix.", A.nRows(), A.nCols());
			LUDecomposition<Scalar> lu(A);
			RowMatrix<Scalar> I;
      getId(A.nRows(), I);
			copy(lu.solve(I), O);
		}

		/**
		 * @param A [in] The matrix to transpose.
		 * @param O [out] The transposition of A.
		 */
		template<class MatrixA, class MatrixO>
		static void transpose(const MatrixA & A, MatrixO & O)
		{
			O.resize(A.nCols(), A.nRows());
			for(unsigned int i = 0; i < A.nCols(); i++)
      {
				for(unsigned int j = 0; j < A.nRows(); j++)
        {
					O(i, j) = A(j, i);
				}
			}
		}
	
    /**
     * @brief Compute the Kronecker product of two row matrices.
     *
     * @param A [in] The first row matrix.
     * @param B [in] The second row matrix.
     * @param O [out] The product \f$A \otimes B\f$.
     */
    template<class Scalar>
    static void kroneckerMult(const Matrix<Scalar> & A, const Matrix<Scalar> & B, Matrix<Scalar> & O)
    {
			unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
      O.resize(nrA*nrB, ncA*ncB);
      for(unsigned int ia = 0; ia < nrA; ia++)
      {
        for(unsigned int ja = 0; ja < ncA; ja++)
        {
          Scalar aij = A(ia, ja);
          for(unsigned int ib = 0; ib < nrB; ib++)
          {
            for(unsigned int jb = 0; jb < ncB; jb++)
            {
              O(ia * nrB + ib, ja * ncB + jb) = aij * B(ib,jb);
            }
          }
        }
      }
    }

    /**
     * @brief Compute the Kronecker sum of two row matrices.
     *
     * @param A [in] The first row matrix.
     * @param B [in] The second row matrix.
     * @param O [out] The product \f$A \oplus B\f$.
     */
    template<class Scalar>
    static void kroneckerSum(const Matrix<Scalar> & A, const Matrix<Scalar> & B, Matrix<Scalar> & O)
    {
			unsigned int ncA = A.nCols();
			unsigned int nrA = A.nRows();
			unsigned int nrB = B.nRows();
			unsigned int ncB = B.nCols();
      O.resize(nrA + nrB, ncA + ncB);
      for(unsigned int ia = 0; ia < nrA; ia++)
      {
        for(unsigned int ja = 0; ja < ncA; ja++)
        {
          O(ia, ja) = A(ia,ja);
        }
      }
      for(unsigned int ib = 0; ib < nrB; ib++)
      {
        for(unsigned int jb = 0; jb < nrB; jb++)
        {
          O(nrA + ib, ncA + jb) = B(ib,jb);
        }
      }
    }

    /**
     * @brief Compute the Kronecker sum of n row matrices.
     *
     * @param vA [in] A vector of row matrices of any size.
     * @param O [out] The product \f$\bigoplus_i A_i\f$.
     */
    template<class Scalar>
    static void kroneckerSum(const vector< Matrix<Scalar> *> & vA, Matrix<Scalar> & O)
    {
			unsigned int nr = 0;
			unsigned int nc = 0;
      for(unsigned int k = 0; k < vA.size(); k++)
      {
        nr += vA[k]->nRows();
        nc += vA[k]->nCols();
      }
      O.resize(nr, nc);
      unsigned int rk = 0; //Row counter
      unsigned int ck = 0; //Col counter
      for(unsigned int k = 0; k < vA.size(); k++)
      {
        const Matrix<Scalar> * Ak = vA[k];
        for(unsigned int i = 0; i < Ak->nRows(); i++)
        {
          for(unsigned int j = 0; j < Ak->nCols(); j++)
          {
            O(rk + i, ck + j) = (*Ak)(i,j);
          }
        }
        rk += Ak->nRows();
        ck += Ak->nCols();
      }
    }

	
};

/* DEPRECATED 
namespace MatrixOperators {
	
	MatrixB operator*(const MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
		return MatrixTools::mult<MatrixA, MatrixB>(A, B);
	}
	
	template<class MatrixA, class MatrixB>
	MatrixA operator+(const MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
		MatrixA C = A;
		MatrixTools::add<MatrixA, MatrixB>(C, B);
		return C;
	}

	template<class MatrixA, class MatrixB>
	MatrixA operator+=(MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
		MatrixTools::add<MatrixA, MatrixB>(A, B);
		return A;
	}

	template<class Matrix>
	Matrix operator-(const Matrix A)
	{
		Matrix B(A.nRows(), A.nCols());
		for(unsigned int i = 0; i < B.nRows(); i++) {
			for(unsigned int j = 0; j < B.nCols(); j++) {
				B(i, j) = -A(i, j);
			}
		}
		return B;
	}

	template<class MatrixA, class MatrixB>
	MatrixA operator-(const MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
//	  unsigned int ncA = A.nCols();
//		unsigned int nrA = A.nRows();
//		unsigned int nrB = B.nRows();
//		unsigned int ncB = B.nCols();
//		if(ncA != ncB) throw DimensionException("MatrixTools::operator-(). A and B must have the same number of colums.", ncB, ncA); 
//		if(nrA != nrB) throw DimensionException("MatrixTools::operator-(). A and B must have the same number of rows.", nrB, nrA); 
//		MatrixB C(A.nRows(), A.nCols());
//		for(unsigned int i = 0; i < A.nRows(); i++) {
//			for(unsigned int j = 0; j < A.nCols(); j++) {
//				C(i, j) = A(i, j) - B(i, j);
//			}
//		}
//		return C;
		MatrixA C = A;
		MatrixTools::add<MatrixA, MatrixB>(C, -B);
		return C;
	}
	
	template<class MatrixA, class MatrixB>
	MatrixA operator-=(MatrixA & A, const MatrixB & B) throw (DimensionException)
	{
		MatrixTools::add<MatrixA, MatrixB>(A, -B);
		return A;
	}

};
*/

} //end of namespace bpp.

#endif	//_MATRIXTOOLS_H_
