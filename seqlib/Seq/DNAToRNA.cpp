//
// File: DNAToRNA.cpp
// Created by: Julien Dutheil
// Created on: Sun Oct 12 14:39:29 2003
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#include "DNAToRNA.h"

using namespace bpp;

/******************************************************************************/

DNAToRNA::DNAToRNA() : AbstractReverseTranslator()
{
	dna = new DNA();
	rna = new RNA();
}

/******************************************************************************/

DNAToRNA::~DNAToRNA()
{
	delete dna;
	delete rna;
}

/******************************************************************************/

int DNAToRNA::translate(int state) const throw (BadIntException)
{
	dna -> intToChar(state);
	return state;
}

/******************************************************************************/

string DNAToRNA::translate(const string & state) const throw (BadCharException)
{
	int i = dna -> charToInt(state);
	return rna -> intToChar(i);
}

/******************************************************************************/

int DNAToRNA::reverse(int state) const throw (BadIntException) 
{
	rna -> intToChar(state);
	return state;
}

/******************************************************************************/

string DNAToRNA::reverse(const string & state) const throw (BadCharException)
{
	int i = rna -> charToInt(state);
	return dna -> intToChar(i);
}

/******************************************************************************/

