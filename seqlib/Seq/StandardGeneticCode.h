//
// File: StandardGeneticCode.h
// Created by: Julien Dutheil
// Created on: Mon Oct 13 15:39:17 2003
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

#ifndef _STANDARDGENETICCODE_H_
#define _STANDARDGENETICCODE_H_

#include "GeneticCode.h"
#include "NucleicAlphabet.h"

using namespace std;

namespace bpp
{

/**
 * @brief This class implements the standard genetic code as describe on the NCBI 
 *        web site: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=t#SG1
 */

class StandardGeneticCode:
  public GeneticCode
{
	public:
		StandardGeneticCode(const NucleicAlphabet * alpha);
		virtual ~StandardGeneticCode();
	
	public:
		int    translate(           int state) const throw (Exception);
		string translate(const string & state) const throw (Exception);
		Sequence * translate(const Sequence & sequence) const throw (Exception)
    {
			return GeneticCode::translate(sequence);	
		}
};

} //end of namespace bpp.

#endif	//_STANDARDGENETICCODE_H_

