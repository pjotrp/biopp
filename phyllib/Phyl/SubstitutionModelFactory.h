//
// File: SubstitutionModelFactory.h
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Fri apr 14 11:11 2006
//

/*
Copyright or � or Copr. CNRS, (November 16, 2004, 2005, 2006)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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

#ifndef _SUBSTITUTIOMODELFACTORY_H_
#define _SUBSTITUTIOMODELFACTORY_H_

#include "models"
#include "Tree.h"

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/AlphabetExceptions.h>

// From the STL:
#include <string>

using namespace std;

namespace bpp
{

/**
 * @brief Utilitary class for creating substitution models.
 */
class SubstitutionModelFactory
{
  public:
    static const string JUKES_CANTOR;
    static const string KIMURA_2P;
    static const string HASEGAWA_KISHINO_YANO;
    static const string TAMURA_NEI;
    static const string GENERAL_TIME_REVERSIBLE;
    static const string TAMURA;
    static const string FELSENSTEIN84;
    static const string JOHN_TAYLOR_THORNTON;
    static const string DAYHOFF_SCHWARTZ_ORCUTT;
  
  protected:
    const Alphabet * _alphabet;
  
  public:
    /**
     * @brief Creates a new factory object with the given alphabet.
     *
     * @param alphabet The alphabet for wich models must be instanciated.
     *
     * Example:
     * @code
     * const Alphabet * alphabet = new DNA();
     * SubstitutionModel * model = SubstitutionModelFactory(alphabet)
     *     .createModel(SubstitutionModelFactory::TAMURA);
     * // model can be used in any object dealing with a nucleotide substitution models.
     * @endcode
     */
    SubstitutionModelFactory(const Alphabet * alphabet): _alphabet(alphabet) {}
    virtual ~SubstitutionModelFactory() {}

  public:
    /**
     * @brief Get a new dynamically created SubstitutionModel object.
     *
     * @param modelName The name of the model to use.
     * @return A pointer toward a new substitution model, with default parameter values.
     * @throw AlphabetException If the model is not compatible with the given alphabet.
     * @throw Exception If the model name do not match any available model.
     */
    virtual SubstitutionModel * createModel(const string& modelName) const throw (AlphabetException, Exception);

};

} //end of namespace bpp.

#endif //_SUBSTITUTIOMODELFACTORY_H_

