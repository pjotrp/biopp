//
// File: NonHomogeneousSequenceSimulator.h
// Created by: Julien Dutheil
//             Bastien Boussau
// Created on: Wed Aug  24 15:20 2005
//

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

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

#ifndef _NONHOMOGENEOUSSEQUENCESIMULATOR_H_
#define _NONHOMOGENEOUSSEQUENCESIMULATOR_H_

#include "DetailedSiteSimulator.h"
#include "SequenceSimulator.h"
#include "TreeTemplate.h"
#include "NodeTemplate.h"
#include "SubstitutionModel.h"

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/Site.h>
#include <Seq/SiteContainer.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/RandomTools.h>

// From the STL:
#include <map>
#include <vector>

using namespace std;

#include "SubstitutionModelSet.h"

namespace bpp
{

class SimData
{
  public:
    int state;
    vector<int> states;
    VVVdouble cumpxy;
    const SubstitutionModel* model;
};

typedef NodeTemplate<SimData> SNode;

/**
 * @brief Site and sequences simulation under non-homogeneous models.
 *
 * Rate across sites variation is supported, using a DiscreteDistribution object or by specifying explicitely the rate of the sites to simulate.
 */
class NonHomogeneousSequenceSimulator:
  public DetailedSiteSimulator,
  public SequenceSimulator
{
	protected:
		const SubstitutionModelSet *_modelSet;
		const Alphabet             *_alphabet;
		const DiscreteDistribution *_rate;
		const TreeTemplate<Node>* _templateTree;
		mutable TreeTemplate<SNode> _tree;
    bool _ownModelSet;
	
		/**
		 * @brief This stores once for all all leaves in a given order.
		 * This order will be used during site creation.
		 */
		vector<SNode *> _leaves;
	
		vector<string> _seqNames;

		unsigned int _nbNodes;
		unsigned int _nbClasses;
		unsigned int _nbStates;

    bool _continuousRates;
	
		/**
		 * @name Stores intermediate results.
     *
     * @{
		 */

	public:		
		NonHomogeneousSequenceSimulator(
			const SubstitutionModelSet * modelSet,
			const DiscreteDistribution * rate,
			const TreeTemplate<Node> * tree
		) throw (Exception);

    NonHomogeneousSequenceSimulator(
			const SubstitutionModel * model,
			const DiscreteDistribution * rate,
			const TreeTemplate<Node> * tree
		);
			
		virtual ~NonHomogeneousSequenceSimulator()
    {
      if(_ownModelSet && _modelSet != NULL) delete _modelSet;
    }

#ifndef NO_VIRTUAL_COV
    NonHomogeneousSequenceSimulator *
#else
    Clonable *
#endif
    clone() const { return new NonHomogeneousSequenceSimulator(*this); }

  private:
    /**
     * @brief Init all probabilities.
     *
     * Method called by constructors.
     */
    void init();

	public:
	
		/**
		 * @name The SiteSimulator interface
		 *
		 * @{
		 */
		Site * simulate() const;
		Site * simulate(int ancestralState) const;
		Site * simulate(int ancestralState, double rate) const;
		Site * simulate(double rate) const;
		vector<string> getSequencesNames() const { return _seqNames; }
		/** @} */
    
		/**
		 * @name The PreciseSiteSimulator interface.
		 *
		 * @{
		 */
    RASiteSimulationResult * dSimulate() const;
    
    RASiteSimulationResult * dSimulate(int ancestralState) const;
    
    RASiteSimulationResult * dSimulate(int ancestralState, double rate) const;

    RASiteSimulationResult * dSimulate(double rate) const;
		/** @} */

    /**
		 * @name The SequenceSimulator interface
		 *
		 * @{
		 */
		SiteContainer * simulate(unsigned int numberOfSites) const;
		/** @} */
    
		/**
		 * @name SiteSimulator and SequenceSimulator interface
		 *
		 * @{
		 */
		const Alphabet * getAlphabet() const { return _alphabet; }
		/** @} */

    /**
     * @name Functions with rate classes instead of absolute rates.
     *
     * @{
     */
		virtual Site * simulate(int ancestralState, unsigned int rateClass) const;
    virtual	RASiteSimulationResult * dSimulate(int ancestralState, unsigned int rateClass) const;
    /** @} */
	
		/**
		 * @brief Get the mutation process associated to this instance.
		 *
		 * @return The MutationProcess object associated to this instance.
		 */
		const SubstitutionModelSet * getSubstitutionModelSet() const { return _modelSet; }
		
	
		
		/**
		 * @brief Get the rate distribution associated to this instance.
		 *
		 * @return The DiscreteDistribution object associated to this instance.
		 */
		const DiscreteDistribution * getRateDistribution() const { return _rate; }

		/**
		 * @brief Get the tree associated to this instance.
		 *
		 * @return The Tree object associated to this instance.
		 */
		const TreeTemplate<Node> * getTree() const { return _templateTree; }

    /**
     * @brief Enable the use of continuous rates instead of discrete rates.
     *
     * To work, the DiscreteDistribution object used should implement the randC method.
     *
     * @param yn Tell if we should use continuous rates.
     */
    void enableContinuousRates(bool yn) { _continuousRates = yn; }
	
	protected:
		
		/**
		 * @brief Evolve from an initial state along a branch, knowing the evolutionary rate class.
		 *
		 * This method is fast since all pijt have been computed in the constructor of the class.
     * This method is used for the implementation of the SiteSimulator interface.
		 */
		int evolve(const SNode * node, int initialState, unsigned int rateClass) const;
		
		/**
		 * @brief Evolve from an initial state along a branch, knowing the evolutionary rate.
		 *
		 * This method is slower than the previous one since exponential terms must be computed.
     * This method is used for the implementation of the SiteSimulator interface.
		 */
		int evolve(const SNode * node, int initialState, double rate) const;
		
    /**
     * @brief The same as the evolve(initialState, rateClass) function, but for several sites at a time.
     *
     * This method is used for the implementation of the SequenceSimulator interface.
     */
		void multipleEvolve(const SNode * node, const Vint & initialState, const vector<unsigned int> & rateClasses, Vint & finalStates) const;
		SiteContainer * multipleEvolve(const Vint & initialStates, const vector<unsigned int> & rateClasses) const;
		
    void dEvolve(int initialState, double rate, RASiteSimulationResult & rassr) const;
		
    /**
     * @name The 'Internal' methods.
     *
     * @{
     */

    /**
     * This method uses the _states variable for saving ancestral states.
     */
    void evolveInternal(SNode * node, unsigned int rateClass) const;
    /**
     * This method uses the _states variable for saving ancestral states.
     */
		void evolveInternal(SNode * node, double rate) const;
    /**
     * This method uses the _multipleStates variable for saving ancestral states.
     */
 		void multipleEvolveInternal(SNode * node, const vector<unsigned int> & rateClasses) const;

    /**
     * This method uses the _states variable for saving ancestral states.
     */
		void dEvolveInternal(SNode * node, double rate, RASiteSimulationResult & rassr) const;
    /** @} */

};

} //end of namespace bpp.

#endif //_NONHOMOGENEOUSSEQUENCESIMULATOR_H_

