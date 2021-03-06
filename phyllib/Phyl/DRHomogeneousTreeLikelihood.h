//
// File: DRHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 18:14:51 2003
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

#ifndef _DRHOMOGENEOUSTREELIKELIHOOD_H_
#define _DRHOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractHomogeneousTreeLikelihood.h"
#include "DRTreeLikelihood.h"
#include "DRASDRTreeLikelihoodData.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/DiscreteDistribution.h>

namespace bpp
{

/**
 * @brief This class implements the likelihood computation for a tree using the double-recursive
 * algorithm.
 *
 * The substitution model is the same over the tree (homogeneous model).
 * A non-uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * This class uses an instance of the DRASDRTreeLikelihoodData for conditionnal likelihood storage.
 *
 * All nodes share the same site patterns.
 */
class DRHomogeneousTreeLikelihood:
  public AbstractHomogeneousTreeLikelihood,
  public DRTreeLikelihood
{
  protected:
    mutable DRASDRTreeLikelihoodData *_likelihoodData;
    
  public:
    /**
     * @brief Build a new DRHomogeneousTreeLikelihood object without data.
     *
     * This constructor only initialize the parameters.
     * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
     *
     * @param tree The tree to use.
     * @param model The substitution model to use.
     * @param rDist The rate across sites distribution to use.
     * @param checkRooted Tell if we have to check for the tree to be unrooted.
     * If true, any rooted tree will be unrooted before likelihood computation.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
    DRHomogeneousTreeLikelihood(
      const Tree & tree,
      SubstitutionModel* model,
      DiscreteDistribution* rDist,
      bool checkRooted = true,
      bool verbose = true)
      throw (Exception);
  
    /**
     * @brief Build a new DRHomogeneousTreeLikelihood object and compute the corresponding likelihood.
     *
     * This constructor initializes all parameters, data, and likelihood arrays.
     *
     * @param tree The tree to use.
     * @param data Sequences to use.
     * @param model The substitution model to use.
     * @param rDist The rate across sites distribution to use.
     * @param checkRooted Tell if we have to check for the tree to be unrooted.
     * If true, any rooted tree will be unrooted before likelihood computation.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
    DRHomogeneousTreeLikelihood(
      const Tree & tree,
      const SiteContainer & data,
      SubstitutionModel* model,
      DiscreteDistribution* rDist,
      bool checkRooted = true,
      bool verbose = true)
      throw (Exception);

    /**
     * @brief Copy constructor.
     */ 
    DRHomogeneousTreeLikelihood(const DRHomogeneousTreeLikelihood & lik);
    
    DRHomogeneousTreeLikelihood & operator=(const DRHomogeneousTreeLikelihood & lik);

    virtual ~DRHomogeneousTreeLikelihood();

#ifndef NO_VIRTUAL_COV
    DRHomogeneousTreeLikelihood*
#else
    Clonable*
#endif
    clone() const { return new DRHomogeneousTreeLikelihood(*this); }

  private:

    /**
     * @brief Method called by constructors.
     */
    void _init() throw (Exception);

  public:

    /**
     * @name The TreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractTreeLikelihood class.
     *
     * @{
     */
    void setData(const SiteContainer & sites) throw (Exception);
    double getLikelihood () const;
    double getLogLikelihood() const;
    double getLikelihoodForASite (unsigned int site) const;
    double getLogLikelihoodForASite(unsigned int site) const;
    /** @} */

    void computeTreeLikelihood();

    
    /**
     * @name The DiscreteRatesAcrossSites interface implementation:
     *
     * @{
     */
    double getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
    double getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
    double getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
    double getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
    /** @} */
  
    /**
     * @brief Implements the Function interface.
     *
     * Update the parameter list and call the fireParameterChanged() method.
     *
     * If a subset of the whole parameter list is passed to the function,
     * only these parameters are updated and the other remain constant (i.e.
     * equal to their last value).
     *
     * @param parameters The parameter list to pass to the function.
     */
    void setParameters(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException);
    
    /**
     * @brief Function and NNISearchable interface.
     */
    double getValue() const throw (Exception);
    
    /**
     * @name DerivableFirstOrder interface.
     *
     * @{
     */
    double getFirstOrderDerivative(const string & variable) const throw (Exception);
    /** @{ */

    /**
     * @name DerivableSecondOrder interface.
     *
     * @{
     */
    double getSecondOrderDerivative(const string & variable) const throw (Exception);
    double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) { return 0; } // Not implemented for now.
    /** @} */
    
  public:  // Specific methods:

    DRASDRTreeLikelihoodData * getLikelihoodData() { return _likelihoodData; }
    const DRASDRTreeLikelihoodData * getLikelihoodData() const { return _likelihoodData; }
  
    virtual void computeLikelihoodAtNode(int nodeId, VVVdouble& likelihoodArray) const;

    /**
     * @brief Retrieves all Pij(t) for a particular node.
     *
     * These intermediate results may be used by other methods.
     */
    virtual const VVVdouble & getTransitionProbabilitiesForNode(const Node * node) const { return pxy_[node->getId()]; }
       
  protected:
  
    /**
     * Initialize the arrays corresponding to each son node for the node passed as argument.
     * The method is called for each son node and the result stored in the corresponding array.
     */
    virtual void computeSubtreeLikelihoodPostfix(const Node * node); //Recursive method.
    /**
     * This method initilize the remaining likelihood arrays, corresponding to father nodes.
     * It must be called after the postfix method because it requires that the arrays for
     * son nodes to be be computed.
     */
    virtual void computeSubtreeLikelihoodPrefix(const Node * node); //Recursive method.

    virtual void computeRootLikelihood();

    virtual void computeTreeDLikelihoodAtNode(const Node * node);
    virtual void computeTreeDLikelihoods();
    
    virtual void computeTreeD2LikelihoodAtNode(const Node * node);
    virtual void computeTreeD2Likelihoods();

    void fireParameterChanged(const ParameterList & params);

    void resetLikelihoodArrays(const Node * node);
  
    /**
     * @brief This method is mainly for debugging purpose.
     *
     * @param node The node at which likelihood values must be displayed.
     */
    virtual void displayLikelihood(const Node * node);

    /**
     * @brief Compute conditional likelihoods.
     *
     * This method is the "core" likelihood computation function, performing all the product uppon all nodes, the summation for each ancestral state and each rate class.
     * It is designed for inner usage, and a maximum efficiency, so no checking is performed on the input parameters.
     * Use with care!
     * 
     * @param iLik A vector of likelihood arrays, one for each conditional node.
     * @param tProb A vector of transition probabilities, one for each node.
     * @param oLik The likelihood array to store the computed likelihoods.
     * @param nbNodes The number of nodes = the size of the input vectors.
     * @param nbDistinctSites The number of distinct sites (the first dimension of the likelihood array).
     * @param nbClasses The number of rate classes (the second dimension of the likelihood array).
     * @param nbStates The number of states (the third dimension of the likelihood array).
     * @param reset Tell if the output likelihood array must be initalized prior to computation.
     * If true, the resetLikelihoodArray method will be called.
     */
    static void computeLikelihoodFromArrays(const vector<const VVVdouble *> & iLik, const vector<const VVVdouble *> & tProb, VVVdouble & oLik, unsigned int nbNodes, unsigned int nbDistinctSites, unsigned int nbClasses, unsigned int nbStates, bool reset = true);

    /**
     * @brief Compute conditional likelihoods.
     *
     * This method is the "core" likelihood computation function, performing all the product uppon all nodes, the summation for each ancestral state and each rate class.
     * This function is specific to non-reversible models: the subtree containing the root is specified separately.
     * It is designed for inner usage, and a maximum efficiency, so no checking is performed on the input parameters.
     * Use with care!
     * 
     * @param iLik A vector of likelihood arrays, one for each conditional node.
     * @param tProb A vector of transition probabilities, one for each node.
     * @param iLikR The likelihood array for the subtree containing the root of the tree.
     * @param tProbR The transition probabilities for thr subtree containing the root of the tree.
     * @param oLik The likelihood array to store the computed likelihoods.
     * @param nbNodes The number of nodes = the size of the input vectors.
     * @param nbDistinctSites The number of distinct sites (the first dimension of the likelihood array).
     * @param nbClasses The number of rate classes (the second dimension of the likelihood array).
     * @param nbStates The number of states (the third dimension of the likelihood array).
     * @param reset Tell if the output likelihood array must be initalized prior to computation.
     * If true, the resetLikelihoodArray method will be called.
     */
    static void computeLikelihoodFromArrays(
        const vector<const VVVdouble *> & iLik,
        const vector<const VVVdouble *> & tProb,
        const VVVdouble * iLikR,
        const VVVdouble * tProbR,
        VVVdouble & oLik,
        unsigned int nbNodes,
        unsigned int nbDistinctSites,
        unsigned int nbClasses,
        unsigned int nbStates,
        bool reset = true);
};

} //end of namespace bpp.

#endif  //_DRHOMOGENEOUSTREELIKELIHOOD_H_

