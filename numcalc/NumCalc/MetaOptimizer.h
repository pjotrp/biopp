//
// File: MetaOptimizer.h
// Created by: Julien Dutheil
// Created on: Fri Oct 12 16:05 2007
// From file: NewtonBrentMetaOptimizer.h
// Created on: Tue Nov 17 17:22 2004
// 
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#ifndef _METAOPTIMIZER_H_
#define _METAOPTIMIZER_H_

#include "AbstractOptimizer.h"

// From the STL:
#include <vector>
using namespace std;

namespace bpp
{

/**
 * @brief Provide a list of optimizer and corresponding options to be used with the MetaOptimizer class.
 */
class MetaOptimizerInfos:
  public virtual Clonable
{
  public:
    static string IT_TYPE_STEP;
    static string IT_TYPE_FULL;

  protected:
    vector<string> _names;
    vector<Optimizer *> _optimizers;
    vector< vector<string> > _parameterNames;
    vector<unsigned short> _derivatives;
    vector<string> _itTypes;

  public:
    MetaOptimizerInfos() {}
    MetaOptimizerInfos(const MetaOptimizerInfos & infos)
    {
      _optimizers     = infos._optimizers;
      _parameterNames = infos._parameterNames;
      _derivatives    = infos._derivatives;
      _itTypes        = infos._itTypes;
      for(unsigned int i = 0; i < _optimizers.size(); i++)
        _optimizers[i] = dynamic_cast<Optimizer *>(infos._optimizers[i]->clone());
    }

    MetaOptimizerInfos& operator=(const MetaOptimizerInfos & infos)
    {
      _optimizers     = infos._optimizers;
      _parameterNames = infos._parameterNames;
      _derivatives    = infos._derivatives;
      _itTypes        = infos._itTypes;
      for(unsigned int i = 0; i < _optimizers.size(); i++)
        _optimizers[i] = dynamic_cast<Optimizer *>(infos._optimizers[i]->clone());
      return *this;
    }

    virtual ~MetaOptimizerInfos()
    {
      for(unsigned int i = 0; i < _optimizers.size(); i++)
        delete _optimizers[i];
    }

  public:
#ifndef NO_VIRTUAL_COV
    MetaOptimizerInfos * clone() const { return new MetaOptimizerInfos(*this); }
#else
    Clonable * clone() const { return new MetaOptimizerInfos(*this); }
#endif

  public:
    /**
     * @brief Add a new optimizer to the set.
     *
     * @param name the name of the optimizer. It is used for display only.
     * @param optimizer A pointer toward the optimizer to add. The set will own the underlying object, which will be destroyed together with the set.
     * @param params A list of parameter names to optimize with this optimizer.
     * @param derivatives 0, 1 or 2: does this parameter use no, first order or second order derivatives?
     * @param type For each optimization step, shall we perform a full optimization with this optimizer or only one step?
     */
    virtual void addOptimizer(const string & name, Optimizer * optimizer, const vector<string> & params, const short derivatives = 0, const string & type = IT_TYPE_STEP)
    {
      _names.push_back(name);
      _optimizers.push_back(optimizer);
      _parameterNames.push_back(params);
      _derivatives.push_back(derivatives);
      _itTypes.push_back(type);
    }

    /**
     * @return The display name of the ith optimizer in the set.
     */
    virtual const string & getName(unsigned int i) const { return _names[i]; }
    
    /**
     * @return The ith optimizer in the set.
     */
    virtual Optimizer * getOptimizer(unsigned int i) { return _optimizers[i]; }
    /**
     * @return The ith optimizer in the set.
     */
    virtual const Optimizer * getOptimizer(unsigned int i) const { return _optimizers[i]; }

    /**
     * @return The parameter names associated to the ith optimizer in the set.
     */
    virtual vector<string> & getParameterNames(unsigned int i) { return _parameterNames[i]; }
    /**
     * @return The parameter names associated to the ith optimizer in the set.
     */
    virtual const vector<string> & getParameterNames(unsigned int i) const { return _parameterNames[i]; }

    /**
     * @return The type of iteration to perform for the ith optimizer in the set.
     */
    virtual string & getIterationType(unsigned int i) { return _itTypes[i]; }
    /**
     * @return The type of iteration to perform for the ith optimizer in the set.
     */
    virtual const string & getIterationType(unsigned int i) const { return _itTypes[i]; }

    /**
     * @return True if the ith optimizer in the set requires first order derivatives.
     */
    virtual bool requiresFirstOrderDerivatives(unsigned int i) const { return _derivatives[i] > 0; }  
    /**
     * @return True if the ith optimizer in the set requires second order derivatives.
     */
    virtual bool requiresSecondOrderDerivatives(unsigned int i) const { return _derivatives[i] > 1; }  

    /**
     * @return The number of optimizers in the set.
     */
    virtual unsigned int getNumberOfOptimizers() const { return _optimizers.size(); }
};

/**
 * @brief Meta-optimizer.
 *
 * This optimizer uses a set of optimizers to applyied sequentially on distinct parameters.
 * The set of optimizers is fully specified by a MetaOptimizerInfos object.
 * 
 * To decrease the optimization time, the precision of the optimizers can be increased progressively:
 * if @f$\varepsilon@f$ is the final precision required, one may consider using a precision increment of @f$\sigma=\log_10(\varepsilon/n)@f$, where @f$n@f$ is the number of progressive steps.
 * During the first step optimization step, the precisions of type 1 and 2 optimizers are set to @f$10^{\sigma}@f$, @f$10^{2\sigma}@f$ for step 2, ... until precision @f$10^{n\sigma}=\varepsilon@f$ at step @f$n@f$ and later.
 * This saves some time spending in the first steps of the estimation.
 * The number of steps @f$n@f$ is set in the constructor of the optimizer.
 *
 * This optimizer can be used with numerical derivatives.
 * 
 * @see MetaOptimizerInfos.
 */
class MetaOptimizer:
  public AbstractOptimizer
{
	protected:
    MetaOptimizerInfos * _optDesc;
    vector<ParameterList> _optParameters;
    vector<unsigned int> _nbParameters;
    unsigned int _n;
    double _precisionStep;
    unsigned int _stepCount;
    double _initialValue;
		
	public:
    /**
     * @brief Build a new MetaOptimizer object.
     *
     * @param function The function to be optimized.
     * @param desc     A MetaOptimizerInfos object describing the optimizers to use.
     *                 The optimizer will own the instance of the MetaOptimizerInfos object.
     * @param n        The number of progressive steps to use in optimization).
     */
		MetaOptimizer(Function * function, MetaOptimizerInfos * desc, unsigned int n = 1);

		virtual ~MetaOptimizer();

    MetaOptimizer(const MetaOptimizer & opt);

    MetaOptimizer & operator=(const MetaOptimizer & opt);

#ifndef NO_VIRTUAL_COV
    MetaOptimizer*
#else
    Clonable*
#endif
    clone() const { return new MetaOptimizer(*this); }

	public:
		
    void setFunction(Function * function)
    { 
      AbstractOptimizer::setFunction(function);
      for(unsigned int i = 0; i < _optDesc->getNumberOfOptimizers(); i++)
        _optDesc->getOptimizer(i)->setFunction(function);
    }

		void doInit(const ParameterList & parameters) throw (Exception);
    
    double doStep() throw (Exception);

    /**
     * @return The MetaOptimizerInfos object associated to this optimizer.
     */
    MetaOptimizerInfos * getOptimizers() { return _optDesc; }
    
    /**
     * @return The MetaOptimizerInfos object associated to this optimizer.
     */
    const MetaOptimizerInfos * getOptimizers() const { return _optDesc; }

};

} //end of namespace bpp.

#endif //_METAOPTIMIZER_H_

