//
// File: AbstractOptimizer.h
// Created by: Julien Dutheil
// Created on: Mon Dec 22 12:18:09 2003
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

#ifndef _ABSTRACTOPTIMIZER_H_
#define _ABSTRACTOPTIMIZER_H_

#include "Optimizer.h"

namespace bpp
{

/**
 * @brief Partial implementation of the Optimizer interface.
 *
 * This implementation is designed for unconstrained or simple-bounded optimization.
 * You should not use it with global contraints.
 */
class AbstractOptimizer:
  public Optimizer
{
  protected:
    
    /**i
     * @brief The function to optimize.
     */
    Function * _function;
  
    /**
     * @brief The parameters that will be optimized.
     */
    ParameterList _parameters;
  
    /**
     * @brief The message handler.
     */
    ostream * _messageHandler;
  
    /**
     * @brief The profiler.
     */
    ostream * _profiler;
    
    /**
     * @brief The constraint policy.
     *
     * Must be one the following:
     * - CONSTRAINTS_KEEP: keep the constraint associated to the parameters (default).
     * - CONSTRAINTS_IGNORE: remove all constraints.
     * - CONSTRAINTS_AUTO: use AutoParameters to deal with constraints.
     *
     * @see AutoParameter
     */      
    string _constraintPolicy;

    /**
     * @brief Tell if the tolerance level has been reached.
     *
     * This field is initilaised by the init() method, maintained by the
     * step() method and used in the optimize() method.
     */
    bool _tolIsReached;
    
    /**
     * @brief The stoping condition to use while optimizing.
     */
    OptimizationStopCondition * _stopCondition;
    
    /**
     * @brief The default stoping condition to use while optimizing.
     */
    OptimizationStopCondition * _defaultStopCondition;

    /**
     * @brief The maximum number of function evaluations allowed.
     */
    unsigned int _nbEvalMax;
    
    /**
     * @brief The current number of function evaluations achieved.
     */
    unsigned int _nbEval;

    /**
     * @brief State of the verbose mode: > 0 = enabled.
     *
     * This may not be used by the Optimizer.
     */
    unsigned int _verbose;

    /**
     * @brief The current value of the function.
     */
    double _currentValue;

    /**
     * @brief Check if the optimizer have been feeded with initial parameters values.
     */
    bool _isInitialized;

    time_t _startTime;

    vector<OptimizationListener *> _listeners;

    bool _updateParameters;

    string _stepChar;

  public:
    AbstractOptimizer(Function * function = NULL);

    AbstractOptimizer(const AbstractOptimizer & opt);
    
    AbstractOptimizer & operator=(const AbstractOptimizer & opt);

    virtual ~AbstractOptimizer()
    {
      delete _stopCondition;
      delete _defaultStopCondition;
    }
  
  public:
    
    /**
     * @name The Optimizer interface.
     *
     * @{
     */
    /**
     * @brief Basic implementation.
     *
     * Store all parameters, call the doInit method, print to profiler, initialize timer and notify all listeners.
     */
    void init(const ParameterList & params) throw (Exception);
    /**
     * @brief Basic implementation.
     *
     * Check if the optimizer is initialized, check if parameters need update because of listeners, call the doStep method, print the current point to the profiler, notify all listeners and return the current value of the function.
     */
    double step() throw (Exception);
    /**
     * @brief Basic implementation.
     *
     * Call the step method untill tolerance is reached.
     */
    double optimize() throw (Exception);
    bool isInitialized() const { return _isInitialized; }
    ParameterList getParameters() const { return _parameters; }
    void setFunction(Function * function)
    { 
      _function = function;
      if(function != NULL) _stopCondition->init();
    }
    const Function * getFunction() const { return _function; }
    Function * getFunction() { return _function; }
    double getFunctionValue() const throw (NullPointerException)
    {
      if(_function == NULL) throw NullPointerException("AbstractOptimizer::getFunctionValue. No function associated to this optimizer.");
      return _currentValue;
    }
    void setMessageHandler(ostream * mh) { _messageHandler = mh; }
    void setProfiler(ostream * profiler) { _profiler = profiler; }
    int getNumberOfEvaluations() const { return _nbEval; }
    void setStopCondition(const OptimizationStopCondition & stopCondition)
    {
      _stopCondition = dynamic_cast<OptimizationStopCondition *>(stopCondition.clone());
    }
    OptimizationStopCondition * getStopCondition() { return _stopCondition; }
    const OptimizationStopCondition * getStopCondition() const { return _stopCondition; }
    OptimizationStopCondition * getDefaultStopCondition() { return _defaultStopCondition; }
    const OptimizationStopCondition * getDefaultStopCondition() const { return _defaultStopCondition; }
    bool isToleranceReached() const { return _tolIsReached; }
    bool isMaximumNumberOfEvaluationsReached() const { return _nbEvalMax >= _nbEvalMax; }
    void setMaximumNumberOfEvaluations(unsigned int max) { _nbEvalMax = max; }
    void setVerbose(unsigned int v) { _verbose = v; }
    unsigned int getVerbose() const { return _verbose; }
    void setConstraintPolicy(const string & constraintPolicy) { _constraintPolicy = constraintPolicy; }
    string getConstraintPolicy() const { return _constraintPolicy; }
    void addOptimizationListener(OptimizationListener * listener) { if(listener) _listeners.push_back(listener); }
    /** @} */

    /**
     * @brief Tell if we shall update all parameters after one optimization step.
     *
     * This is required only for functions that have non-independent parameters,
     * which means that setting one parameter value may modify one or several other parameters.
     * Depending on the optimizer, this may have no effect.
     *
     * @param yn true/false
     */
    void updateParameters(bool yn) { _updateParameters = yn; }

    /**
     * @brief Tell if we shall update all parameters after one optimization step.
     *
     * This is required only for functions that have non-independent parameters,
     * which means that setting one parameter value may modify one or several other parameters.
     * Depending on the optimizer, this may have no effect.
     *
     * @return yn true/false
     */
    bool updateParameters() const { return _updateParameters; }

    /**
     * @brief Set the character to be displayed during optimization.
     *
     * @param c A character.
     */
    void setOptimizationProgressCharacter(const string & c) { _stepChar = c; }
    /**
     * @return The character to be displayed during optimization.
     */
    string getOptimizationProgressCharacter() const { return _stepChar; }
  
  protected:

    /**
     * @brief This function is called by the init() method and contains all calculations.
     *
     * @param params The parameters to use for initialization.
     */
    virtual void doInit(const ParameterList & params) throw (Exception) = 0;
    
    /**
     * @brief This function is called by the step() method and contains all calculations.
     *
     * @return The value of the function after the optimization step.
     */
    virtual double doStep() throw (Exception) = 0;
    
    /**
     * @name Inner utilitary functions
     *
     * @{
     */
    
    /**
     * @brief Build a list of AutoParameter instead of Parameter.
     */
    void autoParameter();
  
    /**
     * @brief Remove the constraints of all the arguments.
     */
    void ignoreConstraints();
  
    /**
     * @brief Print to the profile if there is one.
     *
     * @param v The double value to print.
     */
    void profile(double v);
  
    /**
     * @brief Print to the profile if there is one.
     *
     * @param s The string to print to the profile.
     */
    void profile(const string & s);
  
    /**
     * @brief Print to the profile if there is one and end line.
     *
     * @param v The double value to print.
     */
    void profileln(double v);
  
    /**
     * @brief Print to the profile if there is one and end line.
     *
     * @param s The string to print to the profile.
     */
    void profileln(const string & s);
  
    /**
     * @brief Print parameters and corresponding function evaluation to profiler.
     *
     * @param params The parameters to print.
     * @param value  The function evaluation.
     */
    void printPoint(const ParameterList & params, double value);
    
    /**
     * @brief Give a message to print to the message handler.
     *
     * @param message The message to print.
     */
    void printMessage(const string & message);

    /**
     * @brief Notify all listeners that optimizer initialization was performed.
     *
     * This method should be called by the init method.
     *
     * @param event An OptimizationEvent object.
     */
    void fireOptimizationInitializationPerformed(const OptimizationEvent & event);

    /**
     * @brief Notify all listeners that an optimization step was performed.
     *
     * This method should be called by the step method.
     *
     * @param event An OptimizationEvent object.
     */
    void fireOptimizationStepPerformed(const OptimizationEvent & event);

    bool listenerModifiesParameters() const;
    /** @} */
  
};

} //end of namespace bpp.

#endif  //_ABSTRACTOPTIMIZER_H_

