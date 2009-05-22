#Copyright or © or Copr. CNRS, (November 16, 2004)
#
#This software is a computer program whose purpose is to provide classes
#for numerical calculus. This file is part of the Bio++ project.
#
#This software is governed by the CeCILL  license under French law and
#abiding by the rules of distribution of free software.  You can  use, 
#modify and/ or redistribute the software under the terms of the CeCILL
#license as circulated by CEA, CNRS and INRIA at the following URL
#"http://www.cecill.info". 
#
#As a counterpart to the access to the source code and  rights to copy,
#modify and redistribute granted by the license, users are provided only
#with a limited warranty  and the software's author,  the holder of the
#economic rights,  and the successive licensors  have only  limited
#liability. 
#
#In this respect, the user's attention is drawn to the risks associated
#with loading,  using,  modifying and/or developing or reproducing the
#software by the user in light of its specific status of free software,
#that may mean  that it is complicated to manipulate,  and  that  also
#therefore means  that it is reserved for developers  and  experienced
#professionals having in-depth computer knowledge. Users are therefore
#encouraged to load and test the software's suitability as regards their
#requirements in conditions enabling the security of their systems and/or 
#data to be ensured and,  more generally, to use and operate it in the 
#same conditions as regards security. 
#
#The fact that you are presently reading this means that you have had
#knowledge of the CeCILL license and that you accept its terms.

H_FILES = \
AbstractOptimizer.h\
DownhillSimplexMethod.h\
NumTools.h\
Parametrizable.h\
AbstractParametrizable.h\
ParameterAliasable.h\
AbstractParameterAliasable.h\
Functions.h\
FunctionTools.h\
OneDimensionOptimizationTools.h\
PowellMultiDimensions.h\
BrentOneDimension.h\
GammaDiscreteDistribution.h\
ExponentialDiscreteDistribution.h\
OptimizationStopCondition.h\
RandomTools.h\
ConstantDistribution.h\
GoldenSectionSearch.h\
Optimizer.h\
VectorTools.h\
VectorExceptions.h\
Constraints.h\
IntervalData.h\
ParameterExceptions.h\
DiscreteDistribution.h\
AbstractDiscreteDistribution.h\
SimpleDiscreteDistribution.h\
SimpleMultiDimensions.h\
Matrix.h\
MatrixTools.h\
EigenValue.h\
LUDecomposition.h\
Parameter.h\
AutoParameter.h\
TransformedParameter.h\
Domain.h\
ParameterList.h\
RandomFactory.h\
Uniform01K.h\
Uniform01QD.h\
Uniform01WH.h\
DataTable.h\
DataTableExceptions.h\
parameters\
optimizers\
vectors\
matrix\
distributions\
random\
AbstractNumericalDerivative.h\
TwoPointsNumericalDerivative.h\
ThreePointsNumericalDerivative.h\
FivePointsNumericalDerivative.h\
DirectionFunction.h\
NewtonOneDimension.h\
SimpleNewtonMultiDimensions.h\
ConjugateGradientMultiDimensions.h\
MetaOptimizer.h\
InvariantMixedDiscreteDistribution.h\
StatTest.h\
NumCalcApplicationTools.h\
functions\
ReparametrizationFunctionWrapper.h\
NumConstants.h\
TestUnit.h


CPP_FILES = \
AbstractOptimizer.cpp\
AutoParameter.cpp\
GammaDiscreteDistribution.cpp\
ExponentialDiscreteDistribution.cpp\
Parameter.cpp\
BrentOneDimension.cpp\
GoldenSectionSearch.cpp\
ParameterExceptions.cpp\
ConstantDistribution.cpp\
IntervalData.cpp\
ParameterList.cpp\
AbstractDiscreteDistribution.cpp\
SimpleDiscreteDistribution.cpp\
SimpleMultiDimensions.cpp\
PowellMultiDimensions.cpp\
Domain.cpp\
OneDimensionOptimizationTools.cpp\
RandomTools.cpp\
DownhillSimplexMethod.cpp\
OptimizationStopCondition.cpp\
VectorTools.cpp\
Uniform01K.cpp\
Uniform01QD.cpp\
Uniform01WH.cpp\
DataTable.cpp\
TwoPointsNumericalDerivative.cpp\
ThreePointsNumericalDerivative.cpp\
FivePointsNumericalDerivative.cpp\
DirectionFunction.cpp\
NewtonOneDimension.cpp\
SimpleNewtonMultiDimensions.cpp\
ConjugateGradientMultiDimensions.cpp\
MetaOptimizer.cpp\
NumTools.cpp\
InvariantMixedDiscreteDistribution.cpp\
ReparametrizationFunctionWrapper.cpp\
NumConstants.cpp\
AbstractParametrizable.cpp\
AbstractParameterAliasable.cpp\
NumCalcApplicationTools.cpp\
FunctionTools.cpp\
TestUnit.cpp


lib_LTLIBRARIES = libbpp-numcalc.la
libbpp_numcalc_la_SOURCES = $(H_FILES) $(CPP_FILES)
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current: The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age: The difference between the newest and oldest interfaces that this library implements.
# In other words, the library implements all the interface numbers in the range from number current - age to current.
libbpp_numcalc_la_LDFLAGS = -version-info 7:0:2

library_includedir=$(includedir)/NumCalc
library_include_HEADERS = $(H_FILES)

