// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include "Combinatorics.hpp"
#include <cmath>

using namespace SphericalFunctions;
using std::vector;

const FactorialSingleton* FactorialSingleton::FactorialInstance = NULL;
const BinomialCoefficientSingleton* BinomialCoefficientSingleton::BinomialCoefficientInstance = NULL;
const LadderOperatorFactorSingleton* LadderOperatorFactorSingleton::LadderOperatorFactorInstance = NULL;
