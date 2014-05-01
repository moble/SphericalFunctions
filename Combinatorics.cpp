// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include "Combinatorics.hpp"
#include <cmath>

using namespace SphericalFunctions;
using std::vector;

const FactorialSingleton* FactorialSingleton::FactorialInstance = 0;
const BinomialCoefficientSingleton* BinomialCoefficientSingleton::BinomialCoefficientInstance = 0;
const LadderOperatorFactorSingleton* LadderOperatorFactorSingleton::LadderOperatorFactorInstance = 0;
