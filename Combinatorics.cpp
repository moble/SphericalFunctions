// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#include "Combinatorics.hpp"
#include <cmath>

using namespace SphericalFunctions;
using std::vector;

/// Class to create an object returning the factorial of an argument.
vector<double> FactorialTableCalculator() {
  /// Note that because a double is returned, only values up to 28!
  /// will be exact; higher values will be accurate to machine
  /// precision.  Values up to 170! only are allowed because higher
  /// values overflow.
  vector<double> FactorialTable(171);
  FactorialTable[0] = 1.0;
  for (int i=1;i<171;i++) {
    FactorialTable[i] = i*FactorialTable[i-1];
  }
  return FactorialTable;
}
const std::vector<double> FactorialFunctor::FactorialTable = FactorialTableCalculator();

vector<double> BinomialCoefficientCalculator() {
  /// We need (n+1) coefficients for each value of n from 0 (for
  /// completeness) up to 2*ellMax (hard coded in the header file).
  /// That's a total of
  ///   >>> from sympy import summation, symbols
  ///   >>> ellMax, n, k = symbols('ellMax n k', integer=True)
  ///   >>> summation(n+1, (n, 0, 2*ellMax))
  ///  2*ellMax**2 + 3*ellMax + 1
  /// With a similar calculation, we can see that the associated access
  /// operator needs element (n*(n+1)/2 + k) of the array.
  vector<double> BinomialCoefficientTable(2*ellMax*ellMax + 3*ellMax + 1);
  unsigned int i=0;
  FactorialFunctor Factorial;
  for(unsigned int n=0; n<=2*ellMax; ++n) {
    for(unsigned int k=0; k<=n; ++k) {
      BinomialCoefficientTable[i++] = std::floor(0.5+Factorial(n)/(Factorial(k)*Factorial(n-k)));
    }
  }
  return BinomialCoefficientTable;
}
const std::vector<double> BinomialCoefficientFunctor::BinomialCoefficientTable = BinomialCoefficientCalculator();

vector<double> LadderOperatorFactorCalculator() {
  /// We need (2*ell+1) coefficients for each value of ell from 0 (for
  /// completeness) up to ellMax (hard coded in the header file).
  /// That's a total of
  ///   >>> from sympy import summation, symbols
  ///   >>> ell, ellMax, m, mp = symbols('ell ellMax m mp', integer=True)
  ///   >>> summation(2*ell+1, (ell, 0, ellMax))
  ///   ellMax**2 + 2*ellMax + 1
  /// With a similar calculation, we can see that the associated access
  /// operator needs element
  ///   >>> summation(2*ell+1, (ell, 0, ell-1)) + ell + m
  ///   ell**2 + ell + m
  std::vector<double> FactorTable(ellMax*ellMax + 2*ellMax + 1);
  unsigned int i=0;
  for(int ell=0; ell<=ellMax; ++ell) {
    for(int m=-ell; m<=ell; ++m) {
      FactorTable[i++] = std::sqrt(ell*(ell+1)-m*(m+1));
    }
  }
  return FactorTable;
}
const std::vector<double> LadderOperatorFactorFunctor::FactorTable = LadderOperatorFactorCalculator();
