// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#ifndef COMBINATORICS_HPP
#define COMBINATORICS_HPP

#include <vector>
#include <cmath>

#ifdef DEBUG
#include <iostream>
#include "Errors.hpp"
#endif

namespace SphericalFunctions {

  const int ellMax = 32;
  const double epsilon = 1.0e-14;

  /// Object for pre-computing and retrieving factorials
  class FactorialFunctor {
  private:
    std::vector<double> FactorialTable;
  public:
    FactorialFunctor()
      : FactorialTable(171)
    {
      FactorialTable[0] = 1.0;
      for (int i=1;i<171;i++) {
        FactorialTable[i] = i*FactorialTable[i-1];
      }
    };
    inline double operator()(const unsigned int i) const {
      #ifdef DEBUG
      if(i>171) {
        std::cerr << "\n\ni = " << i << "\tiMax = 171"
                  << "\nFactorialFunctor is only implemented up to 171!; larger values overflow."
                  << std::endl;
        throw(IndexOutOfBounds);
      }
      #endif
      return FactorialTable[i];
    }
  };

  /// Object for pre-computing and retrieving binomials
  class BinomialCoefficientFunctor {
  private:
    std::vector<double> BinomialCoefficientTable;
  public:
    BinomialCoefficientFunctor()
      : BinomialCoefficientTable(2*ellMax*ellMax + 3*ellMax + 1)
    {
      unsigned int i=0;
      FactorialFunctor Factorial;
      for(unsigned int n=0; n<=2*ellMax; ++n) {
        for(unsigned int k=0; k<=n; ++k) {
          BinomialCoefficientTable[i++] = std::floor(0.5+Factorial(n)/(Factorial(k)*Factorial(n-k)));
        }
      }
    };
    inline double operator()(const unsigned int n, const unsigned int k) const {
      #ifdef DEBUG
      if(n>2*ellMax_Utilities || k>n) {
        std::cerr << "\n\n(n, k) = (" << n << ", " << k << ")\t2*ellMax_Utilities = " << 2*ellMax_Utilities
                  << "\nBinomialCoefficientFunctor is only implemented up to n=2*ellMax_Utilities=" << 2*ellMax_Utilities
                  << ".\nTo increase this bound, edit 'ellMax_Utilities' in " << __FILE__ << " and recompile." << std::endl;
        throw(IndexOutOfBounds);
      }
      #endif
      return BinomialCoefficientTable[(n*(n+1))/2+k];
    }
  };

  /// Object for pre-computing and retrieving values of the ladder operators
  class LadderOperatorFactorFunctor {
  private:
    std::vector<double> FactorTable;
  public:
    LadderOperatorFactorFunctor()
      : FactorTable(ellMax*ellMax + 2*ellMax + 1)
    {
      unsigned int i=0;
      for(int ell=0; ell<=ellMax; ++ell) {
        for(int m=-ell; m<=ell; ++m) {
          FactorTable[i++] = std::sqrt(ell*(ell+1)-m*(m+1));
        }
      }
    };
    inline double operator()(const int ell, const int m) const {
      #ifdef DEBUG
      if(ell>ellMax_Utilities || std::abs(m)>ell) {
        std::cerr << "\n\n(ell, m) = (" << ell << ", " << m << ")\tellMax_Utilities = " << ellMax_Utilities
                  << "\nLadderOperatorFactorFunctor is only implemented up to ell=" << ellMax_Utilities
                  << ".\nTo increase this bound, edit 'ellMax_Utilities' in " << __FILE__ << " and recompile." << std::endl;
        throw(IndexOutOfBounds);
      }
      #endif
      return FactorTable[ell*ell+ell+m];
    }
  };

} // namespace SphericalFunctions

#endif // COMBINATORICS_HPP
