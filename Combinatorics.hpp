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

  class FactorialSingleton {
  private:
    static const FactorialSingleton* FactorialInstance;
    std::vector<double> FactorialTable;
    FactorialSingleton() : FactorialTable(171) {
      FactorialTable[0] = 1.0;
      for (int i=1;i<171;i++) {
        FactorialTable[i] = i*FactorialTable[i-1];
      }
    }
    FactorialSingleton(const FactorialSingleton& that) {
      FactorialInstance = that.FactorialInstance;
    }
    FactorialSingleton& operator=(const FactorialSingleton& that) {
      if(this!=&that) FactorialInstance = that.FactorialInstance;
      return *this;
    }
    ~FactorialSingleton() { }
  public:
    static const FactorialSingleton& Instance() {
      static const FactorialSingleton Instance;
      FactorialInstance = &Instance;
      return *FactorialInstance;
    }
    inline double operator()(const unsigned int i) const {
      return FactorialTable[i];
    }
  }; // class FactorialSingleton

  /// Object for pre-computing and retrieving binomials
  class BinomialCoefficientSingleton {
  private:
    static const BinomialCoefficientSingleton* BinomialCoefficientInstance;
    std::vector<double> BinomialCoefficientTable;
    BinomialCoefficientSingleton()
      : BinomialCoefficientTable(2*ellMax*ellMax + 3*ellMax + 1)
    {
      unsigned int i=0;
      const FactorialSingleton& Factorial = FactorialSingleton::Instance();
      for(unsigned int n=0; n<=2*ellMax; ++n) {
        for(unsigned int k=0; k<=n; ++k) {
          BinomialCoefficientTable[i++] = std::floor(0.5+Factorial(n)/(Factorial(k)*Factorial(n-k)));
        }
      }
    }
    BinomialCoefficientSingleton(const BinomialCoefficientSingleton& that) {
      BinomialCoefficientInstance = that.BinomialCoefficientInstance;
    }
    BinomialCoefficientSingleton& operator=(const BinomialCoefficientSingleton& that) {
      if(this!=&that) BinomialCoefficientInstance = that.BinomialCoefficientInstance;
      return *this;
    }
    ~BinomialCoefficientSingleton() { }
  public:
    static const BinomialCoefficientSingleton& Instance() {
      static const BinomialCoefficientSingleton Instance;
      BinomialCoefficientInstance = &Instance;
      return *BinomialCoefficientInstance;
    }
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

  class LadderOperatorFactorSingleton {
  private:
    static const LadderOperatorFactorSingleton* LadderOperatorFactorInstance;
    std::vector<double> FactorTable;
    LadderOperatorFactorSingleton()
      : FactorTable(ellMax*ellMax + 2*ellMax + 1)
    {
      unsigned int i=0;
      for(int ell=0; ell<=ellMax; ++ell) {
        for(int m=-ell; m<=ell; ++m) {
          FactorTable[i++] = std::sqrt(ell*(ell+1)-m*(m+1));
        }
      }
    }
    LadderOperatorFactorSingleton(const LadderOperatorFactorSingleton& that) {
      LadderOperatorFactorInstance = that.LadderOperatorFactorInstance;
    }
    LadderOperatorFactorSingleton& operator=(const LadderOperatorFactorSingleton& that) {
      if(this!=&that) LadderOperatorFactorInstance = that.LadderOperatorFactorInstance;
      return *this;
    }
    ~LadderOperatorFactorSingleton() { }
  public:
    static const LadderOperatorFactorSingleton& Instance() {
      static const LadderOperatorFactorSingleton Instance;
      LadderOperatorFactorInstance = &Instance;
      return *LadderOperatorFactorInstance;
    }
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
