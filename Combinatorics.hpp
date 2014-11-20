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
    inline double operator[](const unsigned int i) const {
      #ifdef DEBUG
      if(i>170) {
        std::cerr << "\n\ni = " << i
                  << "\nIn factorials, 170 is the largest `double` that does not overflow." << std::endl;
        throw(IndexOutOfBounds);
      }
      #endif
      return FactorialTable[i];
    }
    inline double operator()(const unsigned int i) const {
      #ifdef DEBUG
      if(i>170) {
        std::cerr << "\n\ni = " << i
                  << "\nIn factorials, 170 is the largest `double` that does not overflow." << std::endl;
        throw(IndexOutOfBounds);
      }
      #endif
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
      if(n>2*ellMax || k>n) {
        std::cerr << "\n\n(n, k) = (" << n << ", " << k << ")\t2*ellMax = " << 2*ellMax
                  << "\nBinomialCoefficientFunctor is only implemented up to n=2*ellMax=" << 2*ellMax
                  << ".\nTo increase this bound, edit 'ellMax' in " << __FILE__ << " and recompile." << std::endl;
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
      if(ell>ellMax || std::abs(m)>ell) {
        std::cerr << "\n\n(ell, m) = (" << ell << ", " << m << ")\tellMax = " << ellMax
                  << "\nLadderOperatorFactorFunctor is only implemented up to ell=" << ellMax
                  << ".\nTo increase this bound, edit 'ellMax' in " << __FILE__ << " and recompile." << std::endl;
        throw(IndexOutOfBounds);
      }
      #endif
      return FactorTable[ell*ell+ell+m];
    }
  };

  // GSL is about 3 times slower than my function (which I just translated from sympy)
  // #ifdef USE_GSL
  // extern "C" {
  //   #include <gsl/gsl_sf_coupling.h>
  // }
  // inline double Wigner3j(int j_1, int j_2, int j_3, int m_1, int m_2, int m_3) {
  //   return gsl_sf_coupling_3j(2*j_1, 2*j_2, 2*j_3, 2*m_1, 2*m_2, 2*m_3);
  // }
  // #else
  double Wigner3j(int j_1, int j_2, int j_3, int m_1, int m_2, int m_3);
  // #endif

  // class Wigner3jSingleton {
  //   /// CAUTION!!!
  //   ///
  //   /// Preliminary testing shows that this is actually SLOWER than
  //   /// the direct evaluation of the Wigner3j function, at least for
  //   /// ell<=16, which is all this is currently implemented for.
  //   /// Presumably, this is because the indexing function has to
  //   /// manipulate a magic square to translate the arguments into an
  //   /// index.  In particular, it is substantially slower for small
  //   /// ell values, but seems to catch up for larger ones.
  //   ///
  //   /// Also, even if you decide to use this, this is NOT a thread-safe
  //   /// singleton.  Because of various optimizations I've added to make
  //   /// this faster, there are mutables.
  //   ///
  //   /// This class is intended to implement the storage algorithm
  //   /// introduced by Rasch and Yu [SIAM J. Sci. Comput. (2003) 25, 4,
  //   /// 1416--1428].  I may not be doing the magic-square
  //   /// manipulations in the most clever way, which might improve the
  //   /// speed of this algorithm.
  // private:
  //   static const Wigner3jSingleton* Wigner3jInstance;
  //   std::vector<double> FactorTable;
  //   mutable std::vector<std::vector<int> > MagicSquare;
  //   mutable int S,i_S,j_S,L,i_L,j_L,i_T,j_T,B,X,T,j,i;
  //   Wigner3jSingleton();
  //   Wigner3jSingleton(const Wigner3jSingleton& that) {
  //     Wigner3jInstance = that.Wigner3jInstance;
  //   }
  //   Wigner3jSingleton& operator=(const Wigner3jSingleton& that) {
  //     if(this!=&that) Wigner3jInstance = that.Wigner3jInstance;
  //     return *this;
  //   }
  //   ~Wigner3jSingleton() { }
  // public:
  //   static const Wigner3jSingleton& Instance() {
  //     static const Wigner3jSingleton Instance;
  //     Wigner3jInstance = &Instance;
  //     return *Wigner3jInstance;
  //   }
  //   double operator()(int j_1, int j_2, int j_3, int m_1, int m_2, int m_3) const;
  // };


} // namespace SphericalFunctions

#endif // COMBINATORICS_HPP
