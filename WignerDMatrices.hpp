// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#ifndef WIGNERDMATRICES_HPP
#define WIGNERDMATRICES_HPP

#include "Quaternions.hpp"
#include "Combinatorics.hpp"

namespace SphericalFunctions {

  /// Object for pre-computing and retrieving coefficients for the Wigner D matrices
  class WignerCoefficientSingleton {
  private:
    static const WignerCoefficientSingleton* WignerCoefficientInstance;
    std::vector<double> CoefficientTable;
    WignerCoefficientSingleton()
      : CoefficientTable(ellMax*(ellMax*(4*ellMax + 12) + 11)/3 + 1)
    {
      const FactorialSingleton& Factorial = FactorialSingleton::Instance();
      unsigned int i=0;
      for(int ell=0; ell<=ellMax; ++ell) {
        for(int mp=-ell; mp<=ell; ++mp) {
          for(int m=-ell; m<=ell; ++m) {
            CoefficientTable[i++] =
              std::sqrt( Factorial(ell+m)*Factorial(ell-m)
                         / double(Factorial(ell+mp)*Factorial(ell-mp)) );
          }
        }
      }
    }
    WignerCoefficientSingleton(const WignerCoefficientSingleton& that) {
      WignerCoefficientInstance = that.WignerCoefficientInstance;
    }
    WignerCoefficientSingleton& operator=(const WignerCoefficientSingleton& that) {
      if(this!=&that) WignerCoefficientInstance = that.WignerCoefficientInstance;
      return *this;
    }
    ~WignerCoefficientSingleton() { }
  public:
    static const WignerCoefficientSingleton& Instance() {
      static const WignerCoefficientSingleton Instance;
      WignerCoefficientInstance = &Instance;
      return *WignerCoefficientInstance;
    }
    inline double operator()(const int ell, const int mp, const int m) const {
      #ifdef DEBUG
      if(ell>ellMax_Utilities || std::abs(mp)>ell || std::abs(m)>ell) {
        std::cerr << "\n\n(ell, mp, m) = (" << ell << ", " << mp << ", " << m << ")\tellMax_Utilities = " << ellMax_Utilities
                  << "\nWignerCoefficientSingleton is only implemented up to ell=" << ellMax_Utilities
                  << ".\nTo increase this bound, edit 'ellMax_Utilities' in " << __FILE__ << " and recompile." << std::endl;
        throw(IndexOutOfBounds);
      }
      #endif
      return CoefficientTable[ell*(ell*(4*ell + 6) + 5)/3 + mp*(2*ell + 1) + m];
    }
  };

  /// Object for computing the Wigner D matrices as functions of quaternion rotors
  class WignerDMatrix {
  private:
    const BinomialCoefficientSingleton& BinomialCoefficient;
    const WignerCoefficientSingleton& WignerCoefficient;
    std::complex<double> Ra, Rb;
    double absRa, absRb, absRRatioSquared;
  public:
    WignerDMatrix(const Quaternions::Quaternion& iR=Quaternions::Quaternion(1,0,0,0));
    WignerDMatrix& SetRotation(const Quaternions::Quaternion& iR);
    std::complex<double> operator()(const int ell, const int mp, const int m) const;
  };

} // namespace SphericalFunctions

#endif // WIGNERDMATRICES_HPP
