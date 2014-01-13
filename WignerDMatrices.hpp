// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#ifndef WIGNERDMATRICES_HPP
#define WIGNERDMATRICES_HPP

#include "Quaternions.hpp"
#include "Combinatorics.hpp"

namespace SphericalFunctions {

  /// Object for pre-computing and retrieving coefficients for the Wigner D matrices
  class WignerCoefficientFunctor {
  private:
    static const std::vector<double> CoefficientTable;
  public:
    WignerCoefficientFunctor() { };
    inline double operator()(const int ell, const int mp, const int m) const {
      #ifdef DEBUG
      if(ell>ellMax_Utilities || std::abs(mp)>ell || std::abs(m)>ell) {
	std::cerr << "\n\n(ell, mp, m) = (" << ell << ", " << mp << ", " << m << ")\tellMax_Utilities = " << ellMax_Utilities
		  << "\nWignerCoefficientFunctor is only implemented up to ell=" << ellMax_Utilities
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
    BinomialCoefficientFunctor BinomialCoefficient;
    WignerCoefficientFunctor WignerCoefficient;
    std::complex<double> Ra, Rb;
    double absRa, absRb, absRRatioSquared;
  public:
    WignerDMatrix(const Quaternions::Quaternion& iR=Quaternions::Quaternion(1,0,0,0));
    WignerDMatrix& SetRotation(const Quaternions::Quaternion& iR);
    std::complex<double> operator()(const int ell, const int mp, const int m) const;
  };

} // namespace SphericalFunctions

#endif // WIGNERDMATRICES_HPP
