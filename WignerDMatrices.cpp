// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include "WignerDMatrices.hpp"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "Errors.hpp"

using namespace SphericalFunctions;
using Quaternions::Quaternion;
using std::vector;
using std::complex;
using std::cerr;
using std::endl;

#define INFOTOCERR std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": "


const WignerCoefficientSingleton* WignerCoefficientSingleton::WignerCoefficientInstance = NULL;


/// Construct the D matrix object given the (optional) rotor.
WignerDMatrix::WignerDMatrix(const Quaternion& R)
  : ErrorOnBadIndices(true),
    BinomialCoefficient(BinomialCoefficientSingleton::Instance()),
    WignerCoefficient(WignerCoefficientSingleton::Instance()),
    Ra(R[0], R[3]), Rb(R[2], R[1]),
    absRa(abs(Ra)), absRb(abs(Rb)),
    absRRatioSquared(absRa>epsilon ? absRb*absRb/(absRa*absRa) : 0.0),
    intlog10absRa(absRa>epsilon ? std::log10(absRa) : DBL_MIN_10_EXP),
    intlog10absRb(absRb>epsilon ? std::log10(absRb) : DBL_MIN_10_EXP)
{ }

/// Reset the rotor for this object to the given value.
WignerDMatrix& WignerDMatrix::SetRotation(const Quaternion& R) {
  Ra = std::complex<double>(R[0], R[3]);
  Rb = std::complex<double>(R[2], R[1]);
  absRa = abs(Ra);
  absRb = abs(Rb);
  absRRatioSquared = absRa>epsilon ? absRb*absRb/(absRa*absRa) : 0.0;
  intlog10absRa = absRa>epsilon ? std::log10(absRa) : DBL_MIN_10_EXP;
  intlog10absRb = absRb>epsilon ? std::log10(absRb) : DBL_MIN_10_EXP;
  return *this;
}

/// Evaluate the D matrix element for the given (ell, mp, m) indices.
std::complex<double> WignerDMatrix::operator()(const int ell, const int mp, const int m) const {
  // If either sub-rotor, when raised to the exponent (mp-m), and
  // multiplied by anything within machine precision of 1, is not
  // representable as a floating-point number, just treat it as zero.
  if(std::abs(mp)>ell || std::abs(m)>ell) {
    if(ErrorOnBadIndices) {
      INFOTOCERR << "(" << ell << ", " << mp << ", " << m << ") is not a valid set of indices.\n"
                 << "If you want this object (let's call it `D`) to return 0.0 when invalid\n"
                 << " indices are requested, you can set `D.ErrorOnBadIndices = false`.\n" << std::endl;
      throw(ValueError);
    } else {
      return std::complex<double>(0.0, 0.0);
    }
  }
  if(absRa < epsilon || 2*intlog10absRa*(mp-m)<DBL_MIN_10_EXP+17) {
    return (mp!=-m ? 0.0 : ((ell+mp)%2==0 ? 1.0 : -1.0) * std::pow(Rb, 2*m) );
  }
  if(absRb < epsilon || 2*intlog10absRb*(mp-m)<DBL_MIN_10_EXP+17) {
    return (mp!=m ? 0.0 : std::pow(Ra, 2*m) );
  }
  const int rhoMin = std::max(0,mp-m);
  const int rhoMax = std::min(ell+mp,ell-m);
  if(absRa < 1.e-3) { // Deal with NANs in certain cases
    const std::complex<double> Prefactor =
      WignerCoefficient(ell, mp, m) * std::pow(Ra, m+mp) * std::pow(Rb, m-mp);
    const double absRaSquared = absRa*absRa;
    const double absRbSquared = absRb*absRb;
    double Sum = 0.0;
    for(int rho=rhoMax; rho>=rhoMin; --rho) {
      const double aTerm = std::pow(absRaSquared, ell-m-rho);
      if(aTerm != aTerm || aTerm<1.e-100) { // This assumes --fast-math is off
        Sum *= absRbSquared;
        continue;
      }
      Sum = ( (rho%2==0 ? 1 : -1) * BinomialCoefficient(ell+mp,rho) * BinomialCoefficient(ell-mp, ell-rho-m) * aTerm )
        + ( Sum * absRbSquared );
    }
    return Prefactor * Sum * std::pow(absRbSquared, rhoMin);
  } else {
    const std::complex<double> Prefactor =
      (WignerCoefficient(ell, mp, m) * std::pow(absRa, 2*ell-2*m))
      * std::pow(Ra, m+mp) * std::pow(Rb, m-mp);
    double Sum = 0.0;
    for(int rho=rhoMax; rho>=rhoMin; --rho) {
      Sum = ( (rho%2==0 ? 1 : -1) * BinomialCoefficient(ell+mp,rho) * BinomialCoefficient(ell-mp, ell-rho-m) )
        + ( Sum * absRRatioSquared );
    }
    return Prefactor * Sum * std::pow(absRRatioSquared, rhoMin);
  }
}
