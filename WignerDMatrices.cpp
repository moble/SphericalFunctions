// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include "WignerDMatrices.hpp"
#include "Errors.hpp"

using namespace SphericalFunctions;
using Quaternions::Quaternion;
using std::vector;
using std::complex;
using std::cerr;
using std::endl;

const WignerCoefficientSingleton* WignerCoefficientSingleton::WignerCoefficientInstance = NULL;


/// Construct the D matrix object given the (optional) rotor.
WignerDMatrix::WignerDMatrix(const Quaternion& R)
  : BinomialCoefficient(BinomialCoefficientSingleton::Instance()),
    WignerCoefficient(WignerCoefficientSingleton::Instance()),
    Ra(R[0], R[3]), Rb(R[2], R[1]),
    absRa(abs(Ra)), absRb(abs(Rb)), absRRatioSquared(absRb*absRb/(absRa*absRa))
{ }

/// Reset the rotor for this object to the given value.
WignerDMatrix& WignerDMatrix::SetRotation(const Quaternion& R) {
  Ra = std::complex<double>(R[0], R[3]);
  Rb = std::complex<double>(R[2], R[1]);
  absRa = abs(Ra);
  absRb = abs(Rb);
  absRRatioSquared = absRb*absRb/(absRa*absRa);
  return *this;
}

/// Evaluate the D matrix element for the given (ell, mp, m) indices.
std::complex<double> WignerDMatrix::operator()(const int ell, const int mp, const int m) const {
  if(absRa < epsilon) {
    return (mp!=-m ? 0.0 : ((ell+mp)%2==0 ? 1.0 : -1.0) * std::pow(Rb, 2*m) );
  }
  if(absRb < epsilon) {
    return (mp!=m ? 0.0 : std::pow(Ra, 2*m) );
  }
  if(absRa < 1.e-3) { // Deal with NANs in certain cases
    const std::complex<double> Prefactor =
      WignerCoefficient(ell, mp, m) * std::pow(Ra, m+mp) * std::pow(Rb, m-mp);
    const int rhoMin = std::max(0,mp-m);
    const int rhoMax = std::min(ell+mp,ell-m);
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
  }
  const std::complex<double> Prefactor =
    (WignerCoefficient(ell, mp, m) * std::pow(absRa, 2*ell-2*m))
    * std::pow(Ra, m+mp) * std::pow(Rb, m-mp);
  const int rhoMin = std::max(0,mp-m);
  const int rhoMax = std::min(ell+mp,ell-m);
  double Sum = 0.0;
  for(int rho=rhoMax; rho>=rhoMin; --rho) {
    Sum = ( (rho%2==0 ? 1 : -1) * BinomialCoefficient(ell+mp,rho) * BinomialCoefficient(ell-mp, ell-rho-m) )
      + ( Sum * absRRatioSquared );
  }
  return Prefactor * Sum * std::pow(absRRatioSquared, rhoMin);
}
