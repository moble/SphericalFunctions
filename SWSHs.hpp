// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#ifndef SWSHS_HPP
#define SWSHS_HPP

#include "WignerDMatrices.hpp"

namespace SphericalFunctions {

  /// Object for computing values of the spin-weighted spherical harmonics (SWSHs)
  class SWSH {
    /// Note that this object is a functor taking a quaternion
    /// argument.  Ordinarily, we think of this quaternion as being
    /// the rotor taking the unit \f$z\f$ vector into a point
    /// \f$(\vartheta, \varphi)\f$, which gives us the usual form of
    /// SWSHs.  However, more general arguments are possible; no
    /// checking is done to ensure that the argument has the simple
    /// form of a minimal rotation to the spherical coordinate.
  private:
    WignerDMatrix D;
    int spin;
    double sign;
  public:
    // / \@cond
    SWSH(const int s, const Quaternions::Quaternion& iR=Quaternions::Quaternion(1,0,0,0))
      : D(iR), spin(s), sign(s%2==0 ? 1.0 : -1.0)
    { }
    // / \endcond
    void RaiseErrorOnBadIndices(const bool ErrorOnBadIndices=true) { D.ErrorOnBadIndices = ErrorOnBadIndices; }
    inline SWSH& SetRotation(const Quaternions::Quaternion& iR) { D.SetRotation(iR); return *this; }
    inline SWSH& SetAngles(const double vartheta, const double varphi) { D.SetRotation(Quaternions::Quaternion(vartheta, varphi)); return *this; }
    inline std::complex<double> operator()(const int ell, const int m) const {
      return sign * std::sqrt((2*ell+1)/(4*M_PI)) * D(ell, m, -spin);
    }
    std::complex<double> Evaluate(const std::vector<std::complex<double> >& Modes) const;
  };

} // namespace SphericalFunctions

#endif // SWSHS_HPP
