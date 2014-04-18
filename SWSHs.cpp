// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include "SWSHs.hpp"

int abs(const int a) {
  return (a>0 ? -a : a);
}

/// Return the size of the array needed to express this ell
inline int N_ellm(const int ell) {
  return 1 + 2*ell + ell*ell;
}

/// Evaluate Modes at the point given by the rotor.
std::complex<double> SphericalFunctions::SWSH::Evaluate(const std::vector<std::complex<double> >& Modes) const {
  ///
  /// \param Modes vector<complex<double> > in spinsfast order (which include ell=0, etc.)
  ///
  /// Assuming the current rotor (or angles), evaluate the given data
  /// at that point.  (See the `SWSHs` class documentation for more
  /// detail regarding the rotor.)
  ///
  /// Note that the input data are assumed to start with ell=0,
  /// incrementing m from -ell to ell, then incrementing ell.  This is
  /// assumed even if the spin weight is nonzero, so ell=0 modes must
  /// be nonzero, to agree with the order of data required for the
  /// `spinsfast` module.  Any data in these "impossible" slots will
  /// just be ignored.
  ///
  /// In many cases, evaluating with spinsfast itself should be
  /// preferable.  However, the major advantage of this code is that
  /// it allows for evaluation with respect to a rotor.  For functions
  /// of nonzero spin weight, this automatically incorporates the spin
  /// orientation of the relevant tetrad, which is necessary, for
  /// example when imposing a boost.

  std::complex<double> r(0.0);

  int i=N_ellm(abs(spin)-1);
  for(int ell=abs(spin); i<Modes.size(); ++ell) {
    for(int m=-ell; (m<=ell && i<Modes.size()); ++m, ++i) {
      r += Modes[i]*(*this)(ell,m);
    }
  }

  return r;
}

