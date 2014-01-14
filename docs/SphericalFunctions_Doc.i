%feature("docstring") SphericalFunctions::BinomialCoefficientFunctor::operator() """


  Parameters
  ----------
    const unsigned int n
    const unsigned int k
  
  Returns
  -------
    double
  
"""

%feature("docstring") SphericalFunctions::SWSH::SetRotation """


  Parameters
  ----------
    const Quaternions::Quaternion& iR
  
  Returns
  -------
    SWSH&
  
"""

%feature("docstring") SphericalFunctions::WignerCoefficientFunctor::operator() """


  Parameters
  ----------
    const int ell
    const int mp
    const int m
  
  Returns
  -------
    double
  
"""

%feature("docstring") SphericalFunctions::LadderOperatorFactorFunctor::LadderOperatorFactorFunctor """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    LadderOperatorFactorFunctor
  
"""

%feature("docstring") SphericalFunctions::LadderOperatorFactorFunctor """
class SphericalFunctions::LadderOperatorFactorFunctor
=====================================================
  Object for pre-computing and retrieving values of the ladder operators.
  
  Member variables
  ----------------
    const vector<double> FactorTable
  
"""

%feature("docstring") SphericalFunctions::SWSH::operator() """


  Parameters
  ----------
    const int ell
    const int m
  
  Returns
  -------
    complex<double>
  
"""

%feature("docstring") SphericalFunctions::WignerCoefficientFunctor::WignerCoefficientFunctor """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    WignerCoefficientFunctor
  
"""

%feature("docstring") FactorialTableCalculator """
Class to create an object returning the factorial of an argument.
=================================================================
  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    Note that because a double is returned, only values up to 28! will be
    exact; higher values will be accurate to machine precision. Values up to
    170! only are allowed because higher values overflow.
  
"""

%feature("docstring") SphericalFunctions::BinomialCoefficientFunctor::BinomialCoefficientFunctor """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    BinomialCoefficientFunctor
  
"""

%feature("docstring") SphericalFunctions::FactorialFunctor::operator() """


  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    double
  
"""

%feature("docstring") BinomialCoefficientCalculator """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    We need (n+1) coefficients for each value of n from 0 (for completeness) up
    to 2*ellMax (hard coded in the header file). That's a total of from sympy
    import summation, symbols ellMax, n, k = symbols('ellMax n k',
    integer=True) summation(n+1, (n, 0, 2*ellMax))
    
    2*ellMax**2 + 3*ellMax + 1 With a similar calculation, we can see that the
    associated access operator needs element (n*(n+1)/2 + k) of the array.
  
"""

%feature("docstring") SphericalFunctions::WignerCoefficientFunctor """
class SphericalFunctions::WignerCoefficientFunctor
==================================================
  Object for pre-computing and retrieving coefficients for the Wigner D
  matrices.
  
  Member variables
  ----------------
    const vector<double> CoefficientTable
  
"""

%feature("docstring") SphericalFunctions::BinomialCoefficientFunctor """
class SphericalFunctions::BinomialCoefficientFunctor
====================================================
  Object for pre-computing and retrieving binomials.
  
  Member variables
  ----------------
    const vector<double> BinomialCoefficientTable
  
"""

%feature("docstring") WignerDMatrix::SetRotation """
Reset the rotor for this object to the given value.
===================================================
  Parameters
  ----------
    const Quaternions::Quaternion& iR
  
  Returns
  -------
    WignerDMatrix&
  
"""

%feature("docstring") WignerDMatrix::operator() """
Evaluate the D matrix element for the given (ell, mp, m) indices.
=================================================================
  Parameters
  ----------
    const int ell
    const int mp
    const int m
  
  Returns
  -------
    complex<double>
  
"""

%feature("docstring") WignerCoefficientCalculator """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    We need (2*ell+1)*(2*ell+1) coefficients for each value of ell from 0 (for
    completenes) up to ellMax (hard coded in the header file). That's a total
    of from sympy import summation, symbols, simplify from
    sympy.polys.polyfuncs import horner ell, ellMax, m, mp = symbols('ell
    ellMax m mp', integer=True) horner(simplify(summation((2*ell+1)**2, (ell,
    0, ellMax))))
    
    ellMax*(ellMax*(4*ellMax/3 + 4) + 11/3) + 1 With a similar calculation, we
    can see that the associated access operator needs element
    horner(summation((2*ell+1)**2, (ell, 0, ell-1)) + (2*ell+1)*(ell+mp) + ell
    + m)
    
    ell*(ell*(4*ell/3 + 2) + 5/3) + mp*(2*ell + 1) + m of the array.
  
"""

%feature("docstring") SphericalFunctions::FactorialFunctor """
class SphericalFunctions::FactorialFunctor
==========================================
  Object for pre-computing and retrieving factorials.
  
  Member variables
  ----------------
    const vector<double> FactorialTable
  
"""

%feature("docstring") SphericalFunctions::SWSH """
class SphericalFunctions::SWSH
==============================
  Object for computing values of the spin-weighted spherical harmonics (SWSHs)
  
  Member variables
  ----------------
    WignerDMatrix D

            Note that this object is a functor taking a quaternion argument.
      Ordinarily, we think of this quaternion as being the rotor taking the
      unit $z$ vector into a point $(\\vartheta, \\varphi)$, which gives us the
      usual form of SWSHs. However, more general arguments are possible; no
      checking is done to ensure that the argument has the simple form of a
      minimal rotation to the spherical coordinate.    int spin
    double sign
  
"""

%feature("docstring") std """
namespace std
=============
  STL namespace.
  
"""

%feature("docstring") LadderOperatorFactorCalculator """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    vector<double>
  
  Description
  -----------
    We need (2*ell+1) coefficients for each value of ell from 0 (for
    completeness) up to ellMax (hard coded in the header file). That's a total
    of from sympy import summation, symbols ell, ellMax, m, mp = symbols('ell
    ellMax m mp', integer=True) summation(2*ell+1, (ell, 0, ellMax))
    
    ellMax**2 + 2*ellMax + 1 With a similar calculation, we can see that the
    associated access operator needs element summation(2*ell+1, (ell, 0,
    ell-1)) + ell + m
    
    ell**2 + ell + m
  
"""

%feature("docstring") SphericalFunctions::LadderOperatorFactorFunctor::operator() """


  Parameters
  ----------
    const int ell
    const int m
  
  Returns
  -------
    double
  
"""

%feature("docstring") SphericalFunctions """
namespace SphericalFunctions
============================
"""

%feature("docstring") SphericalFunctions::WignerDMatrix """
class SphericalFunctions::WignerDMatrix
=======================================
  Object for computing the Wigner D matrices as functions of quaternion rotors.
  
  Member variables
  ----------------
    BinomialCoefficientFunctor BinomialCoefficient
    WignerCoefficientFunctor WignerCoefficient
    complex<double> Ra
    complex<double> Rb
    double absRa
    double absRb
    double absRRatioSquared
  
"""

%feature("docstring") SphericalFunctions::SWSH::SetAngles """


  Parameters
  ----------
    const double vartheta
    const double varphi
  
  Returns
  -------
    SWSH&
  
"""

%feature("docstring") WignerDMatrix::WignerDMatrix """
Construct the D matrix object given the (optional) rotor.
=========================================================
  Parameters
  ----------
    const Quaternions::Quaternion& iR = Quaternions::Quaternion(1, 0, 0, 0)
  
  Returns
  -------
    WignerDMatrix
  
"""

%feature("docstring") SphericalFunctions::FactorialFunctor::FactorialFunctor """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    FactorialFunctor
  
"""

