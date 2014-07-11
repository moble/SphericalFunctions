%feature("docstring") abs """


  Parameters
  ----------
    const int a
  
  Returns
  -------
    int
  
"""

%feature("docstring") SphericalFunctions::LadderOperatorFactorSingleton::operator() """


  Parameters
  ----------
    const int ell
    const int m
  
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

%feature("docstring") SphericalFunctions::BinomialCoefficientSingleton::Instance """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    const BinomialCoefficientSingleton&
  
"""

%feature("docstring") SphericalFunctions::SWSH::Evaluate """
Evaluate Modes at the point given by the rotor.
===============================================
  Parameters
  ----------
    const vector<complex<double>>& Modes
      vector<complex<double>> in spinsfast order (which include ell=0, etc.)
  
  Returns
  -------
    complex<double>
  
  Description
  -----------
    Assuming the current rotor (or angles), evaluate the given data at that
    point. (See the SWSHs class documentation for more detail regarding the
    rotor.)
    
    Note that the input data are assumed to start with ell=0, incrementing m
    from -ell to ell, then incrementing ell. This is assumed even if the spin
    weight is nonzero, so ell=0 modes must be nonzero, to agree with the order
    of data required for the spinsfast module. Any data in these 'impossible'
    slots will just be ignored.
    
    In many cases, evaluating with spinsfast itself should be preferable.
    However, the major advantage of this code is that it allows for evaluation
    with respect to a rotor. For functions of nonzero spin weight, this
    automatically incorporates the spin orientation of the relevant tetrad,
    which is necessary, for example when imposing a boost.
  
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

%feature("docstring") N_ellm """
Return the size of the array needed to express this ell.
========================================================
  Parameters
  ----------
    const int ell
  
  Returns
  -------
    int
  
"""

%feature("docstring") SphericalFunctions::WignerCoefficientSingleton """
class SphericalFunctions::WignerCoefficientSingleton
====================================================
  Object for pre-computing and retrieving coefficients for the Wigner D
  matrices.
  
  Member variables
  ----------------
    const WignerCoefficientSingleton * WignerCoefficientInstance
    vector<double> CoefficientTable
  
  Non-public member functions
  ---------------------------
     WignerCoefficientSingleton
     WignerCoefficientSingleton
    WignerCoefficientSingleton& operator=
     ~WignerCoefficientSingleton
  
"""

%feature("docstring") SphericalFunctions::BinomialCoefficientSingleton """
class SphericalFunctions::BinomialCoefficientSingleton
======================================================
  Object for pre-computing and retrieving binomials.
  
  Member variables
  ----------------
    const BinomialCoefficientSingleton * BinomialCoefficientInstance
    vector<double> BinomialCoefficientTable
  
  Non-public member functions
  ---------------------------
     BinomialCoefficientSingleton
     BinomialCoefficientSingleton
    BinomialCoefficientSingleton& operator=
     ~BinomialCoefficientSingleton
  
"""

%feature("docstring") SphericalFunctions::FactorialSingleton::operator() """


  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    double
  
"""

%feature("docstring") SphericalFunctions::Wigner3j """
Evaluate Wigner's 3-j symbol.
=============================
  Parameters
  ----------
    int j_1
    int j_2
    int j_3
    int m_1
    int m_2
    int m_3
  
  Returns
  -------
    double
  
  Description
  -----------
    This function is just a translation of the function in
    sympy.physics.wigner, written by Jens Rasch.
  
"""

%feature("docstring") SphericalFunctions::FactorialSingleton::Instance """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    const FactorialSingleton&
  
"""

%feature("docstring") SphericalFunctions::WignerCoefficientSingleton::operator() """


  Parameters
  ----------
    const int ell
    const int mp
    const int m
  
  Returns
  -------
    double
  
"""

%feature("docstring") SphericalFunctions::SWSH::RaiseErrorOnBadIndices """


  Parameters
  ----------
    const bool ErrorOnBadIndices = true
  
  Returns
  -------
    void
  
"""

%feature("docstring") SphericalFunctions::WignerDMatrix::SetRotation """


  Parameters
  ----------
    const Quaternions::Quaternion& iR
  
  Returns
  -------
    WignerDMatrix&
  



  Parameters
  ----------
    const double alpha
    const double beta
    const double gamma
  
  Returns
  -------
    WignerDMatrix&
  
"""

%feature("docstring") SphericalFunctions::WignerCoefficientSingleton::Instance """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    const WignerCoefficientSingleton&
  
"""

%feature("docstring") SphericalFunctions::LadderOperatorFactorSingleton::Instance """


  Parameters
  ----------
    (none)
  
  Returns
  -------
    const LadderOperatorFactorSingleton&
  
"""

%feature("docstring") SphericalFunctions::BinomialCoefficientSingleton::operator() """


  Parameters
  ----------
    const unsigned int n
    const unsigned int k
  
  Returns
  -------
    double
  
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

%feature("docstring") SphericalFunctions::LadderOperatorFactorSingleton """
class SphericalFunctions::LadderOperatorFactorSingleton
=======================================================
  Member variables
  ----------------
    const LadderOperatorFactorSingleton * LadderOperatorFactorInstance
    vector<double> FactorTable
  
  Non-public member functions
  ---------------------------
     LadderOperatorFactorSingleton
     LadderOperatorFactorSingleton
    LadderOperatorFactorSingleton& operator=
     ~LadderOperatorFactorSingleton
  
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
    bool ErrorOnBadIndices

            Note that this object is a functor. The rotation should be set, and then
      components can be taken by calling the object with arguments (ell,mp,m).
      The rotation can then be set to another value, and the process repeated.
      Evaluation in this order is more efficient than the other way around.    const BinomialCoefficientSingleton& BinomialCoefficient
    const WignerCoefficientSingleton& WignerCoefficient
    complex<double> Ra
    complex<double> Rb
    double absRa
    double absRb
    double absRRatioSquared
    int intlog10absRa
    int intlog10absRb
  
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

%feature("docstring") SphericalFunctions::FactorialSingleton """
class SphericalFunctions::FactorialSingleton
============================================
  Member variables
  ----------------
    const FactorialSingleton * FactorialInstance
    vector<double> FactorialTable
  
  Non-public member functions
  ---------------------------
     FactorialSingleton
     FactorialSingleton
    FactorialSingleton& operator=
     ~FactorialSingleton
  
"""

%feature("docstring") SphericalFunctions::FactorialSingleton::operator[] """


  Parameters
  ----------
    const unsigned int i
  
  Returns
  -------
    double
  
"""

