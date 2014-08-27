// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include "Combinatorics.hpp"
#include <iostream>
#include <complex>
#include <cstdlib>
#include <cmath>
#include "Errors.hpp"

// This macro is useful for debugging
#define INFOTOCERR std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": "

using namespace SphericalFunctions;
using std::vector;

const FactorialSingleton* FactorialSingleton::FactorialInstance = 0;
const BinomialCoefficientSingleton* BinomialCoefficientSingleton::BinomialCoefficientInstance = 0;
const LadderOperatorFactorSingleton* LadderOperatorFactorSingleton::LadderOperatorFactorInstance = 0;

// #ifndef USE_GSL
/// Evaluate Wigner's 3-j symbol
double SphericalFunctions::Wigner3j(int j_1, int j_2, int j_3, int m_1, int m_2, int m_3) {
  /// This function is just a translation of the function in
  /// `sympy.physics.wigner`, written by Jens Rasch.
  #ifdef DEBUG
  if(int(j_1 * 2) != j_1 * 2 || int(j_2 * 2) != j_2 * 2 || int(j_3 * 2) != j_3 * 2) {
    INFOTOCERR << "\n\n(j_1,j_2,j_3,m_1,m_2,m_3) = (" << j_1 << ","  << j_2 << ","  << j_3 << "," << m_1 << "," << m_2 << "," << m_3
               << "); j values must be integer or half integer" << std::endl;
    throw(ValueError);
  }
  if(int(m_1 * 2) != m_1 * 2 || int(m_2 * 2) != m_2 * 2 || int(m_3 * 2) != m_3 * 2) {
    INFOTOCERR << "\n\n(j_1,j_2,j_3,m_1,m_2,m_3) = (" << j_1 << ","  << j_2 << ","  << j_3 << "," << m_1 << "," << m_2 << "," << m_3
               << "); m values must be integer or half integer" << std::endl;
    throw(ValueError);
  }
  #endif

  if(m_1 + m_2 + m_3 != 0) {
    return 0;
  }
  // const int prefid = std::pow(-1,int(j_1 - j_2 - m_3));
  const int prefid = ((j_1 - j_2 - m_3)%2==0 ? 1 : -1);
  m_3 = -m_3;
  const int a1 = j_1 + j_2 - j_3;
  if(a1 < 0) {
    return 0;
  }
  const int a2 = j_1 - j_2 + j_3;
  if(a2 < 0) {
    return 0;
  }
  const int a3 = -j_1 + j_2 + j_3;
  if(a3 < 0) {
    return 0;
  }
  if( (std::abs(m_1) > j_1) || (std::abs(m_2) > j_2) || (std::abs(m_3) > j_3) ) {
    return 0;
  }

  const FactorialSingleton& _Factlist = FactorialSingleton::Instance();

  const double argsqrt = double(_Factlist[int(j_1 + j_2 - j_3)] *
                                _Factlist[int(j_1 - j_2 + j_3)] *
                                _Factlist[int(-j_1 + j_2 + j_3)] *
                                _Factlist[int(j_1 - m_1)] *
                                _Factlist[int(j_1 + m_1)] *
                                _Factlist[int(j_2 - m_2)] *
                                _Factlist[int(j_2 + m_2)] *
                                _Factlist[int(j_3 - m_3)] *
                                _Factlist[int(j_3 + m_3)])
    / double(_Factlist[int(j_1 + j_2 + j_3 + 1)]);

  const double ressqrt = std::sqrt(argsqrt); // std::sqrt(std::complex<double>(argsqrt,0.0)).real();

  int imin = std::max(-j_3 + j_1 + m_2, std::max(-j_3 + j_2 - m_1, 0));
  int imax = std::min(j_2 + m_2, std::min(j_1 - m_1, j_1 + j_2 - j_3));
  double sumres = 0;
  for(int ii=imin; ii<imax+1; ++ii) {
    const double den = _Factlist[ii] *
      _Factlist[int(ii + j_3 - j_1 - m_2)] *
      _Factlist[int(j_2 + m_2 - ii)] *
      _Factlist[int(j_1 - ii - m_1)] *
      _Factlist[int(ii + j_3 - j_2 + m_1)] *
      _Factlist[int(j_1 + j_2 - j_3 - ii)];
    // sumres = sumres + pow(-1, ii) / den;
    if(ii%2==0) {
      sumres = sumres + 1.0 / den;
    } else {
      sumres = sumres - 1.0 / den;
    }
  }

  return ressqrt * sumres * prefid;
}
// #endif

// const int W3jEllMax = ellMax/2;

// const Wigner3jSingleton* Wigner3jSingleton::Wigner3jInstance = 0;

// SphericalFunctions::Wigner3jSingleton::Wigner3jSingleton()
//   : FactorTable(ellMax*(274+ellMax*(225+ellMax*(85+ellMax*(15+ellMax))))/120+1), // Eq. (2.14)
//     MagicSquare(3, std::vector<int>(3))
// {
//   // std::cerr << "FactorTable.size()=" << FactorTable.size() << std::endl;
//   // The following logic (and equation references) come from Rasch and Yu (2004)
//   int i=0;
//   for(int L=0; L<=ellMax; ++L) {
//     for(int X=0; X<=L; ++X) {
//       for(int T=0; T<=X; ++T) {
//         for(int B=0; B<=T; ++B) {
//           for(int S=0; S<=B; ++S) {
//             const int j_3 = (S+L)/2;
//             const int j_2 = (S+(X+B-T))/2;
//             const int j_1 = (L+(X+B-T))/2;
//             const int m_1 = j_1-X;
//             const int m_2 = j_2-B;
//             const int m_3 = T-j_3;
//             // std::cerr << "(" << L << " " << X << " " << T << " " << B << " " << S << "); "
//             //           << "(" << j_1 << " " << j_2 << " " << j_3 << " " << m_1 << " " << m_2 << " " << m_3 << ") " << std::endl;
//             FactorTable[i++] = Wigner3j(j_1, j_2, j_3, m_1, m_2, m_3);
//           }
//         }
//       }
//     }
//   }
// }

// double SphericalFunctions::Wigner3jSingleton::operator()(int j_1, int j_2, int j_3, int m_1, int m_2, int m_3) const {
//   // INFOTOCERR << "(" << j_1 << " " << j_2 << " " << j_3 << " " << m_1 << " " << m_2 << " " << m_3 << ") " << std::endl;
//   #ifdef DEBUG
//   if(j_1>W3jEllMax || j_2>W3jEllMax || j_3>W3jEllMax) {
//     std::cerr << "\n\n(j_1,j_2,j_3) = (" << j_1 << ","  << j_2 << ","  << j_3
//               << "); j values must be less than or equal to W3jEllMax=" << W3jEllMax << std::endl;
//     throw(ValueError);
//   }
//   if(std::abs(m_1)>j_1 || std::abs(m_2)>j_2 || std::abs(m_3)>j_3) {
//     std::cerr << "\n\n(j_1,j_2,j_3) = (" << j_1 << ","  << j_2 << ","  << j_3 << "; (m_1,m_2,m_3) = (" << m_1 << ","  << m_2 << ","  << m_3
//               << "); j and m values do not make sense" << std::endl;
//     throw(ValueError);
//   }
//   #endif

//   if(m_1+m_2+m_3!=0) {
//     return 0.;
//   }

//   // The following logic (and equation references) come from Rasch and Yu (2004)
//   MagicSquare[0][0] = -j_1+j_2+j_3;
//   MagicSquare[0][1] = +j_1-j_2+j_3;
//   MagicSquare[0][2] = +j_1+j_2-j_3;
//   MagicSquare[1][0] = j_1-m_1;
//   MagicSquare[1][1] = j_2-m_2;
//   MagicSquare[1][2] = j_3-m_3;
//   MagicSquare[2][0] = j_1+m_1;
//   MagicSquare[2][1] = j_2+m_2;
//   MagicSquare[2][2] = j_3+m_3;

//   // Find the smallest element, and its indices
//   S=MagicSquare[0][0], i_S=0, j_S=0;
//   L=MagicSquare[0][0], i_L=0, j_L=0;
//   for(i=0; i<3; ++i) {
//     for(j=0; j<3; ++j) {
//       if(MagicSquare[i][j]<S) {
//         S=MagicSquare[i][j];
//         i_S=i;
//         j_S=j;
//       }
//       if(MagicSquare[i][j]>L) {
//         L=MagicSquare[i][j];
//         i_L=i;
//         j_L=j;
//       }
//     }
//   }
//   if(i_L==i_S && j_L==j_S) {
//     // The magic square is a constant, presumably because j_1=j_2=j_3
//     // and m_1=m_2=m_3=0.  Just randomly change one index
//     j_L = (j_L+1)%3;
//   }

//   if(i_L!=i_S && j_L!=j_S) {
//     // Find the largest element somewhere on S's row or column.  There
//     // may be other largest elements (equal to the one we find here)
//     // on diagonals from S, but there must be one on S's row or
//     // column.  If some values are equal, let's just randomly choose
//     // whichever other one it is, in case the magic square is
//     // constant.  That way we end up with distinct rows and/or
//     // columns.
//     i_L=i_S, j_L=j_S, L=MagicSquare[i_L][j_L];
//     if(MagicSquare[i_L][(j_S+1)%3]>=L) {
//       j_L=(j_S+1)%3;
//       L=MagicSquare[i_L][j_L];
//     }
//     if(MagicSquare[i_L][(j_S+2)%3]>=L) {
//       j_L=(j_S+2)%3;
//       L=MagicSquare[i_L][j_L];
//     }
//     if(MagicSquare[(i_S+1)%3][j_L]>=L) {
//       i_L=(i_S+1)%3;
//       L=MagicSquare[i_L][j_L];
//     }
//     if(MagicSquare[(i_S+2)%3][j_L]>=L) {
//       i_L=(i_S+2)%3;
//       L=MagicSquare[i_L][j_L];
//     }
//   }

//   // Either i_S==i_L or j_S==j_L, per discussion above Eq. (2.11)
//   if(i_S==i_L) { // S and L are in the same row
//     // Get the other column
//     j_T = ((j_S==0||j_L==0) ? ((j_S==1||j_L==1) ? 2 : 1) : 0);
//     if(MagicSquare[(i_L+1)%3][j_L]<MagicSquare[(i_L+2)%3][j_L]) {
//       B = MagicSquare[(i_L+1)%3][j_L];
//       X = MagicSquare[(i_L+1)%3][j_S];
//       T = MagicSquare[(i_L+2)%3][j_T];
//     } else if(MagicSquare[(i_L+2)%3][j_L]<MagicSquare[(i_L+1)%3][j_L]) {
//       B = MagicSquare[(i_L+2)%3][j_L];
//       X = MagicSquare[(i_L+2)%3][j_S];
//       T = MagicSquare[(i_L+1)%3][j_T];
//     } else {
//       if(MagicSquare[(i_L+1)%3][j_T]<=MagicSquare[(i_L+2)%3][j_T]) {
//         B = MagicSquare[(i_L+1)%3][j_L];
//         X = MagicSquare[(i_L+1)%3][j_S];
//         T = MagicSquare[(i_L+2)%3][j_T];
//       } else {
//         B = MagicSquare[(i_L+2)%3][j_L];
//         X = MagicSquare[(i_L+2)%3][j_S];
//         T = MagicSquare[(i_L+1)%3][j_T];
//       }
//     }
//   } else { // S and L are in the same column
//     // Get the other row
//     i_T = ((i_S==0||i_L==0) ? ((i_S==1||i_L==1) ? 2 : 1) : 0);
//     if(MagicSquare[i_L][(j_L+1)%3]<MagicSquare[i_L][(j_L+2)%3]) {
//       B = MagicSquare[i_L][(j_L+1)%3];
//       X = MagicSquare[i_S][(j_L+1)%3];
//       T = MagicSquare[i_T][(j_L+2)%3];
//     } else if(MagicSquare[i_L][(j_L+2)%3]<MagicSquare[i_L][(j_L+1)%3]) {
//       B = MagicSquare[i_L][(j_L+2)%3];
//       X = MagicSquare[i_S][(j_L+2)%3];
//       T = MagicSquare[i_T][(j_L+1)%3];
//     } else {
//       if(MagicSquare[i_T][(j_L+1)%3]<=MagicSquare[i_T][(j_L+2)%3]) {
//         B = MagicSquare[i_L][(j_L+1)%3];
//         X = MagicSquare[i_S][(j_L+1)%3];
//         T = MagicSquare[i_T][(j_L+2)%3];
//       } else {
//         B = MagicSquare[i_L][(j_L+2)%3];
//         X = MagicSquare[i_S][(j_L+2)%3];
//         T = MagicSquare[i_T][(j_L+1)%3];
//       }
//     }
//   }

//   // Eq. (2.13)
//   // INFOTOCERR << int(L*(24+L*(50+L*(35+L*(10+L))))/120 + X*(6+X*(11+X*(6+X)))/24 + T*(2+T*(3+T))/6 + B*(B+1)/2 + S + 1) << "/" << FactorTable.size() << std::endl;
//   return FactorTable[int(L*(24+L*(50+L*(35+L*(10+L))))/120 + X*(6+X*(11+X*(6+X)))/24 + T*(2+T*(3+T))/6 + B*(B+1)/2 + S + 1)];
// }
