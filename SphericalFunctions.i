// -*- c++ -*-

// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

%module SphericalFunctions

 // Quiet warnings about overloaded operators being ignored.
#pragma SWIG nowarn=362,389,401,509
%include <typemaps.i>
%include <stl.i>

#ifndef SWIGIMPORTED
// Use numpy below
%{
  #define SWIG_FILE_WITH_INIT
  %}
%include "Quaternions/numpy.i"
%init %{
  import_array();
%}
%pythoncode %{
  import numpy;
%}
#endif

%import "Quaternions/Quaternions.i"
%import "Quaternions/Quaternions_typemaps.i"

%include "docs/SphericalFunctions_Doc.i"

///////////////////////////////////
//// Handle exceptions cleanly ////
///////////////////////////////////

// The following will appear in the header of the `_wrap.cpp` file.
%{
  // The following allows us to elegantly fail in python from manual
  // interrupts and floating-point exceptions found in the c++ code.
  // The setjmp part of this was inspired by the post
  // <http://stackoverflow.com/a/12155582/1194883>.  The code for
  // setting the csr flags is taken from SpEC.
  #include <csetjmp>
  #include <csignal>
  #ifdef __APPLE__
    #include <xmmintrin.h>
    int fegetexcept() { return _mm_getcsr(); }
  #else
    #include <fenv.h>     // For feenableexcept. Doesn't seem to be a <cfenv>
    #ifndef __USE_GNU
      extern "C" int feenableexcept (int EXCEPTS);
    #endif
  #endif
  namespace SphericalFunctions {
    static sigjmp_buf FloatingPointExceptionJumpBuffer;
    void FloatingPointExceptionHandler(int sig) {
      siglongjmp(FloatingPointExceptionJumpBuffer, sig);
    }
    class ExceptionHandlerSwitcher {
    private:
      int OriginalExceptionFlags;
      void (*OriginalFloatingPointExceptionHandler)(int);
    public:
      ExceptionHandlerSwitcher()
        : OriginalExceptionFlags(fegetexcept()),
          OriginalFloatingPointExceptionHandler(signal(SIGFPE, FloatingPointExceptionHandler))
      {
        #ifdef __APPLE__
          _mm_setcsr( _MM_MASK_MASK &~
                     (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO));
        #else
          feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
        #endif
      }
      ~ExceptionHandlerSwitcher() {
        #ifdef __APPLE__
          _mm_setcsr(OriginalExceptionFlags);
        #else
          feenableexcept(OriginalExceptionFlags);
        #endif
        signal(SIGFPE, OriginalFloatingPointExceptionHandler);
      }
    };
  }

  const char* const SphericalFunctionsErrors[] = {
    "This function is not yet implemented.",
    "Unknown exception",// "Failed system call.",
    "Unknown exception",// "Bad file name.",
    "Unknown exception",// "Failed GSL call.",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Bad value.",
    "Unknown exception",// "Bad switches; we should not have gotten here.",
    "Index out of bounds.",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",// "Vector size mismatch.",
    "Unknown exception",// "Matrix size mismatch.",
    "Unknown exception",// "Matrix size is assumed to be 3x3 in this function.",
    "Unknown exception",// "Not enough points to take a derivative.",
    "Unknown exception",// "Empty intersection requested.",
  };
  const int SphericalFunctionsNumberOfErrors = 20;
  PyObject* const SphericalFunctionsExceptions[] = {
    PyExc_NotImplementedError, // Not implemented
    PyExc_RuntimeError, // PyExc_SystemError, // Failed system call
    PyExc_RuntimeError, // PyExc_IOError, // Bad file name
    PyExc_RuntimeError, // PyExc_RuntimeError, // GSL failed
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_ValueError, // Bad value
    PyExc_RuntimeError, // PyExc_ValueError, // Bad switches
    PyExc_IndexError, // Index out of bounds
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // PyExc_AssertionError, // Mismatched vector size
    PyExc_RuntimeError, // PyExc_AssertionError, // Mismatched matrix size
    PyExc_RuntimeError, // PyExc_AssertionError, // 3x3 matrix assumed
    PyExc_RuntimeError, // PyExc_AssertionError, // Not enough points for derivative
    PyExc_RuntimeError, // PyExc_AssertionError, // Empty intersection
  };
%}

// This will go inside every python wrapper for any function I've
// included; the code of the function itself will replace `$action`.
// It's a good idea to try to keep this part brief, just to cut down
// the size of the wrapper file.
%exception {
  if (!sigsetjmp(SphericalFunctions::FloatingPointExceptionJumpBuffer, 1)) {
    try {
      // const SphericalFunctions::ExceptionHandlerSwitcher Switcher;
      $action;
    } catch(int i) {
      std::stringstream s;
      if(i>-1 && i<SphericalFunctionsNumberOfErrors) { s << "$fulldecl: " << SphericalFunctionsErrors[i]; }
      else  { s << "$fulldecl: Unknown exception number {" << i << "}"; }
      PyErr_SetString(SphericalFunctionsExceptions[i], s.str().c_str());
      return 0;
    } catch(...) {
      PyErr_SetString(PyExc_RuntimeError, "$fulldecl: Unknown exception; default handler");
      return 0;
    }
  } else {
    PyErr_SetString(PyExc_RuntimeError, "$fulldecl: Caught a floating-point exception in the c++ code.");
    return 0;
  }
}


/////////////////////////////////////////////////
//// These will be needed by the c++ wrapper ////
/////////////////////////////////////////////////
%{
  #include <vector>
  #include <iostream>
  #include <string>
  #include <sstream>
  #include <iomanip>
  #include <complex>
  #include "Quaternions.hpp"
  #include "Combinatorics.hpp"
  #include "WignerDMatrices.hpp"
  #include "SWSHs.hpp"
%}


%pythoncode %{
  ## We must be able to import numpy
  import numpy

  ## We must be able to import Quaternions
  import Quaternions

  ## We might be able to get away without spinsfast
  try :
    import spinsfast
  except ImportError :
    pass
%}

//////////////////////////////////////////////////////////////////////
//// The following translates between c++ and python types nicely ////
//////////////////////////////////////////////////////////////////////
//// Make sure std::strings are dealt with appropriately
%include <std_string.i>
//// Make sure std::complex numbers are dealt with appropriately
%include <std_complex.i>
//// Make sure std::vectors are dealt with appropriately
%include <std_vector.i>
namespace std {
  // %template(complexd) complex<double>; // Don't use this line!!!
  %template(vectori) vector<int>;
  %template(vectorvectori) vector<vector<int> >;
  %template(vectorc) vector<std::complex<double> >;
  %template(vectorvectorc) vector<vector<std::complex<double> > >;
  %template(vectorq) vector<Quaternions::Quaternion>;
  %template(vectors) vector<string>;
  %template(vectorvectors) vector<vector<std::string> >;
};


////////////////////////////////////////////////////////////
//// Import the various functions for spherical harmonics //
////////////////////////////////////////////////////////////
%include "Combinatorics.hpp"
%include "WignerDMatrices.hpp"
%include "SWSHs.hpp"


/// Add utility functions that are specific to python.  Note that
/// these are defined in the SphericalFunctions namespace.
%insert("python") %{


%}
