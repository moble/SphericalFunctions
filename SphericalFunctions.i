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
%{const char* const SphericalFunctionsErrors[] = {
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
    "Unknown exception",// "Bad value.",
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
    PyExc_RuntimeError, // PyExc_ValueError, // Bad value
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
  try {
    $action;
  } catch(int i) {
    std::stringstream s;
    if(i>-1 && i<SphericalFunctionsNumberOfErrors) { s << "SphericalFunctions exception: " << SphericalFunctionsErrors[i]; }
    else  { s << "SphericalFunctions: Unknown exception number {" << i << "}"; }
    PyErr_SetString(SphericalFunctionsExceptions[i], s.str().c_str());
    return NULL;
  }
}


/////////////////////////////////////////////////
//// These will be needed by the c++ wrapper ////
/////////////////////////////////////////////////
%{
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
