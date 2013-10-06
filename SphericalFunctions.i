// -*- c++ -*-

// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details


%module SphericalFunctions

 // Quiet warnings about overloaded operators being ignored.
#pragma SWIG nowarn=362,389,401,509
%include <typemaps.i>
%include <stl.i>

// %{
//   #include "/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h"
//   #include <vector>
// %}

// %init %{
//   import_array();
// %}

// %typemap(out) std::vector<int> {
//   npy_intp result_size = $1.size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_INT);
//   int* dat = (int*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = $1[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<int>& {
//   npy_intp result_size = $1->size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_INT);
//   int* dat = (int*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = (*$1)[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<std::vector<int> >& {
//   npy_intp result_size = $1->size();
//   npy_intp result_size2 = (result_size>0 ? (*$1)[0].size() : 0);
//   npy_intp dims[2] = { result_size, result_size2 };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_INT);
//   int* dat = (int*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { for (size_t j = 0; j < result_size2; ++j) { dat[i*result_size2+j] = (*$1)[i][j]; } }
//   $result = PyArray_Return(npy_arr);
// }

// %typemap(out) std::vector<double> {
//   npy_intp result_size = $1.size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
//   double* dat = (double*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = $1[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<double>& {
//   npy_intp result_size = $1->size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
//   double* dat = (double*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = (*$1)[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<std::vector<double> >& {
//   npy_intp result_size = $1->size();
//   npy_intp result_size2 = (result_size>0 ? (*$1)[0].size() : 0);
//   npy_intp dims[2] = { result_size, result_size2 };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
//   double* dat = (double*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { for (size_t j = 0; j < result_size2; ++j) { dat[i*result_size2+j] = (*$1)[i][j]; } }
//   $result = PyArray_Return(npy_arr);
// }

%include "docs/SphericalFunctions_Doc.i"

///////////////////////////////////
//// Handle exceptions cleanly ////
///////////////////////////////////
%exception {
  try {
    $action;
  } catch(int i) {
    if(i==0) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Not yet implemented.");
    } else if(i==1) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Index out of bounds.");
    } else if(i==2) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Infinitely many solutions.");
    } else if(i==3) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Input vector size mismatch.");
    } else if(i==4) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Cannot extrapolate quaternions.");
    } else if(i==5) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Matrix size mismatch.");
    } else if(i==6) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Matrix size is assumed to be 3x3 in this function.");
    } else if(i==7) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Quaternion constructor's vector size not understood; should be 3 or 4.");
    } else if(i==8) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Waveform is missing requested l,m component.");
    } else if(i==9) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Bad file name.");
    } else if(i==10) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Not enough points to take a derivative.");
    } else if(i==11) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Empty intersection requested.");
    } else if(i==12) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Failed system call.");
    } else if(i==13) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Wrong FrameType for this operation.  Maybe you forgot to `SetFrameType`?");
    } else if(i==14) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: GSL failed.");
    } else if(i==15) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Bad Waveform information.");
    } else if(i==16) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Bad switches; we should not have gotten here.");
    } else if(i==17) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Bad value.");
    } else  {
      PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
    }
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
  #include "SWSHs.hpp"
%}


%pythoncode %{
  try :
    import Quaternions
  except ImportError :
    pass
%}


// This imports spinsfast.so, so we don't need to pull any tricks
%pythoncode %{
  try :
    import spinsfast
  except ImportError :
    pass
  %}

//////////////////////////////////////////////////////////////////////
//// The following translates between c++ and python types nicely ////
//////////////////////////////////////////////////////////////////////
//// This lets me use numpy.array in the code below
%pythoncode %{
  import numpy;
  %}
//// Make sure std::strings are dealt with appropriately
%include <std_string.i>
//// Make sure std::complex numbers are dealt with appropriately
%include <std_complex.i>
// namespace std {
//   %template(complexd) complex<double>; // Don't use this line!!!
// };
//// Make sure std::vectors are dealt with appropriately
%include <std_vector.i>
namespace Quaternions {
  class Quaternion;
 };
namespace std {
  %template(vectori) vector<int>;
  %template(vectorvectori) vector<vector<int> >;
  %template(vectord) vector<double>;
  %template(vectorvectord) vector<vector<double> >;
  %template(vectorc) vector<std::complex<double> >;
  %template(vectorvectorc) vector<vector<std::complex<double> > >;
  %template(vectorq) vector<Quaternions::Quaternion>;
  %template(vectors) vector<string>;
  %template(vectorvectors) vector<vector<std::string> >;
};


////////////////////////////////////////////////////////////
//// Import the various functions for spherical harmonics //
////////////////////////////////////////////////////////////
%include "SWSHs.hpp"


/// Add utility functions that are specific to python.  Note that
/// these are defined in the SphericalFunctions namespace.
%insert("python") %{


%}
