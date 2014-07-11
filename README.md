SphericalFunctions
==================

C++ library for computing Wigner D matrices, spin-weighted spherical
harmonics, and associated quantities needed for functions on the
sphere.  Python bindings are supplied via SWIG.

These functions are written in terms of quaternions, but have been
overloaded to also take Euler angles as arguments.  The advantages of
using quaternions directly in these functions are:

 * Quaternions are completely free of the singularities present in the
   Euler angle representation of rotations.

 * No conversion to Euler angles is necessary when quaternions are
   used for other purposes (which they generally should be, because
   quaternions are far superior to other representations of rotations
   in almost every way).

 * Conversions from Euler angles are trivial, so there is that type of
   backwards compatibility.


Downloading
===========

To get this code, run

    git clone --recursive https://github.com/MOBle/SphericalFunctions.git

The `recursive` flag must be passed because this project uses
[`Quaternions`](https://github.com/MOBle/Quaternions) as a git
submodule, which also needs to be downloaded.  Note, however, that the
GSL dependency in the `Quaternions` module is not needed by this
module.


Compiling with other software
=============================

Compilation is fairly standard.  The only caveat is that the files
`Quaternions/Quaternions.{c,h}pp` from the submodule need to be
compiled and included in any compilation, respectively.


Installing the python module
============================

Though this code can be included as a library in other code, it can
also be used on its own as a python module.  Just run

    python setup.py install --user

The `--user` flag installs the module to the user's home directory,
which means that no root permissions are necessary.

If the build succeeds, just open an python session and type

    import SphericalHarmonics

In *ipython* (not just ordinary python), you can then see your options
using tab completion by typing

    SphericalHarmonics.

(including the dot) and then hitting Tab.  Help is available on any
documented function in ipython by typing a question mark after the
function name.  For example:

    SphericalHarmonics.WignerDMatrix?

In plain python, essentially the same thing can be achieved by
entering `help(SphericalHarmonics.WignerDMatrix)`.  But if you're
using plain python interactively, you should really give ipython a
try.
