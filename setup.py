#! /usr/bin/env python

"""
Installation file for python code associated with the paper "Angular
velocity of gravitational radiation from precessing binaries and the
corotating frame".

To build this code and run it in place, run
    python setup.py build_ext --inplace
then open python and type 'import SphericalFunctions' in
the current directory.

To install in the user's directory, run
    python setup.py install --user
Now, 'import SphericalFunctions' may be run from a python
instance started in any directory on the system.
"""

from os.path import isdir, isfile
from subprocess import check_call
from sys import argv, executable

# Make sure `Quaternions` is either checked out here or at the level above
if isdir('Quaternions') and isfile('Quaternions/numpy.i'):
    QuaternionsPath = 'Quaternions'
elif isdir('../Quaternions') and isfile('../Quaternions/numpy.i'):
    QuaternionsPath = '../Quaternions'
else:
    raise EnvironmentError("Can't find `Quaternions` module.  Did you forget to `git submodule init` and `git submodule update`?")

## Build Quaternions first
print("\nInstalling Quaternions")
cmd = ' '.join(['cd {0} && {1}'.format(QuaternionsPath, executable),]+argv)
print(cmd)
check_call(cmd, shell=True)
print("Finished installing Quaternions\n")

## Check for `--no-GSL` option
if '--no-GSL' in argv:
    GSL=False
    GSLDef = ''
    argv.remove('--no-GSL')
else:
    GSL=True
    GSLDef = '-DUSE_GSL'

## If PRD won't let me keep a subdirectory, make one
from os.path import exists
from os import makedirs
if not exists('SphericalFunctions') :
    makedirs('SphericalFunctions')

## distutils doesn't build swig modules in the correct order by
## default -- the python module is installed first.  This will pop
## 'build_ext' to the beginning of the command list.
from distutils.command.build import build
build.sub_commands = sorted(build.sub_commands, key=lambda sub_command: int(sub_command[0]!='build_ext'))

## We also need to copy the SWIG-generated python script SphericalFunctions.py
## to SphericalFunctions/__init__.py so that it gets installed correctly.
from distutils.command.build_ext import build_ext as _build_ext
from distutils.file_util import copy_file
class build_ext(_build_ext):
    """Specialized Python source builder for moving SWIG module."""
    def run(self):
        _build_ext.run(self)
        copy_file('SphericalFunctions.py', 'SphericalFunctions/__init__.py')

## Now import the basics
from distutils.core import setup, Extension
from subprocess import check_output, CalledProcessError
from os import devnull, environ

# Add directories for numpy and other inclusions
from numpy import get_include
IncDirs = [get_include(), QuaternionsPath]
LibDirs = []

# If /opt/local directories exist, use them
if isdir('/opt/local/include'):
    IncDirs += ['/opt/local/include']
if isdir('/opt/local/lib'):
    LibDirs += ['/opt/local/lib']

# Add directories for GSL, if needed
if GSL :
    SourceFiles = [QuaternionsPath+'/Quaternions.cpp',
                   QuaternionsPath+'/QuaternionUtilities.cpp',
                   QuaternionsPath+'/IntegrateAngularVelocity.cpp',
                   'Combinatorics.cpp',
                   'WignerDMatrices.cpp',
                   'SWSHs.cpp',
                   'SphericalFunctions.i']
    Dependencies = [QuaternionsPath+'/Quaternions.hpp',
                    QuaternionsPath+'/QuaternionUtilities.hpp',
                    QuaternionsPath+'/IntegrateAngularVelocity.hpp',
                    'Combinatorics.hpp',
                    'WignerDMatrices.hpp',
                    'SWSHs.hpp',
                    'Errors.hpp']
    Libraries = ['gsl', 'gslcblas']
    ## See if GSL_HOME is set; if so, use it
    if "GSL_HOME" in environ :
        IncDirs = [environ["GSL_HOME"]+'/include'] + IncDirs
        LibDirs = [environ["GSL_HOME"]+'/lib'] + IncDirs
else :
    SourceFiles = [QuaternionsPath+'/Quaternions.cpp',
                   QuaternionsPath+'/Utilities.cpp',
                   'Combinatorics.cpp',
                   'WignerDMatrices.cpp',
                   'SWSHs.cpp',
                   'SphericalFunctions.i']
    Dependencies = [QuaternionsPath+'/Quaternions.hpp',
                    QuaternionsPath+'/Utilities.hpp',
                    'Combinatorics.hpp',
                    'WignerDMatrices.hpp',
                    'SWSHs.hpp',
                    'Errors.hpp']
    Libraries = []

## See if FFTW3_HOME is set; if so, use it
if "FFTW3_HOME" in environ :
    IncDirs += [environ["FFTW3_HOME"]+'/include']
    LibDirs += [environ["FFTW3_HOME"]+'/lib']

## Remove a compiler flag that doesn't belong there for C++
import distutils.sysconfig as ds
cfs=ds.get_config_vars()
for key, value in cfs.items() :
    if(type(cfs[key])==str) :
        cfs[key] = value.replace('-Wstrict-prototypes', '').replace('-Wshorten-64-to-32', '')

## Read in the license
try :
    with open('LICENSE', 'r') as myfile :
        License=myfile.read()
except IOError :
    License = 'See LICENSE file in the source code for details.'

## Add -py3 if this is python3
swig_opts=['-globals', 'constants', '-c++', '-builtin']
if '../' in QuaternionsPath:
    swig_opts+=['-I..',]
try:
    import sys
    python_major = sys.version_info.major
    if(python_major==3) :
        swig_opts += ['-py3']
except AttributeError:
    pass # This should probably be an error, because python is really old, but let's keep trying...

## This does the actual work
setup(name="SphericalFunctions",
      # version=PackageVersion,
      description='Library implementing Wigner D matrices and spin-weighted spherical harmonics in C++, with python bindings via SWIG',
      #long_description=""" """,
      author='Michael Boyle',
      author_email='boyle@astro.cornell.edu',
      url='https://github.com/MOBle',
      license=License,
      packages = ['SphericalFunctions'],
      # py_modules = ['SphericalFunctions'],
      # scripts = ['Scripts/RunExtrapolations.py', 'Scripts/ConvertGWDatToH5.py'],
      ext_modules = [
        Extension('_SphericalFunctions',
                  sources=SourceFiles,
                  depends=Dependencies,
                  include_dirs=IncDirs,
                  library_dirs=LibDirs,
                  libraries=Libraries,
                  #define_macros = [('CodeRevision', CodeRevision)],
                  language='c++',
                  swig_opts=swig_opts,
                  extra_link_args=['-fPIC'],
                  extra_compile_args=['-Wno-deprecated', '-ffast-math', '-O3', GSLDef],
                  # extra_compile_args=['-fopenmp']
              ),
      ],
      # classifiers = ,
      # distclass = ,
      # script_name = ,
      # script_args = ,
      # options = ,
      # license = ,
      # keywords = ,
      # platforms = ,
      # cmdclass = ,
      cmdclass={'build_ext': build_ext},
      # data_files = ,
      # package_dir =
      )
