#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.command.install import install as Install
from distutils.command.build_ext import build_ext as BuildExt
from setuptools.extension import Extension

import os
import os.path
import shutil

class BuildFromMakefile(BuildExt):
    """ Required for correctly building and copying the whole nonsense """
    def run_command(self, command):
        if os.system(command):
            raise OSError(command + " failed")
    
    def run(self):
        self.run_command("make -C ctqmc_fortran")
        self.run_command("make -C maxent")
        self.extensions = [] # clear the list
        BuildExt.run(self)
    
class InstallHack(Install):
    def run(self):
        Install.run(self)
        shutil.copy("auxiliaries/CTQMC.so", 
                    os.path.join(self.install_lib, "auxiliaries"))

fsources = ["ctqmc_fortran/Makefile",
            "ctqmc_fortran/makedependf90.pl",
            "ctqmc_fortran/AngularMomentum.F90",
            "ctqmc_fortran/CTQMC.F90",
            "ctqmc_fortran/H5_lookup.F90",
            "ctqmc_fortran/Lanczos.F90",
            "ctqmc_fortran/LegendrePoly.F90",
            "ctqmc_fortran/MatrixUpdate.F90",
            "ctqmc_fortran/MersenneTwister.F90",
            "ctqmc_fortran/Nfft_base.F90",
            "ctqmc_fortran/Nfft_z.F90",
            "ctqmc_fortran/Operator.F90",
            "ctqmc_fortran/Parameters.F90",
            "ctqmc_fortran/Progress.F90",
            "ctqmc_fortran/Signals.F90",
            "ctqmc_fortran/SparseMatrix.F90",
            "ctqmc_fortran/States.F90",
            "ctqmc_fortran/Trace.F90",
            "ctqmc_fortran/VertexComponent.F90",
            "maxent/Makefile",
            "maxent/makedependf90.pl",
            "maxent/MaximumEntropy.F90",
            "maxent/MersenneTwister.F90",
            "Parameters.in"
            ]

setup(name='w2dynamics',
      version='0.75',
      
      # Description of the package
      description='Wuerzburg/Vienna strong-coupling impurity solver',
      author='N. Parragh, M. Wallerberger, G. Sangiovanni',
      author_email='nico.parragh@gmail.com',
      keywords='ctqmc dmft impurity strong coupling',
      url='http://www.ifp.tuwien.ac.at',
      classifiers=["Development Status :: 4 - Beta",
                   "Environment :: Console",
                   "Intended Audience :: Science/Research",
                   "License :: Other/Proprietary License",
                   "Natural Language :: English",
                   "Programming Language :: Fortran",
                   "Programming Language :: Python :: 2.6",
                   "Topic :: Scientific/Engineering :: Chemistry",
                   "Topic :: Scientific/Engineering :: Physics",
                   ],
      install_requires=['numpy >=1.0', 
                        'scipy', 
                        'h5py >=1.3', 
                        'configobj'],
      
      # Contents, build and deployment instructions
      packages=['auxiliaries', 'dmft', 'maxent'],
      package_data={'auxiliaries': ['configspec']},
      ext_modules=[Extension('auxiliaries/CTQMC', fsources),
                   Extension('maxent/MaximumEntropy', fsources)],
      scripts=['DMFT.py', 'Maxent.py', 'hgrep', 'cprun', 'cthyb'],
      cmdclass={'build_ext': BuildFromMakefile, 
                'install': InstallHack,
                }
     )
