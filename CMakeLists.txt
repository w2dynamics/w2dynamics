# CMake project file for ctqmc_fortran

##################################################
# Define the project and the dependencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(ctqmc Fortran CXX C)

# Set the version
SET(VERSION 1.0)

# Set the module path, so that CMake also considers our modules
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
OPTION(USE_MPI "Use the MPI library for parallelization" ON)
#OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

SET( CMAKE_Fortran_FLAGS             "${CMAKE_Fortran_FLAGS} -DLAPACK77_Interface" )

############################################################
# Define the actual files and folders that make up the build
############################################################

find_package(PythonInterp 2.4)# Systems with python that old won't have up-to-date cmake....
find_package(PackageHandleStandardArgs)
if(NOT DEFINED ${PYTHON_VERSION_MAJOR})
#cmake-2.4 is so old that the python module does not set the Python versions... we have to 
# find them ourselves
#message(STATUS "trying to guess python")
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%s.%s' % sys.version_info[0:2])"
 RESULT_VARIABLE PYVER)
endif()


INCLUDE(${CMAKE_MODULE_PATH}/installviapip.cmake)

option(USE_NFFT "Compile with support for routines requiring the NFFT library for non-equispaced FFTs" ON)

if (USE_NFFT)
  find_package(NFFT "3.3.0")
  if(NOT NFFT_FOUND)
    find_package(FFTW REQUIRED)
    include(ExternalProject)
#   set(ExternalProjectCMakeArgs
#     -DHDF5_BUILD_CPP_LIB=ON
#     )
    set(_nfft_version "3.4.1")
#   file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/nfft_local/src)
    ExternalProject_Add(nfft_local
      URL http://www-user.tu-chemnitz.de/~potts/nfft/download/nfft-${_nfft_version}.tar.gz
      PREFIX nfft_local
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --libdir=<INSTALL_DIR>/lib
      BUILD_COMMAND make
      INSTALL_COMMAND make install
      )
# in the previous findPackage(NFFT) call an old version of NFFT could have been found. Here we override it with our up-to-date self-compiled version. We don't use the include directories.
    set(NFFT_LIBRARIES ${PROJECT_BINARY_DIR}/nfft_local/lib/libnfft3.a) 
  endif(NOT NFFT_FOUND)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DUSE_NFFT=1")
else ()
  set(NFFT_LIBRARIES "")
endif ()

find_package(NUMPY 1.4)
IF(NOT NUMPY_FOUND)
    MESSAGE(FATAL_ERROR "NumPy NOT found in a recent enough version!")
ENDIF()
find_package(HDF5 1.6 COMPONENTS Fortran)# enforce the availability of the Fortran bindings. FIXME: Test on a suitable computer

if(NOT HDF5_FOUND)
include(ExternalProject)
    set(hdf5_local_major_version "1.8")
    set(hdf5_local_minor_version "20")
    find_package(ZLIB)
    #HDF5 also requires in its configure step a C++ compiler
    ExternalProject_Add(HDF5_LOCAL
        URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_local_major_version}/hdf5-${hdf5_local_major_version}.${hdf5_local_minor_version}/src/hdf5-${hdf5_local_major_version}.${hdf5_local_minor_version}.tar.gz
        PREFIX hdf5_local
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-fortran
        BUILD_COMMAND make
        INSTALL_COMMAND make install
    )
    set(HDF5_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/hdf5_local/include)
    set(HDF5_Fortran_LIBRARIES ${PROJECT_BINARY_DIR}/hdf5_local/lib/libhdf5_fortran.so;${ZLIB_LIBRARIES})
    set(HDF5_LIBRARY_DIRS ${PROJECT_BINARY_DIR}/hdf5_local/lib/)
    set(HDF5_Fortran_COMPILER_EXECUTABLE ${PROJECT_BINARY_DIR}/hdf5_local/bin/h5fc)
    #Try again to find HDF5
    set(ENV{HDF5_ROOT} ${PROJECT_BINARY_DIR}/hdf5_local)
    find_package(HDF5 1.6 COMPONENTS Fortran)
endif()

find_package(MPI4PY)
find_package(H5PY)
find_package(SCIPY 0.10)
IF(NOT SCIPY_FOUND)
    MESSAGE(FATAL_ERROR "SciPy NOT found in a recent enough version!")
ENDIF()
find_package(ConfigObj)
if(NOT H5PY_FOUND)
    install_via_pip("h5py" ERRORCODE)
    if (${ERRORCODE})
    MESSAGE(FATAL_ERROR "Couldn't install H5PY -> ABORTING!!!")
    endif()
endif()
if(NOT MPI4PY_FOUND AND MPI_Fortran_FOUND)
    install_via_pip("mpi4py" ERRORCODE)
    if (${ERRORCODE})
    MESSAGE(FATAL_ERROR "Couldn't install MPI4PY -> ABORTING!!!")
    endif()
endif()
if(NOT CONFIGOBJ_FOUND)
    install_via_pip("configobj" ERRORCODE)
    if (${ERRORCODE})
    MESSAGE(FATAL_ERROR "Couldn't install configobj(required for DMFT.py) -> ABORTING!!!")
    endif()
endif()
#set up some common F2PY Variables
INCLUDE(${CMAKE_MODULE_PATH}/SetupF2Py.cmake)

# Define the executable name
SET(CTQMCEXE ctqmc)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB ${CMAKE_SOURCE_DIR}/lib)
#SET(BIN ${CMAKE_SOURCE_DIR}/bin)

SET(SRCCTQMC ${SRC}/ctqmc_fortran)
SET(SRCMAXENT ${SRC}/maxent)
SET(SRCMTRNG ${SRC}/mtrng)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# The source for the BAR library and have it placed in the lib folder
ADD_SUBDIRECTORY(${SRCMTRNG} ${LIB})

# The source for the ctqmc binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRCCTQMC})
ADD_SUBDIRECTORY(${SRCMAXENT})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    ADD_SUBDIRECTORY(testsuite/mtrng.tests)
    ADD_SUBDIRECTORY(testsuite/maxent.tests)
    ADD_SUBDIRECTORY(testsuite/ctqmc.tests)
    enable_testing()
endif()
