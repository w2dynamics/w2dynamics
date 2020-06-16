#Check wether we have already searched the Python interpreter
if(NOT PYTHONINTERP_FOUND)
    find_package(PythonInterp REQUIRED)
endif()

if (NOT F2PY_SUFFIX)
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import distutils.sysconfig; print(distutils.sysconfig.get_config_var('EXT_SUFFIX') or distutils.sysconfig.get_config_var('SO'))"
                  OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
                  RESULT_VARIABLE FOUND_PYTHON_EXT_SUFFIX)
  string(STRIP "${PYTHON_EXT_SUFFIX}" PYTHON_EXT_SUFFIX)
  if (NOT ${FOUND_PYTHON_EXT_SUFFIX} EQUAL 0)
    set (F2PY_SUFFIX "" CACHE STRING "Suffix added by F2Py to the module name to get the output file name." )
    message(FATAL_ERROR "Unable to determine file extension of compiled Python modules - specify it with F2PY_SUFFIX")
  endif (NOT ${FOUND_PYTHON_EXT_SUFFIX} EQUAL 0)
  set (F2PY_SUFFIX ${PYTHON_EXT_SUFFIX} CACHE INTERNAL "the F2PY extension")
endif (NOT F2PY_SUFFIX)
#strip trailing newlines
string(REGEX REPLACE "\n$" "" F2PY_SUFFIX "${F2PY_SUFFIX}")

## Path to the f2py executable
find_program(F2PY_EXECUTABLE NAMES "f2py${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}"
                                   "f2py-${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}"
                                   "f2py"  # if a user-installed distribution comes with 'f2py' only, prefer it to potentially incompatible system-wide 'f2pyx'
                                   "f2py${PYTHON_VERSION_MAJOR}"
                             PATHS ~/.local/bin
                             REQUIRED)

# Get the compiler-id and map it to compiler vendor as used by f2py.
  # Currently, we only check for GNU, but this can easily be extended. 
  # Cache the result, so that we only need to check once.
  if(NOT F2PY_FCOMPILER)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
        set(_fcompiler "gnu95")
      else(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
        set(_fcompiler "gnu")
      endif(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
      
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      # using Intel
      # A small check to distinguish 32bits from 64bits
      if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(_fcompiler "intelem")#We don't support Intel for Itanium
      else(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(_fcompiler "intel")
      endif(CMAKE_SIZEOF_VOID_P EQUAL 8)
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "XL")
        set(_fcompiler "ibm")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "PathScale")
        set(_fcompiler "pathf95")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "SunPro")
        set(_fcompiler "sun")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Absoft")
        set(_fcompiler "absoft")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MIPSPro")
        set(_fcompiler "mips")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "HP-UX")
        set(_fcompiler "hpux")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
        set(_fcompiler "pg")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "G95")
        set(_fcompiler "g95")
    elseif()
        set(_fcompiler "F2PY_FCOMPILER-NOTFOUND")
    endif()


    set(F2PY_FCOMPILER ${_fcompiler} CACHE STRING
      "F2PY: Fortran compiler type by vendor" FORCE)
    if(NOT F2PY_FCOMPILER)
      message(STATUS "[F2PY]: Could not determine Fortran compiler type. "
                     "Troubles ahead!")
    endif(NOT F2PY_FCOMPILER)
  endif(NOT F2PY_FCOMPILER)

   # Set f2py compiler options: compiler vendor and path to Fortran77/90 compiler.
  if(F2PY_FCOMPILER)
    string(REPLACE " " ";" CMAKE_Fortran_FLAGS_LST ${CMAKE_Fortran_FLAGS})
    set(_fcompiler_opts --fcompiler=${F2PY_FCOMPILER})
    list(APPEND _fcompiler_opts --f77exec=${CMAKE_Fortran_COMPILER} --f77flags="${CMAKE_Fortran_FLAGS_LST}" -DF2PY)
    if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
      list(APPEND _fcompiler_opts --f90exec=${CMAKE_Fortran_COMPILER} --f90flags="${CMAKE_Fortran_FLAGS_LST}" -DF2PY)
    endif(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    if(USE_OPENMP)
      list(APPEND CMAKE_Fortran_FLAGS -DF2PY -fopenmp)
      list(APPEND _fcompiler_opts --f90flags="${CMAKE_Fortran_FLAGS_LST}")
    endif(USE_OPENMP)
endif(F2PY_FCOMPILER)
