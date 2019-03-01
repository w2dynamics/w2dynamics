########################################
# Set up how to compile the source files
########################################
  # Get a list of the include directories.
  # The f2py --include_paths option, used when generating a signature file,
  # needs a colon-separated list. The f2py -I option, used when compiling
  # the sources, must be repeated for every include directory.
  get_directory_property(_inc_dirs INCLUDE_DIRECTORIES)
#  string(REPLACE ";" ":" _inc_paths "${_inc_dirs}")
  set(_inc_opts)
  foreach(_dir ${_inc_dirs})
    list(APPEND _inc_opts "-I${_dir}")
  endforeach(_dir)

#separate off the names of the source files
set(CTQMCSRCNAMES AngularMomentum.F90
Ausgabe.F90
CompoundIndex.F90
Lanczos.F90
LegendrePoly.F90
MatrixUpdate.F90
Nfft_base.F90
Nfft_worm.F90
Nfft_z.F90
Operator.F90
Parameters.F90
Progress.F90
Signals.F90
SparseMatrix.F90
States.F90
Trace.F90
)

 # Wrap most of our source files into a small library
  set(CTQMC_src)
  foreach(src ${CTQMCSRCNAMES})
    list(APPEND CTQMC_src "${SRCCTQMC}/${src}")
  endforeach(src)
  
# Define the executable in terms of the source files
ADD_LIBRARY(CTQMCLIB OBJECT ${CTQMC_src})
SET_PROPERTY(TARGET CTQMCLIB PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(CTQMCLIB PROPERTIES COMPILE_FLAGS "-DLAPACK77_Interface")
add_dependencies(CTQMCLIB mtrng)
 
  # create a string that contains the path of the objects.
  set(F2PYOBJS)

  foreach(obj ${CTQMCSRCNAMES})
    list(APPEND F2PYOBJS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CTQMCLIB.dir/${obj}.o)#HACK: only works for UNIX makefiles...cmake has different behaviour for VS.
  endforeach(obj)

  set(_name CTQMC)
  set (CPPSTDLIBRARY "-lstdc++")#placeholder if we decide to better parse platforms in the future...
  set(SWITCHES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#  set(EXTERNAL_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${NFFT_LIBRARIES}  ${FFTW_LIBRARIES})
  if(APPLE)
    set(SWITCHES "--link-lapack_opt" "--link-fftw")
  endif()

  # hack to not pass .dylib libraries directly because f2py does not accept them as arguments
  if (FFTW_LIBRARIES)
    string(TOUPPER ${FFTW_LIBRARIES} FFTW_LIBS_UPPER)
    if (${FFTW_LIBS_UPPER} MATCHES "\\.DYLIB$")
      get_filename_component(FFTW_CONTAINING_DIR ${FFTW_LIBRARIES} DIRECTORY)
      list(APPEND SWITCHES "-L${FFTW_CONTAINING_DIR}" "-lfftw3")
      set(FFTW_LIBRARIES "")
    endif ()
  endif ()

  if (NFFT_LIBRARIES)
    string(TOUPPER ${NFFT_LIBRARIES} NFFT_LIBS_UPPER)
    if (${NFFT_LIBS_UPPER} MATCHES "\\.DYLIB$")
      get_filename_component(NFFT_CONTAINING_DIR ${NFFT_LIBRARIES} DIRECTORY)
      list(APPEND SWITCHES "-L${NFFT_CONTAINING_DIR}" "-lnfft3")
      set(NFFT_LIBRARIES "")
    endif ()
  endif ()

  # Define the command to generate the Fortran to Python interface module. The
  # output will be a shared library that can be imported by python.
    add_custom_command(OUTPUT "${_name}${F2PY_SUFFIX}" COMMAND ${F2PY_EXECUTABLE} --quiet -m ${_name} #-c "${CMAKE_CURRENT_BINARY_DIR}/f2py-${_name}/${_name}.pyf"
              --build-dir "${CMAKE_Fortran_MODULE_DIRECTORY}"
              ${_fcompiler_opts} --opt="-O3"
              ${SWITCHES}
              ${LIB}/libmtrng.a
              ${CPPSTDLIBRARY}
              ${CMAKE_THREAD_LIBS_INIT} -c ${_inc_opts} ${F2PYOBJS} ${SRCCTQMC}/CTQMC.F90 ${NFFT_LIBRARIES} ${FFTW_LIBRARIES}
      DEPENDS ${SRCCTQMC}/CTQMC.F90 
      WORKING_DIRECTORY ${SRCCTQMC}
#      WORKING_DIRECTORY ${BIN}/ctqmc
      COMMENT "[F2PY] Building Fortran to Python interface module ${_name}")
  # Add a custom target <name> to trigger the generation of the python module.
IF(USE_NFFT AND NOT NFFT_FOUND)
  add_custom_target(${_name} ALL DEPENDS "${_name}${F2PY_SUFFIX}" nfft_local)
ELSE()
  add_custom_target(${_name} ALL DEPENDS "${_name}${F2PY_SUFFIX}")
ENDIF(USE_NFFT AND NOT NFFT_FOUND)
  add_dependencies(${_name} CTQMCLIB)  
  # CTQMC.so is generated in the source directory since f2py doesn't allow to specify the output directory of the object file.
  # Therefore we have to move it
  add_custom_command(TARGET ${_name} POST_BUILD 
                    COMMAND ${CMAKE_COMMAND} -E copy
                    ${SRCCTQMC}/${_name}${F2PY_SUFFIX} ${CMAKE_SOURCE_DIR}/w2dyn/auxiliaries/${_name}${F2PY_SUFFIX})


IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
#INSTALL(TARGETS ${CTQMCEXE} RUNTIME DESTINATION bin)
