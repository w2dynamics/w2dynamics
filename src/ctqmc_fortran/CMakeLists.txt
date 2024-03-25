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
set(CTQMCSRCNAMES Accumulators.F90
AngularMomentum.F90
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
ADD_LIBRARY(CTQMCLIB STATIC ${CTQMC_src})
SET_PROPERTY(TARGET CTQMCLIB PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(CTQMCLIB PROPERTIES COMPILE_FLAGS "-DLAPACK77_Interface")
target_link_libraries(CTQMCLIB mtrng)
if (USE_NFFT)
  target_link_libraries(CTQMCLIB nfft)
endif (USE_NFFT)

  set(_name CTQMC)
  # Define the command to generate the Fortran to Python interface module. The
  # output will be a shared library that can be imported by python.
  add_custom_command(
    OUTPUT ${_name}module.c ${_name}-f2pywrappers.f ${_name}-f2pywrappers2.f90
    DEPENDS CTQMC.F90 .f2py_f2cmap
    VERBATIM
    COMMAND "${CMAKE_COMMAND}" -E touch "${_name}module.c" "${_name}-f2pywrappers.f" "${_name}-f2pywrappers2.f90"
    COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
            "${CMAKE_CURRENT_SOURCE_DIR}/CTQMC.F90" -m ${_name} --lower
            --f2cmap "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap")

  python_add_library(${_name} MODULE WITH_SOABI
                     "${CMAKE_CURRENT_BINARY_DIR}/${_name}module.c"
                     "${CMAKE_CURRENT_BINARY_DIR}/${_name}-f2pywrappers.f"
                     "${CMAKE_CURRENT_BINARY_DIR}/${_name}-f2pywrappers2.f90"
                     "${CMAKE_CURRENT_SOURCE_DIR}/CTQMC.F90")
  set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/${_name}module.c" PROPERTIES
    COMPILE_OPTIONS $<$<AND:$<C_COMPILER_ID:Clang>,$<VERSION_GREATER_EQUAL:$<C_COMPILER_VERSION>,16>>:-Wno-error=incompatible-function-pointer-types>)
  target_link_libraries(${_name} PRIVATE fortranobject CTQMCLIB ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${FFTW_LIBRARIES} $<$<BOOL:${USE_NFFT}>:nfft>)
  add_custom_command(TARGET ${_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy -t "${PROJECT_SOURCE_DIR}/w2dyn/auxiliaries" $<TARGET_FILE:${_name}>
    COMMAND_EXPAND_LISTS)

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
#INSTALL(TARGETS ${CTQMCEXE} RUNTIME DESTINATION bin)
