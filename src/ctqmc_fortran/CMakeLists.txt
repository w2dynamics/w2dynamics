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
  set (CPPSTDLIBRARY "-lstdc++")#placeholder if we decide to better parse platforms in the future...

  # hack to not pass .dylib or .framework library file names directly
  # because f2py does not accept them as arguments
  function(rephrase_dylib TARGET_VAR LIBRARIES_STRING)
    set(RESULT_LIST "")
    foreach(entry ${LIBRARIES_STRING})
      string(TOUPPER "${entry}" upcaseentry)
      if(upcaseentry MATCHES "LIB([A-Za-z0-9_]*)\\.DYLIB$")
	get_filename_component(entry_dir "${entry}" DIRECTORY)
	string(LENGTH "${CMAKE_MATCH_1}" libname_length)
	string(FIND "${upcaseentry}" "${CMAKE_MATCH_1}" libname_begin REVERSE)
	string(SUBSTRING "${entry}" "${libname_begin}" "${libname_length}" libname)
	list(APPEND RESULT_LIST "-L${entry_dir}" "-l${libname}")
      elseif((upcaseentry MATCHES ".*ACCELERATE\\.FRAMEWORK.*") OR (upcaseentry MATCHES ".*VECLIB\\.FRAMEWORK.*"))
	# f2py will use numpy.distutils.system_info
	list(APPEND RESULT_LIST "--link-accelerate")
      else()
	list(APPEND RESULT_LIST "${entry}")
      endif()
    endforeach()
    set("${TARGET_VAR}" "${RESULT_LIST}" PARENT_SCOPE)
  endfunction()

  if (FFTW_LIBRARIES)
    rephrase_dylib(FFTW_LIBRARIES "${FFTW_LIBRARIES}")
  endif ()

  if (NFFT_LIBRARIES)
    rephrase_dylib(NFFT_LIBRARIES "${NFFT_LIBRARIES}")
  endif ()

  if (BLAS_LIBRARIES)
    rephrase_dylib(BLAS_LIBRARIES "${BLAS_LIBRARIES}")
  endif ()

  if (LAPACK_LIBRARIES)
    rephrase_dylib(LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
  endif ()

  # Define the command to generate the Fortran to Python interface module. The
  # output will be a shared library that can be imported by python.
  add_custom_command(
    OUTPUT ${_name}module.c ${_name}-f2pywrappers.f
    DEPENDS CTQMC.F90 .f2py_f2cmap
    VERBATIM
    COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
            "${CMAKE_CURRENT_SOURCE_DIR}/CTQMC.F90" -m ${_name} --lower
            --f2cmap "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap")

  python_add_library(${_name} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${_name}module.c"
                     "${CMAKE_CURRENT_SOURCE_DIR}/CTQMC.F90" WITH_SOABI)
  set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/${_name}module.c" PROPERTIES
    COMPILE_OPTIONS $<$<AND:$<CXX_COMPILER_ID:Clang>,$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,16>>:-Wno-error=incompatible-function-pointer-types>)
  target_link_libraries(${_name} PRIVATE fortranobject)
  target_link_libraries(${_name} PUBLIC CTQMCLIB $<$<BOOL:${USE_NFFT}>:nfft>)

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
#INSTALL(TARGETS ${CTQMCEXE} RUNTIME DESTINATION bin)
