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

  set(_name CTQMC)
  set (CPPSTDLIBRARY "-lstdc++")#placeholder if we decide to better parse platforms in the future...

  # hack to not pass .dylib libraries directly because f2py does not
  # accept them as arguments, while it does accept .so ones...
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
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${_name}${F2PY_SUFFIX}"
    COMMAND ${CMAKE_COMMAND} -E copy "${SRCCTQMC}/.f2py_f2cmap" "${CMAKE_CURRENT_BINARY_DIR}/.f2py_f2cmap"
    COMMAND ${PYTHON_EXECUTABLE} -m numpy.f2py -c $<$<CONFIG:RELEASE>:--quiet> -m ${_name}
      --build-dir "${CMAKE_Fortran_MODULE_DIRECTORY}"
      ${_fcompiler_opts}
      ${LIB}/libmtrng.a
      ${CPPSTDLIBRARY}
      ${CMAKE_THREAD_LIBS_INIT} ${_inc_opts} $<TARGET_LINKER_FILE:CTQMCLIB> ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${NFFT_LIBRARIES} ${FFTW_LIBRARIES} ${SRCCTQMC}/CTQMC.F90
    DEPENDS "${SRCCTQMC}/CTQMC.F90" "${SRCCTQMC}/.f2py_f2cmap" CTQMCLIB
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "[F2PY] Building Fortran to Python interface module ${_name}")
  # Add a custom target <name> to trigger the generation of the python module.
IF(USE_NFFT AND NOT NFFT_FOUND)
  add_custom_target(${_name} ALL DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${_name}${F2PY_SUFFIX}" nfft_local)
ELSE()
  add_custom_target(${_name} ALL DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${_name}${F2PY_SUFFIX}")
ENDIF(USE_NFFT AND NOT NFFT_FOUND)
  add_dependencies(${_name} CTQMCLIB)
  # CTQMC.so is generated in the source directory since f2py doesn't allow to specify the output directory of the object file.
  # Therefore we have to move it
  add_custom_command(TARGET ${_name} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy
                    "${CMAKE_CURRENT_BINARY_DIR}/${_name}${F2PY_SUFFIX}" "${CMAKE_SOURCE_DIR}/w2dyn/auxiliaries/${_name}${F2PY_SUFFIX}"
		    COMMAND ${CMAKE_COMMAND} -E remove
                    "${CMAKE_CURRENT_BINARY_DIR}/${_name}${F2PY_SUFFIX}")


IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
#INSTALL(TARGETS ${CTQMCEXE} RUNTIME DESTINATION bin)
