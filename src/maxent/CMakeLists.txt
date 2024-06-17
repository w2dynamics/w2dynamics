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
set(MAXENTSRCNAMES 
MersenneTwister.F90
)
set(MAXENTSRCNAMES ${MAXENTSRCNAMES} PARENT_SCOPE)


# Wrap most of our source files into a small library
  set(MAXENT_src)
  foreach(src ${MAXENTSRCNAMES})
    list(APPEND MAXENT_src "${SRCMAXENT}/${src}")
  endforeach(src)

# Define the executable in terms of the source files
ADD_LIBRARY(MAXENTLIB STATIC ${MAXENT_src})
SET_PROPERTY(TARGET MAXENTLIB PROPERTY POSITION_INDEPENDENT_CODE ON)

set_target_properties(MAXENTLIB PROPERTIES COMPILE_FLAGS "-DLAPACK77_Interface")

  set( _name MAXENT )
  # Define the command to generate the Fortran to Python interface module. The
  # output will be a shared library that can be imported by python.
  set(ADDITIONAL_OPTS -fPIC -DLAPACK77_Interface)
  add_custom_command(
    OUTPUT ${_name}module.c ${_name}-f2pywrappers.f ${_name}-f2pywrappers2.f90
    DEPENDS MaximumEntropy.F90 .f2py_f2cmap
    VERBATIM
    COMMAND "${CMAKE_COMMAND}" -E touch "${_name}module.c" "${_name}-f2pywrappers.f" "${_name}-f2pywrappers2.f90"
    COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
            "${CMAKE_CURRENT_SOURCE_DIR}/MaximumEntropy.F90" -m ${_name} --lower
            --f2cmap "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap")

  python_add_library(${_name} MODULE WITH_SOABI
                     "${CMAKE_CURRENT_BINARY_DIR}/${_name}module.c"
                     "${CMAKE_CURRENT_BINARY_DIR}/${_name}-f2pywrappers.f"
                     "${CMAKE_CURRENT_BINARY_DIR}/${_name}-f2pywrappers2.f90"
                     "${CMAKE_CURRENT_SOURCE_DIR}/MaximumEntropy.F90")
  if (CMAKE_C_COMPILER_ID STREQUAL "Clang" AND CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL 16)
    set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/${_name}module.c" PROPERTIES
      COMPILE_OPTIONS "-Wno-error=incompatible-function-pointer-types")
  elseif (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/${_name}module.c" PROPERTIES
      COMPILE_OPTIONS "-Wno-error=incompatible-function-pointer-types")
  elseif (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
    set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/${_name}module.c" PROPERTIES
      COMPILE_OPTIONS "-Wno-error=incompatible-pointer-types")
  endif()
  target_link_libraries(${_name} PRIVATE fortranobject MAXENTLIB ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  add_custom_command(TARGET ${_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${_name}> ${PROJECT_SOURCE_DIR}/w2dyn/maxent/
    COMMAND_EXPAND_LISTS)

#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
#INSTALL(TARGETS ${CTQMCEXE} RUNTIME DESTINATION bin)
