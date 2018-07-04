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
ADD_LIBRARY(MAXENTLIB OBJECT ${MAXENT_src})
SET_PROPERTY(TARGET MAXENTLIB PROPERTY POSITION_INDEPENDENT_CODE ON)

set_target_properties(MAXENTLIB PROPERTIES COMPILE_FLAGS "-DLAPACK77_Interface")

  # create a string that contains the path of the objects.
  set(MAXENTF2PYOBJS)

  foreach(obj ${MAXENTSRCNAMES})
    list(APPEND MAXENTF2PYOBJS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/MAXENTLIB.dir/${obj}.o)#HACK: only works for UNIX makefiles...cmake has different behaviour for VS.
  endforeach(obj)
set(MAXENTF2PYOBJS ${MAXENTF2PYOBJS} PARENT_SCOPE)
set(EXTERNAL_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(APPLE)
    set(EXTERNAL_LIBRARIES "--link-lapack_opt")
endif()

  set( _name MAXENT )
  # Define the command to generate the Fortran to Python interface module. The
  # output will be a shared library that can be imported by python.
    add_custom_command(OUTPUT "${_name}${F2PY_SUFFIX}"
      COMMAND ${F2PY_EXECUTABLE} --quiet -m ${_name} #-c "${CMAKE_CURRENT_BINARY_DIR}/f2py-${_name}/${_name}.pyf"
              --build-dir "${CMAKE_Fortran_MODULE_DIRECTORY}"
              ${_fcompiler_opts} --opt="-O3 -fPIC -DLAPACK77_Interface" -c ${_inc_opts}
                                 ${EXTERNAL_LIBRARIES}
                                ${CMAKE_THREAD_LIBS_INIT} ${MAXENTF2PYOBJS} ${SRCMAXENT}/MaximumEntropy.F90
      DEPENDS ${SRCMAXENT}/MaximumEntropy.F90
      WORKING_DIRECTORY ${SRCMAXENT}
      COMMENT "[F2PY] Building Fortran to Python interface module ${_name}")
  # Add a custom target <name> to trigger the generation of the python module.
  add_custom_target(${_name} ALL DEPENDS "${_name}${F2PY_SUFFIX}")
  # MAXENT.so is generated in the source directory since f2py doesn't allow to specify the output directory of the object file.
  # Therefore we have to move it
  add_custom_command(TARGET ${_name} POST_BUILD 
                    COMMAND ${CMAKE_COMMAND} -E rename
                    ${SRCMAXENT}/${_name}${F2PY_SUFFIX} ${CMAKE_SOURCE_DIR}/maxent/${_name}${F2PY_SUFFIX})
                    
  add_dependencies(${_name} MAXENTLIB)

#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
#INSTALL(TARGETS ${CTQMCEXE} RUNTIME DESTINATION bin)
