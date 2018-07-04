# - FindMPI4PY
# By Florian Goth, fgoth@physik.uni-wuerzburg.de
# Find mpi4Py includes and library
# This module defines:
# MPI4PY_FOUND
# MPI4PY_VERSION             - the version of mpi4Py found as a string
# MPI4PY_VERSION_MAJOR       - the major version number of mpi4Py
# MPI4PY_VERSION_MINOR       - the minor version number of mpi4Py
# MPI4PY_VERSION_DECIMAL     - e.g. version 1.6.1 is 10601
# MPI4PY_INCLUDE_DIRS        - path to the mpi4Py include files

# Finding mpi4Py involves calling the Python interpreter
#Check wether we have already searched the Python interpreter
if(NOT PYTHONINTERP_FOUND)
    find_package(PythonInterp REQUIRED)
endif()

if(NOT MPI4PY_FOUND AND PYTHONINTERP_FOUND)
    execute_process(COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "exec(\"try:\\n import mpi4py;\\n print(mpi4py.__version__);\\n print(mpi4py.get_include())\\nexcept:\\n exit(1)\")"
      OUTPUT_VARIABLE _MPI4PY_VALUES
      RESULT_VARIABLE MPI4PY_COMMAND_RESULT
      OUTPUT_STRIP_TRAILING_WHITESPACE)

#    if(NOT MPI4PY_COMMAND_RESULT MATCHES 0)
#        message("SciPy import failure:\n${_MPI4PY_ERROR_VALUE}")
    if(MPI4PY_COMMAND_RESULT MATCHES 0)
        # Convert the process output into a list
        string(REGEX REPLACE ";" "\\\\;" _MPI4PY_VALUES ${_MPI4PY_VALUES})
        string(REGEX REPLACE "\n" ";" _MPI4PY_VALUES ${_MPI4PY_VALUES})
        list(GET _MPI4PY_VALUES 0 MPI4PY_VERSION)
        list(GET _MPI4PY_VALUES 1 MPI4PY_INCLUDE_DIRS)

        # Make sure all directory separators are '/'
        string(REGEX REPLACE "\\\\" "/" MPI4PY_INCLUDE_DIRS ${MPI4PY_INCLUDE_DIRS})

        # Get the major and minor version numbers
        string(REGEX REPLACE "\\." ";" _MPI4PY_VERSION_LIST ${MPI4PY_VERSION})
        list(GET _MPI4PY_VERSION_LIST 0 MPI4PY_VERSION_MAJOR)
        list(GET _MPI4PY_VERSION_LIST 1 MPI4PY_VERSION_MINOR)
#        list(GET _MPI4PY_VERSION_LIST 2 MPI4PY_VERSION_PATCH)
        #string(REGEX MATCH "[0-9]*" NUMPY_VERSION_PATCH ${NUMPY_VERSION_PATCH})
        math(EXPR MPI4PY_VERSION_DECIMAL
            "(${MPI4PY_VERSION_MAJOR} * 10000) + (${MPI4PY_VERSION_MINOR} * 100)")
        
    endif(MPI4PY_COMMAND_RESULT MATCHES 0)
endif(NOT MPI4PY_FOUND AND PYTHONINTERP_FOUND)

find_package_handle_standard_args(  MPI4PY
                                    REQUIRED_VARS MPI4PY_INCLUDE_DIRS
                                    VERSION_VAR MPI4PY_VERSION
                                 )
