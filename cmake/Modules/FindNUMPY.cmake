# - FindNumPy
# By Florian Goth, fgoth@physik.uni-wuerzburg.de
# Find NumPy includes and library
# This module defines:
# NUMPY_FOUND
# NUMPY_VERSION             - the version of NumPy found as a string
# NUMPY_VERSION_MAJOR       - the major version number of NumPy
# NUMPY_VERSION_MINOR       - the minor version number of NumPy
# NUMPY_VERSION_PATCH
# NUMPY_VERSION_DECIMAL     - e.g. version 1.6.1 is 10601
# NUMPY_INCLUDE_DIRS        - path to the SciPy include files

# Finding NumPy involves calling the Python interpreter
#Check wether we have already searched the Python interpreter
if(NOT PYTHONINTERP_FOUND)
    find_package(PythonInterp REQUIRED)
endif()

if(NOT NUMPY_FOUND AND PYTHONINTERP_FOUND)
    execute_process(COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "exec(\"try:\\n import numpy;\\n print(numpy.__version__);\\n print(numpy.get_include())\\nexcept:\\n exit(1)\")"
      OUTPUT_VARIABLE _NUMPY_VALUES
      RESULT_VARIABLE NUMPY_COMMAND_RESULT
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NUMPY_COMMAND_RESULT MATCHES 0)
        # Convert the process output into a list
        string(REGEX REPLACE ";" "\\\\;" _NUMPY_VALUES ${_NUMPY_VALUES})
        string(REGEX REPLACE "\n" ";" _NUMPY_VALUES ${_NUMPY_VALUES})
        list(GET _NUMPY_VALUES 0 NUMPY_VERSION)
        list(GET _NUMPY_VALUES 1 NUMPY_INCLUDE_DIRS)

        # Make sure all directory separators are '/'
        string(REGEX REPLACE "\\\\" "/" NUMPY_INCLUDE_DIRS ${NUMPY_INCLUDE_DIRS})

        # Get the major and minor version numbers
        string(REGEX REPLACE "\\." ";" _NUMPY_VERSION_LIST ${NUMPY_VERSION})
        list(GET _NUMPY_VERSION_LIST 0 NUMPY_VERSION_MAJOR)
        list(GET _NUMPY_VERSION_LIST 1 NUMPY_VERSION_MINOR)
        list(GET _NUMPY_VERSION_LIST 2 NUMPY_VERSION_PATCH)
        string(REGEX MATCH "[0-9]*" NUMPY_VERSION_PATCH ${NUMPY_VERSION_PATCH})
        math(EXPR NUMPY_VERSION_DECIMAL
            "(${NUMPY_VERSION_MAJOR} * 10000) + (${NUMPY_VERSION_MINOR} * 100) + ${NUMPY_VERSION_PATCH}")
        
    endif(NUMPY_COMMAND_RESULT MATCHES 0)
endif(NOT NUMPY_FOUND AND PYTHONINTERP_FOUND)

find_package_handle_standard_args(  NUMPY
                                    REQUIRED_VARS NUMPY_INCLUDE_DIRS
                                    VERSION_VAR NUMPY_VERSION
                                 )
