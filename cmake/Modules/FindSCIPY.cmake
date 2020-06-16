# - FindSciPy
# By Florian Goth, fgoth@physik.uni-wuerzburg.de
# Find SciPy includes and library
# This module defines:
# SCIPY_FOUND
# SCIPY_VERSION             - the version of SciPy found as a string
# SCIPY_VERSION_MAJOR       - the major version number of SciPy
# SCIPY_VERSION_MINOR       - the minor version number of SciPy
# SCIPY_VERSION_DECIMAL     - e.g. version 1.6.1 is 10601
# SCIPY_INCLUDE_DIRS        - path to the SciPy include files

# Finding SciPy involves calling the Python interpreter
#Check wether we have already searched the Python interpreter
if(NOT PYTHONINTERP_FOUND)
    find_package(PythonInterp REQUIRED)
endif()

if(NOT SCIPY_FOUND AND PYTHONINTERP_FOUND)
    execute_process(COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "exec(\"try:\\n import numpy;\\n import scipy;\\n print(scipy.__version__);\\n print(numpy.get_include())\\nexcept:\\n exit(1)\")"
      OUTPUT_VARIABLE _SCIPY_VALUES
      RESULT_VARIABLE SCIPY_COMMAND_RESULT
      OUTPUT_STRIP_TRAILING_WHITESPACE)

#    if(NOT SCIPY_COMMAND_RESULT MATCHES 0)
#        message("SciPy import failure:\n${_SCIPY_ERROR_VALUE}")
    if(SCIPY_COMMAND_RESULT MATCHES 0)
        # Convert the process output into a list
        string(REGEX REPLACE ";" "\\\\;" _SCIPY_VALUES ${_SCIPY_VALUES})
        string(REGEX REPLACE "\n" ";" _SCIPY_VALUES ${_SCIPY_VALUES})
        list(GET _SCIPY_VALUES 0 SCIPY_VERSION)
        list(GET _SCIPY_VALUES 1 SCIPY_INCLUDE_DIRS)

        # Make sure all directory separators are '/'
        string(REGEX REPLACE "\\\\" "/" SCIPY_INCLUDE_DIRS ${SCIPY_INCLUDE_DIRS})

        # Get the major and minor version numbers
        string(REGEX REPLACE "\\." ";" _SCIPY_VERSION_LIST ${SCIPY_VERSION})
        list(GET _SCIPY_VERSION_LIST 0 SCIPY_VERSION_MAJOR)
        list(GET _SCIPY_VERSION_LIST 1 SCIPY_VERSION_MINOR)
        list(GET _SCIPY_VERSION_LIST 2 SCIPY_VERSION_PATCH)
        string(REGEX MATCH "[0-9]*" NUMPY_VERSION_PATCH ${NUMPY_VERSION_PATCH})
        math(EXPR SCIPY_VERSION_DECIMAL
            "(${SCIPY_VERSION_MAJOR} * 10000) + (${SCIPY_VERSION_MINOR} * 100) + ${SCIPY_VERSION_PATCH}")
        
    endif(SCIPY_COMMAND_RESULT MATCHES 0)
endif(NOT SCIPY_FOUND AND PYTHONINTERP_FOUND)

find_package_handle_standard_args(  SCIPY
                                    REQUIRED_VARS SCIPY_INCLUDE_DIRS
                                    VERSION_VAR SCIPY_VERSION
                                 )
