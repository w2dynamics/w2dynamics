# - FindH5PY
# By Florian Goth, fgoth@physik.uni-wuerzburg.de
# Find H5Py includes and library
# This module defines:
# H5PY_FOUND               - TRUE if H5PY is found
# H5PY_VERSION             - the version of H5Py found as a string
# H5PY_VERSION_MAJOR       - the major version number of H5Py
# H5PY_VERSION_MINOR       - the minor version number of H5Py
# H5PY_VERSION_DECIMAL     - e.g. version 1.6.1 is 10601
# H5PY_INCLUDE_DIRS        - path to the H5Py include files

# Finding H5Py involves calling the Python interpreter
#Check wether we have already searched the Python interpreter
if(NOT PYTHONINTERP_FOUND)
    find_package(PythonInterp REQUIRED)
endif()

if(NOT H5PY_FOUND AND PYTHONINTERP_FOUND)
    execute_process(COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "exec(\"try:\\n import h5py;\\n print(h5py.version.version);\\nexcept:\\n exit(1)\")"
      OUTPUT_VARIABLE _H5PY_VALUES
      RESULT_VARIABLE H5PY_COMMAND_RESULT
      OUTPUT_STRIP_TRAILING_WHITESPACE)

#    if(NOT H5PY_COMMAND_RESULT MATCHES 0)
#        message("SciPy import failure:\n${_H5PY_ERROR_VALUE}")
    if(H5PY_COMMAND_RESULT MATCHES 0)
        # Convert the process output into a list
        string(REGEX REPLACE ";" "\\\\;" _H5PY_VALUES ${_H5PY_VALUES})
        string(REGEX REPLACE "\n" ";" _H5PY_VALUES ${_H5PY_VALUES})
        list(GET _H5PY_VALUES 0 H5PY_VERSION)
#        list(GET _H5PY_VALUES 1 H5PY_INCLUDE_DIRS)

        # Make sure all directory separators are '/'
#        string(REGEX REPLACE "\\\\" "/" H5PY_INCLUDE_DIRS ${H5PY_INCLUDE_DIRS})

        # Get the major and minor version numbers
        string(REGEX REPLACE "\\." ";" _H5PY_VERSION_LIST ${H5PY_VERSION})
        list(GET _H5PY_VERSION_LIST 0 H5PY_VERSION_MAJOR)
        list(GET _H5PY_VERSION_LIST 1 H5PY_VERSION_MINOR)
        list(GET _H5PY_VERSION_LIST 2 H5PY_VERSION_PATCH)
        string(REGEX MATCH "[0-9]*" H5PY_VERSION_PATCH ${H5PY_VERSION_PATCH})
        math(EXPR H5PY_VERSION_DECIMAL
            "(${H5PY_VERSION_MAJOR} * 10000) + (${H5PY_VERSION_MINOR} * 100) + ${H5PY_VERSION_PATCH}")
        
    endif(H5PY_COMMAND_RESULT MATCHES 0)
endif(NOT H5PY_FOUND AND PYTHONINTERP_FOUND)

find_package_handle_standard_args(  H5PY
                                    REQUIRED_VARS H5PY_VERSION
                                    VERSION_VAR H5PY_VERSION
                                 )
