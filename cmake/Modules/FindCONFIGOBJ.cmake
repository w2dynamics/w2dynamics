# - FindCONFIGOBJ
# By Florian Goth, fgoth@physik.uni-wuerzburg.de
# Find ConfigObj
# This module defines:
# CONFIGOBJ_FOUND
# CONFIGOBJ_VERSION
# CONFIGOBJ_VERSION_MAJOR
# CONFIGOBJ_VERSION_MINOR
# CONFIGOBJ_VERSION_PATCH

# Finding ConfigObj involves calling the Python interpreter
#Check wether we have already searched the Python interpreter
if(NOT PYTHONINTERP_FOUND)
    find_package(PythonInterp REQUIRED)
endif()

  execute_process(COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "exec(\"try:\\n import configobj;\\n print(configobj.__version__);\\nexcept:\\n exit(1)\")"
      OUTPUT_VARIABLE CONFIGOBJ_VERSION
      RESULT_VARIABLE CONFIGOBJ_COMMAND_RESULT
      OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT CONFIGOBJ_FOUND AND PYTHONINTERP_FOUND)
#    if(NOT CONFIGOBJ_COMMAND_RESULT MATCHES 0)
#        message("ConfigObj import failure:\n${_H5PY_ERROR_VALUE}")
#    endif()
    if(CONFIGOBJ_COMMAND_RESULT MATCHES 0)
        set(CONFIGOBJ_FOUND TRUE)
        # Get the major and minor version numbers
        string(REGEX REPLACE "\\." ";" CONFIGOBJ_VERSION_LIST ${CONFIGOBJ_VERSION})
        list(GET CONFIGOBJ_VERSION_LIST 0 CONFIGOBJ_VERSION_MAJOR)
        list(GET CONFIGOBJ_VERSION_LIST 1 CONFIGOBJ_VERSION_MINOR)
        list(GET CONFIGOBJ_VERSION_LIST 2 CONFIGOBJ_VERSION_PATCH)
        string(REGEX MATCH "[0-9]*" CONFIGOBJ_VERSION_PATCH ${CONFIGOBJ_VERSION_PATCH})
        math(EXPR CONFIGOBJ_VERSION_DECIMAL
            "(${CONFIGOBJ_VERSION_MAJOR} * 10000) + (${CONFIGOBJ_VERSION_MINOR} * 100) + ${CONFIGOBJ_VERSION_PATCH}")
    endif(CONFIGOBJ_COMMAND_RESULT MATCHES 0)
endif(NOT CONFIGOBJ_FOUND AND PYTHONINTERP_FOUND)

find_package_handle_standard_args(  CONFIGOBJ
                                    REQUIRED_VARS CONFIGOBJ_VERSION
                                    VERSION_VAR CONFIGOBJ_VERSION
                                 )
