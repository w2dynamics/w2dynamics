if(NOT PYTHONINTERP_FOUND)
        find_package(PythonInterp REQUIRED)
endif()

FUNCTION(INSTALL_VIA_PIP module RESULT_NAME)
#check for pip. The pip installer needs the python-xml module.
execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip freeze OUTPUT_VARIABLE PYTHON_PACKAGE_LIST)
    if ("${PYTHON_PACKAGE_LIST}" STREQUAL "")
        execute_process(COMMAND pip freeze OUTPUT_VARIABLE PYTHON_PACKAGE_LIST)
        if ("${PYTHON_PACKAGE_LIST}" STREQUAL "")
            message(WARNING "Failed to find pip. Pip is required to automatically install ${module}. Will now install pip.")
            file(DOWNLOAD https://bootstrap.pypa.io/get-pip.py ./get-pip.py STATUS SUCCESS)
            LIST(GET ${SUCCESS} 0 PIPDLFAILED) # first element of download status != 0 indicates error
            if (${PIPDLFAILED})
                message(WARNING "CMake was not able to download pip. Trying with a direct call to curl")
                execute_process(COMMAND curl -k https://bootstrap.pypa.io/get-pip.py -o get-pip.py RESULT_VARIABLE result)
            endif()
            execute_process(COMMAND ${PYTHON_EXECUTABLE} ./get-pip.py --user RESULT_VARIABLE RES)
#            message(${RES})
            if(${RES})#FIXME: Check wether this works for other versions...
                message(WARNING "Not able to successfully execute Pip. assuming this is due to an old version of Python.")
                file(DOWNLOAD https://files.pythonhosted.org/packages/45/db/4fb9a456b4ec4d3b701456ef562b9d72d76b6358e0c1463d17db18c5b772/pip-1.5.6.tar.gz ./pip-1.5.6.tar.gz
		  STATUS SUCCESS
		  EXPECTED_HASH SHA256=b1a4ae66baf21b7eb05a5e4f37c50c2706fa28ea1f8780ce8efe14dcd9f1726c
		  )
                LIST(GET ${SUCCESS} 0 PIPDLFAILED)
                if (${PIPDLFAILED})
                    message("CMake was not able to download pip-1.1. Trying with a direct call to curl")
                    execute_process(COMMAND curl -k https://files.pythonhosted.org/packages/45/db/4fb9a456b4ec4d3b701456ef562b9d72d76b6358e0c1463d17db18c5b772/pip-1.5.6.tar.gz -o ./pip-1.5.6.tar.gz RESULT_VARIABLE result)
                    execute_process(COMMAND tar -xzvf ./pip-1.5.6.tar.gz)
                    execute_process(COMMAND python ./pip-1.5.6/setup.py install)  # pip needs setuptools.
                endif()
            endif()
        endif()
    endif()
# pip is present

message(STATUS "Installing ${module}")

execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip "install" "--user" ${module} RESULT_VARIABLE SUCCESS)
if (NOT "${SUCCESS}" STREQUAL "0")
  execute_process(COMMAND pip "install" "--user" ${module} RESULT_VARIABLE SUCCESS)
  if (NOT "${SUCCESS}" STREQUAL "0")
    message(WARNING "Failed to automatically install ${module}. Please install manually")
    SET(${RESULT_NAME}  ${SUCCESS} PARENT_SCOPE)
  endif()
endif()
ENDFUNCTION()
