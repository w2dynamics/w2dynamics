cmake_minimum_required(VERSION 3.18)
Project(w2dyn.tests)
find_package(PythonInterp REQUIRED)

enable_testing()
if (PYTHON_VERSION_MAJOR GREATER 2)
  add_test(NAME 1-syntax_check
    COMMAND ${PYTHON_EXECUTABLE} -m compileall ${PROJECT_SOURCE_DIR})
else ()
  add_test(NAME 1-syntax_check
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compileall.py ${PROJECT_SOURCE_DIR})
endif()
