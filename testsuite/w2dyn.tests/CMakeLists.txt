cmake_minimum_required(VERSION 3.18)
find_package(PythonInterp REQUIRED)

enable_testing()
add_test(NAME 1-syntax_check
         COMMAND ${PYTHON_EXECUTABLE} -m compileall ${PROJECT_SOURCE_DIR})
