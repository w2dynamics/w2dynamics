cmake_minimum_required(VERSION 2.8.9)
Project(maxent.tests CXX Fortran)
find_package(PythonLibs REQUIRED)

add_executable(1-getgridunit 1-getgridunit.F90)
set (maxenttests 1-getgridunit)

foreach(iter ${maxenttests})
target_link_libraries(${iter} -Wl,${CMAKE_SOURCE_DIR}/w2dyn/maxent/MAXENT${F2PY_SUFFIX})
target_link_libraries(${iter} ${PYTHON_LIBRARIES})
endforeach(iter)

enable_testing()
add_test(1-getgridunit 1-getgridunit)
