cmake_minimum_required(VERSION 2.8.9)
Project(ctqmc.tests CXX Fortran)
find_package(PythonLibs REQUIRED)

add_executable(1-get-MatFull 1-get-MatFull.F90)
add_executable(2-time-evolve-in-eigenbasis 2-time-evolve-in-eigenbasis.F95)
set (ctqmctests 1-get-MatFull 2-time-evolve-in-eigenbasis)

foreach(iter ${ctqmctests})
target_link_libraries(${iter} -Wl,${CMAKE_SOURCE_DIR}/auxiliaries/CTQMC${F2PY_SUFFIX}) #exploit that cmake doesn't touch stuff beginning with "-"
target_link_libraries(${iter} ${PYTHON_LIBRARIES})
target_link_libraries(${iter} ${LAPACK_LIBRARIES})
target_link_libraries(${iter} ${NFFT_LIBRARIES})
target_link_libraries(${iter} -lstdc++)
endforeach(iter)

enable_testing()
add_test(1-get-MatFull 1-get-MatFull)
add_test(2-time-evolve-in-eigenbasis 2-time-evolve-in-eigenbasis)
