cmake_minimum_required(VERSION 2.8.9)
Project(mtrng.tests Fortran CXX)

add_executable(1-basic-call 1-basic-call.F90)
set (mtrngtests 1-basic-call)

foreach(iter ${mtrngtests})
set_target_properties(${iter} PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(${iter} mtrng)
endforeach(iter)

enable_testing()
add_test(1-basic-call 1-basic-call)
