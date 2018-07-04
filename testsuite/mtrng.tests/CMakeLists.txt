cmake_minimum_required(VERSION 2.8.9)
Project(mtrng.tests CXX Fortran)

add_executable(1-basic-call 1-basic-call.F90)
set (mtrngtests 1-basic-call)

foreach(iter ${mtrngtests})
target_link_libraries(${iter} mtrng)
endforeach(iter)

enable_testing()
add_test(1-basic-call 1-basic-call)
