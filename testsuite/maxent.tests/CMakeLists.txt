add_executable(1-getgridunit 1-getgridunit.F90)
set (maxenttests 1-getgridunit)

foreach(iter ${maxenttests})
set_target_properties(${iter} PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(${iter} MAXENTLIB)
target_link_libraries(${iter} ${LAPACK_LIBRARIES})
endforeach(iter)

enable_testing()
add_test(1-getgridunit 1-getgridunit)
