! compile with
! gfortran -I /tmp/w2dynamics/lib/ -L /tmp/w2dynamics/lib/  1-getgridunit.F90
! /tmp/w2dynamics/lib/tmp/w2dynamics/src/maxent/MaximumEntropy.o -llapack -lblas
! /tmp/w2dynamics/bin/maxent/CMakeFiles/MAXENTLIB.dir/
PROGRAM MAXENT_GETGRIDUNIT

  USE MMaximumEntropy
  IMPLICIT NONE
  
  INTEGER :: Ngrid, i
  REAL(KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: grid, unit
  
  NGrid = 100
  ALLOCATE(grid(Ngrid), unit(Ngrid))
  integrationmethod = 0
  ! setup dummy init data
  DO i = 1, Ngrid
    grid(i) = i-1
  ENDDO
  CALL getGridUnit(grid, unit, Ngrid)
  ! check results
  DO i = 2, Ngrid
    IF(unit(i) .NE. unit(i-1)) STOP 2 ! We assume identity at the bit level
  ENDDO
  ! at least test the call for the trapezoid integration
  integrationmethod = 1
  DO i = 1, Ngrid
    grid(i) = i-1
  ENDDO
  CALL getGridUnit(grid, unit, Ngrid)
  DEALLOCATE(grid, unit)
  write (*,*) "success"
  
END PROGRAM MAXENT_GETGRIDUNIT
