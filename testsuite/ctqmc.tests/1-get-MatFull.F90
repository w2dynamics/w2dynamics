! compile with
! gfortran -I /tmp/w2dynamics/lib/ -L /tmp/w2dynamics/lib/  1-getgridunit.F90
! /tmp/w2dynamics/lib/tmp/w2dynamics/src/maxent/MaximumEntropy.o -llapack -lblas
! /tmp/w2dynamics/bin/maxent/CMakeFiles/MAXENTLIB.dir/
PROGRAM CTQMCGETMATDETFULL

  USE MMatrixUpdate
  IMPLICIT NONE
  
  INTEGER :: k, i, j
  REAL(KIND=KIND(0.D0)), DIMENSION(:, :), ALLOCATABLE :: mat
  type(TLogDet) :: mydet
  
  k = 10
  ALLOCATE(mat(k,k))
  ! setup dummy init data
  DO i = 1, k
  DO j = 1, k
    mat(i, j) = i+j
  ENDDO
  ENDDO
  CALL get_MatLogDetFull(k, mat, mydet)
  ! check results
  write (*,*) "success"
  
END PROGRAM CTQMCGETMATDETFULL
