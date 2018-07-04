! compile with
! gfortran -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ main.f90 ../../Libraries/Modules/modules_90.a -llapack -lblas 

PROGRAM MTRNGTEST1

  USE MRandom
  IMPLICIT NONE
  
  INTEGER :: seed
  REAL(KIND=KIND(0.D0)) :: res
  
  seed = 42

  CALL sgrnd(seed)
  res = grnd()
  write (*,*) "success"
  
END PROGRAM MTRNGTEST1
