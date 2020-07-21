PROGRAM MTRNGTEST1
  use iso_c_binding
  USE MRandom
  IMPLICIT NONE

  interface
     subroutine abort() bind(c)
     end subroutine abort
  end interface
  
  INTEGER, PARAMETER    :: seed = 42, realiter = 1000, intiter = 1000, min = 8, max = 16
  REAL(KIND=KIND(0.D0)) :: res
  INTEGER               :: ires, i

  CALL sgrnd(seed)
  do i = 1, realiter
     res = grnd()
     if (res < 0.0d0 .or. res > 1.0d0) then
        write (*, *) "grnd() out of range"
        call abort()
     end if
  end do
  do i = 1, intiter
     ires = randint(min, max)
     if (ires < min .or. ires > max) then
        write (*, *) "randint() out of range"
        call abort()
     end if
  end do
  write (*,*) "success"
  
END PROGRAM MTRNGTEST1
