! compile with
! gfortran -I /tmp/w2dynamics/lib/ -L /tmp/w2dynamics/lib/  1-getgridunit.F90
! /tmp/w2dynamics/lib/tmp/w2dynamics/src/maxent/MaximumEntropy.o -llapack -lblas
! /tmp/w2dynamics/bin/maxent/CMakeFiles/MAXENTLIB.dir/
PROGRAM CTQMCTIMEEVOLVEINEIGENBASIS

  USE MStates
  USE MTrace
  IMPLICIT NONE
  
    interface
        subroutine time_evolve_in_eigenbasis_old(substate, state_vector, state_vector_evolved, slen, tau, egs)
            USE MStates
            type(TSubStates), intent(in)        :: substate
            integer, intent(in)                 :: slen
            real(KINDR), intent(inout)          :: state_vector_evolved(0:slen-1)
            real(KINDR), intent(in)             :: state_vector(0:slen-1)
            real(KINDR), intent(in)             :: egs
            real(KINDR), intent(in), value      :: tau
        end subroutine
    end interface
  
  INTEGER :: k
  Type(TSubStates) :: substate
  REAL(KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE         :: bra_in, bra_out, bra_out_old
  REAL(KIND=KIND(0.D0)), dimension(:), ALLOCATABLE, target :: evals
  REAL(KIND=KIND(0.D0)) :: tau, egs
  integer :: slen, i
  
  ! set up input data
  egs = 0
  tau = 0.2
  slen = 10
  allocate(evals(0:slen-1), bra_in(slen), bra_out(slen), bra_out_old(slen))
  do i = 0, slen-1
      evals(i) = 0.1*i
  enddo
  bra_in = 1
  substate%NStates = slen
  substate%Offset = 0
  substate%Eval => evals
  ! see what the current implementation does
  call time_evolve_in_eigenbasis(substate, bra_in, bra_out, slen, tau+0.1, egs)
  ! check results against reference implementation
  call time_evolve_in_eigenbasis_old(substate, bra_in, bra_out_old, slen, tau+0.1, egs)
  
  ! perform FP comparison
  do i = 1, slen
    if (abs(bra_out_old(i) - bra_out(i)) > max(abs(bra_out_old(i)), abs(bra_out(i)))*1E-14) then
        write (*,*) "error at ", i, "got ", bra_out(i), "expected ", bra_out_old(i)
        STOP 2
    endif
  enddo
  
  write (*,*) "success"
  
END PROGRAM CTQMCTIMEEVOLVEINEIGENBASIS

!===============================================================================
subroutine time_evolve_in_eigenbasis_old(substate, state_vector, state_vector_evolved, slen, tau, egs)
!===============================================================================
!input
   USE MStates
   type(TSubStates), intent(in)        :: substate
   integer, intent(in)                 :: slen
   real(KINDR), intent(inout)          :: state_vector_evolved(0:slen-1)
   real(KINDR), intent(in)             :: state_vector(0:slen-1)
   real(KINDR), intent(in)             :: egs
   real(KINDR), intent(in), value      :: tau
   integer :: iS, iS_in_SS

   iS_in_SS=0
   do iS=substate%Offset, substate%Offset+substate%NStates-1
       state_vector_evolved(is)=exp(-(substate%Eval(iS_in_SS)-egs)*tau)*state_vector(is)
       iS_in_SS=iS_in_SS+1
   enddo

end subroutine time_evolve_in_eigenbasis_old
