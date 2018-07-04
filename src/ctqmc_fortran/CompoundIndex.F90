!> This calculates general band spin index patterns of arbitrary length 'N' from
!! a single flavor index 'ind'. A one-to-one mapping requires 'Nbands'.
!! The inverse transform is also included.
!!
!! Authors: Josef Kaufmann, Patrik Gunacker
!! 
!! Compiler Call: ifort CompoundIndex.F90 -DCompoundIndex_Test
!===============================================================================
module MCompoundIndex
!===============================================================================

contains
!===============================================================================
!> converting an index into a band-spin pattern
!> 'N' number of operators
!> 'ind' general flavor index
!> 'bs' band-spin array of length N
!> 'b' band array of length N
!> 's' spin array of length N
subroutine index2component_general(Nbands, N, ind, bs, b, s)
!===============================================================================
  integer,intent(in) :: Nbands,N,ind
  integer,intent(out) :: bs(N), b(N), s(N)
  integer :: tmp(0:N),ind_tmp
  integer :: i

  ! the proposed back conversion assumes the indices are
  ! given form 0 to max-1
  ind_tmp = ind - 1
  do i=0,N
    tmp(i) = (2*Nbands)**(N-i)
  end do

  do i=1,N
    bs(i) = ind_tmp/tmp(i)
    s(i) = mod(bs(i),2)
    b(i) = (bs(i)-s(i))/2
    ind_tmp = ind_tmp - tmp(i)*bs(i)
  end do
  s=s+1
  b=b+1
end subroutine index2component_general

!===============================================================================
!> converting a band-spin pattern into an index
!> 'N' number of operators
!> 'ind' general flavor index
!> 'bs' band-spin array of length N
!> 'b' band array of length N
!> 's' spin array of length N
subroutine component2index_general(Nbands, N, ind, bs, b, s)
!===============================================================================
   integer,intent(in) :: Nbands,N
   integer,intent(in) :: b(N), s(N)
   integer,intent(out) :: ind, bs(N)
   integer :: i

   ind = 1
   do i=1,N
      bs(i) = 2*(b(i)-1) + s(i)
      ind = ind + (2*Nbands)**(N-i)*(bs(i)-1)
   enddo

end subroutine component2index_general

end module MCompoundIndex

#ifdef CompoundIndex_Test

program CompoundIndex
  use MCompoundIndex
  implicit none
  integer :: ind
  integer, allocatable :: bs(:), b(:),s(:)
  integer,parameter :: Nbands=3,N=2
  integer :: i,j
  
  allocate(bs(N),b(N),s(N))
  
  write(*,*) "Testing transform and backtransform"
  write(*,*) "No further output if everything working"

  do i=1,(2*Nbands)**N
    call index2component_general(NBands,N,i,bs,b,s)
    call component2index_general(Nbands,N,ind,bs,b,s)
    if(i .ne. ind) write(*,*), ind, b, s
  end do

  deallocate(bs,b,s)
end program

#endif

