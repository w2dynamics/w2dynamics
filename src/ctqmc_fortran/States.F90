
!> This is the 
!! module for generating all possible configurations for a given number of 
!! bands and electrons - it is possible to sort the states by N_t (total
!! number of electrons) and Sz_t (total spin in z-direction)
!! this module is intended to be used with the module Hloc to create a
!! Hamiltonian and creation, annihilation matrices suitable for a
!! continuous time QMC
module MStates
!===============================================================================
use MParameters
use MSparseMatrix
use MAusgabe

   !> \brief Structure holding the system´s occupancy number basis
   !!
   !! If there are good quantum numbers in this basis (the Hamiltonian takes a
   !! block-diagonal form), then the basis states of the disconnected subspaces
   !! are grouped together using \c TSubStates.
   type :: TStates                          
      !> the number of orbitals
      integer                          :: NBands
      !> the total number of states in the fock space
      integer                          :: NStates
      !> the number of blocks in the hamiltonian
      integer                          :: NSStates
      !> the number of states in the biggest block
      integer                          :: NStatesMax
      !> this lookup array has two functions:
      !! (:, 1) tell you which index belongs to which state
      !! (:, 2) it tells you which state (index) belongs to which block value 
      integer,pointer                  :: StatesSStates(:,:)=>null()
      Type(TSubStates),pointer         :: SubStates(:)=>null()
   end type TStates

   type :: TStates_pointer
      type(TStates),pointer :: ptr
   end type TStates_pointer

   type :: TPsis_EB
      real(KINDR),allocatable :: Op(:,:)
   end type TPsis_EB

   type :: TSubStates
      integer                          :: NStates,Offset
      integer,pointer                  :: Connect(:,:,:)=>null()
      integer,pointer                  :: States(:)=>null()
      type(TSparseMatrix),pointer      :: Psis(:,:,:)=>null()
      type(TPsis_EB), pointer          :: Psis_EB(:,:,:)=>null()
      type(TSparseMatrix),pointer      :: Hamiltonian=>null()
      real(KINDR),pointer              :: EVal(:)=>null()
      real(KINDR),pointer              :: EVec(:,:)=>null()
      real(KINDR),pointer              :: S2(:,:)=>null()
   end type TSubStates

contains
!===============================================================================
subroutine init_States(this)
!===============================================================================
   type(TStates)                       :: this
!local
   integer                             :: iSt,iB
   character(vallength)                :: H
   integer,allocatable                 :: occ(:)

   H=adjustl(get_String_Parameter("Hamiltonian"))

   this%NBands=get_Integer_Parameter("Nd")
   this%NStates=2**(2*this%NBands)
   allocate(this%StatesSStates(0:this%NStates-1,3))

   if(H.eq.'ReadIn')then
      open(500,file='occup.dat')
      allocate(occ(this%NBands*2))
      this%StatesSStates(:,1)=0
      this%StatesSStates(:,2)=-1
      do iSt=0,this%NStates-1
         occ=0
         read(500,*)occ
         do iB=1,size(occ)
            if(occ(iB).eq.1)this%StatesSStates(iSt,1)=ibset(this%StatesSStates(iSt,1),iB-1)
         enddo
      enddo
      close(500)
      deallocate(occ)
   else
      this%StatesSStates(:,:)=-1
      do iSt=0,this%NStates-1
         this%StatesSStates(iSt,1)=iSt
      enddo
   endif

end subroutine init_States


!> this subroutine initializes the substates via an array which tells
!! the subroutine which state (index) belongs to which substate (value)
subroutine init_substates(this, nsubstates, states2substates)
   type(TStates)                       :: this
   integer, intent(in)                 :: nsubstates
   integer, intent(in), allocatable    :: states2substates(:)
   !> loop variables
   integer                             :: st1, st2, isst, n
   !> index of the state in a block
   integer                             :: blockidx
   !> saves the states belonging to one substate
   integer                             :: temp(0:this%NStates)

   this%NSStates = nsubstates
   this%NStatesMax = 0
   this%StatesSStates(:,2) = states2substates(:)

   allocate(this%SubStates(0:this%NSStates-1))
   this%SubStates(0)%Offset=0
   do st1 = 0, this%NSStates-1
      blockidx = 0
      do st2 = 0, this%NStates-1
         if(this%StatesSStates(this%statessstates(st2,1),2).eq.st1)then
            temp(blockidx) = this%StatesSStates(st2,1)
            blockidx = blockidx + 1
         endif
      enddo
      this%SubStates(st1)%NStates = blockidx
      if (blockidx > this%NStatesMax) this%NStatesMax = blockidx
      if(st1.gt.0)&
         this%SubStates(st1)%Offset=this%SubStates(st1-1)%Offset+this%SubStates(st1-1)%NStates
      allocate(this%SubStates(st1)%States(0:blockidx-1))
      this%SubStates(st1)%States(0:blockidx-1)=temp(0:blockidx-1)
   enddo

        !StatesSStates(:,3) gives a integer-representation of a ket ( 0 0 ... 0 0 1 0 0 .. 0 0)
   n=0
   do isst=0,this%NSStates-1
!          write(unit=*,fmt=*) "--> isst", isst
      do st1=0,this%substates(isst)%NStates-1
!        write(unit=*,fmt=*) "n", n
!        write(unit=*,fmt=*) "st1", st1
!        write(unit=*,fmt=*) "zustand", this%substates(isst)%States(st1)
         this%StatesSStates(this%substates(isst)%States(st1),3)=n
         n=n+1
      enddo
   enddo

end subroutine 

!> reads the quantum numbers from the parameters and gives you 
!! the number of blocks of the hamiltonian and a lookup array
!! telling you which state (index) belongs to which block 
subroutine qns2substates(this, nsubstates, states2substates)

   type(TStates)                       :: this
   !> loop variables
   integer                             :: st1, st2
   !> the number of blocks in the hamiltonian = number of substates
   integer                             :: nsubstates
   !> saves the states belonging to one substate
   integer, allocatable                :: states2substates(:)
   !> a string containing the quantum numbers
   character(vallength)                :: qns
   integer                             :: hst1, hst2

   allocate(states2substates(0 : this%nstates-1))
   states2substates(:) = -1
   qns = get_String_Parameter("QuantumNumbers")
   nsubstates = 0

   outer: do st1 = 0, this%nstates-1
      hst1 = this%StatesSStates(st1,1)
      if(states2substates(hst1).ne.-1)cycle outer
      nsubstates = nsubstates + 1
      states2substates(hst1) = nsubstates - 1
      inner: do st2 = st1, this%nstates-1
         hst2 = this%StatesSStates(st2,1)
         if(states2substates(hst2).ne.-1)cycle inner
         if( index(QNs,"Nt").ne.0.and.(get_Nt(this,hst1) .ne. get_Nt(this,hst2))) cycle inner
         if(index(QNs,"Szt").ne.0.and..not.eq(get_Szt(this,hst1),get_Szt(this,hst2)))cycle inner
         if(index(QNs,"Qzt").ne.0.and..not.eq(get_Qzt(this,hst1),get_Qzt(this,hst2)))cycle inner
         if(index(QNs,"Azt").ne.0.and..not.eq(get_Azt(this,hst1),get_Azt(this,hst2)))cycle inner
         if(index(QNs,"Lzt").ne.0.and..not.eq(get_Lzt(this,hst1),get_Lzt(this,hst2)))cycle inner
         if(index(QNs,"Jzt").ne.0.and..not.eq(get_Jzt(this,hst1),get_Jzt(this,hst2)))cycle inner
         states2substates(hst2) = nsubstates-1
      enddo inner
   enddo outer
end subroutine qns2substates


!===============================================================================
integer function get_Nt(this, State)
!===============================================================================
   type(TStates)                       :: this
!input
   integer,intent(in)                  :: State
!local
   integer                             :: iB

   get_Nt = 0
   do iB = 0, this%NBands-1
      if(btest(State,iB              )) get_Nt=get_Nt+1
      if(btest(State,iB + this%NBands)) get_Nt=get_Nt+1
   enddo
end function get_Nt

!===============================================================================
real(KINDR) function get_Szt(this,State)
!===============================================================================
   type(TStates)                       :: this
!input
   integer,intent(in)                  :: State
!local
   integer                             :: iB
   get_Szt=0d0
   do iB=0,this%NBands-1
      if(btest(State,iB))get_Szt=get_Szt-0.5d0
      if(btest(State,iB+this%NBands))get_Szt=get_Szt+0.5d0
   enddo
end function get_Szt

!===============================================================================
real(KINDR) function get_Qzt(this,State)
!===============================================================================
   type(TStates)                       :: this
!input
   integer,intent(in)                  :: State
!local
   integer                             :: iB
   get_Qzt=0d0
   do iB=0,this%NBands-1
      if((btest(State,iB).and..not.btest(State,iB+this%NBands))&
         .or.(.not.btest(State,iB).and.btest(State,iB+this%NBands)))get_Qzt=get_Qzt+10**iB
   enddo
end function get_Qzt

!===============================================================================
real(KINDR) function get_Azt(this,State)
!===============================================================================
   type(TStates)                       :: this
!input
   integer,intent(in)                  :: State
!local
   integer                             :: iB
   get_Azt=0d0
   do iB=0,this%NBands-1
      if(btest(State,iB))get_Azt=get_Azt+10d0**(iB)
      if(btest(State,iB+this%NBands))get_Azt=get_Azt+10d0**(iB+this%NBands)
   enddo
end function get_Azt

!===============================================================================
integer function get_eg_up(State)
!===============================================================================
!input
   integer,intent(in)                  :: State
!local
   integer                             :: iB
   get_eg_up=0
   do iB=3,4
      if(btest(State,iB))get_eg_up=get_eg_up+1
   enddo
end function get_eg_up

!===============================================================================
integer function get_eg_do(this,State)
!===============================================================================
   type(TStates)                       :: this
!input
   integer,intent(in)                  :: State
!local
   integer                             :: iB
   get_eg_do=0
   do iB=3,4
      if(btest(State,iB+this%NBands))get_eg_do=get_eg_do+1
   enddo
end function get_eg_do

!===============================================================================
integer function get_t2g_up(State)
!===============================================================================
!input
   integer,intent(in)                  :: State
!local
   integer                             :: iB
   get_t2g_up=0
   do iB=0,2
      if(btest(State,iB))get_t2g_up=get_t2g_up+1
   enddo
end function get_t2g_up

!===============================================================================
integer function get_t2g_do(this,State)
!===============================================================================
   type(TStates)                       :: this
!input
   integer,intent(in)                  :: State
!local
   integer                             :: iB
   get_t2g_do=0
   do iB=0,2
      if(btest(State,iB+this%NBands))get_t2g_do=get_t2g_do+1
   enddo
end function get_t2g_do

!===============================================================================
real(KINDR) function get_Lzt(this,State)
!===============================================================================
   type(TStates)                       :: this
!input
   integer,intent(in)                  :: State
!local
   integer                             :: iB,iS
   get_Lzt=0d0
   do iS=1,2
      do iB=1,this%NBands
         if(btest(State,(iB-1)+this%NBands*(iS-1)))get_Lzt=get_Lzt+dble(iB)-dble(this%NBands+1)/2d0
      enddo
   enddo
end function get_Lzt

!===============================================================================
real(KINDR) function get_Jzt(this,State)
!===============================================================================
   type(TStates)                       :: this
!input
   integer,intent(in)                  :: State
   get_Jzt=get_Lzt(this,State)+get_Szt(this,State)
end function get_Jzt

!===============================================================================
integer pure function get_occ(this,state,b,sp)
!===============================================================================
   type(TStates),intent(in)            :: this
!input
   integer,intent(in)                  :: state,b,sp
   get_occ=ibits(state,(b-1)+this%nbands*(sp-1),1)
end function get_occ

!===============================================================================
character(50) function get_OccStr(this,State)
!===============================================================================
   type(TStates)                       :: this
!input
   integer                             :: State
!local
   integer                             :: iB,iS

   get_OccStr(:)=" "
   do iS=1,2
   get_OccStr(len_trim(get_OccStr)+1:len_trim(get_OccStr)+1)="|"
   do iB=1,this%NBands
      if(get_Occ(this,State,iB,iS).eq.1)then
         get_OccStr(len_trim(get_OccStr)+1:len_trim(get_OccStr)+1)="1"
      elseif(get_Occ(this,State,iB,iS).eq.0)then
         get_OccStr(len_trim(get_OccStr)+1:len_trim(get_OccStr)+1)="0"
      else
         write(stdout,'("Error in get_OccStr: Got a number .ne. 0/1")')
      endif
   enddo
   get_OccStr(len_trim(get_OccStr)+1:len_trim(get_OccStr)+1)=">"
   enddo
   get_OccStr(len_trim(get_OccStr)+1:len_trim(get_OccStr)+1)=";"
end function get_OccStr

!===============================================================================
subroutine print_States(this)
!===============================================================================
   type(TStates)                       :: this
!local
   integer                             :: iSt
   character(len=30)                   :: chr
   
   do iSt=0,this%NStates-1
      write(stdout,'("State: ",2I5)')this%StatesSStates(iSt,1),iSt
      write(chr,'("|",2f5.1,"> ")')real(get_Nt(this,this%StatesSStates(iSt,1))),get_Szt(this,this%StatesSStates(iSt,1))
      write(stdout,*)chr,get_OccStr(this,this%StatesSStates(iSt,1))
   enddo
end subroutine print_States

!===============================================================================
subroutine print_JanStates(this)
!===============================================================================
   type(TStates)                       :: this
!local
   integer                             :: iSt,iSSt
   character(len=70)                   :: chr

   write(stdout,'("# State            Nt   Szt   t2g_up  t2g_do  eg_up  eg_do SSt  St")')
   do iSSt=0,this%NSStates-1
      do iSt=0,this%SubStates(iSSt)%NStates-1
         write(chr,'(I5,f6.2,7I4)')&
         get_Nt(this,this%SubStates(iSSt)%States(iSt)),get_Szt(this,this%SubStates(iSSt)%States(iSt)),&
         get_t2g_up(this%SubStates(iSSt)%States(iSt)),get_t2g_do(this,this%SubStates(iSSt)%States(iSt)),&
         get_eg_up(this%SubStates(iSSt)%States(iSt)),get_eg_do(this,this%SubStates(iSSt)%States(iSt)),&
         iSSt,iSt,this%SubStates(iSSt)%Offset+iSt
         write(stdout,*)trim(get_OccStr(this,this%SubStates(iSSt)%States(iSt))),chr
      enddo
   enddo
end subroutine print_JanStates


!===============================================================================
subroutine print_SubStates(this)
!===============================================================================
   type(TStates)                       :: this
!local
   integer                             :: iSSt,iE
   character(5)                        :: formatstr


   do iSSt=0,this%NSStates-1
      write(stdout,'("Nt: ",I5,", Szt: ",f5.2,"; NStates: ",I3,"; ISSt: ",I3," ;")',advance="no")&
         get_Nt(this,this%SubStates(iSSt)%States(0)),&
         get_Szt(this,this%SubStates(iSSt)%States(0)),&
         this%SubStates(iSSt)%NStates,&
         iSSt
      do iE=0,size(this%SubStates(iSSt)%States)-1
         formatstr="(A  )"
         write(formatstr(3:4),'(I2)')len_trim(get_OccStr(this,this%SubStates(iSSt)%States(iE)))
         write(stdout,formatstr,advance="no")get_OccStr(this,this%SubStates(iSSt)%States(iE))
      enddo
      write(stdout,*)
   enddo
end subroutine print_SubStates

!===============================================================================
subroutine print_densitymatrix_basis(this, iter_no)
!===============================================================================
   type(TStates)                       :: this
   integer                             :: iter_no
!local
   integer                             :: iSSt,iE
   character(5)                        :: formatstr
   CHARACTER(LEN=:), ALLOCATABLE :: fname

   ALLOCATE(CHARACTER(LEN=32)::fname)
   write(*,*) "iter_no", iter_no
   if(iter_no.lt.10)then
      write(fname,'(A25,2H00I1,A4)')'densitymatrix_basis_iter_', iter_no, '.dat'
   elseif(iter_no.lt.100)then
      write(fname,'(A25,1H0I2,A4)')'densitymatrix_basis_iter_', iter_no, '.dat'
   elseif(iter_no.lt.1000)then
      write(fname,'(A25,I3,A4)')'densitymatrix_basis_iter_', iter_no, '.dat'
   endif

   write(*,*) "fname", fname

   OPEN(345,file=fname,form='formatted',status='replace')

   do iSSt=0,this%NSStates-1
   do iE=0,size(this%SubStates(iSSt)%States)-1

      write(345,'("Nt: ",f5.2," , Szt: ",f5.2," , Qzt: ",f12.2," ; NStates: ",I4," , ISSt: ",I4," , integer repr.: ",I4," , state: ")',advance="no")&
         real(get_Nt(this,this%SubStates(iSSt)%States(iE))),&
         get_Szt(this,this%SubStates(iSSt)%States(iE)),&
         get_Qzt(this,this%SubStates(iSSt)%States(iE)),&
         this%SubStates(iSSt)%NStates,&
         iSSt,&
         this%SubStates(iSSt)%States(iE)
         formatstr="(A  )"
         write(formatstr(3:4),'(I2)')len_trim(get_OccStr(this,this%SubStates(iSSt)%States(iE)))
         write(345,formatstr,advance="yes")get_OccStr(this,this%SubStates(iSSt)%States(iE))

   enddo
   enddo

   close(345)

end subroutine print_densitymatrix_basis

!===============================================================================
subroutine dest_States(this)
!===============================================================================
   type(TStates)                       :: this
   integer                             :: iSS,iB,iS,iCA

   if(associated(this%SubStates))then
      do iSS=0,this%NSStates-1
         if(associated(this%SubStates(iSS)%Connect))&
            deallocate(this%SubStates(iSS)%Connect)
         if(associated(this%SubStates(iSS)%States))&
            deallocate(this%SubStates(iSS)%States)
         if(associated(this%SubStates(iSS)%EVal))&
            deallocate(this%SubStates(iSS)%EVal)
         if(associated(this%SubStates(iSS)%EVec))&
            deallocate(this%SubStates(iSS)%EVec)
         if(associated(this%SubStates(iSS)%S2))&
            deallocate(this%SubStates(iSS)%S2)
         if(associated(this%SubStates(iSS)%Hamiltonian))then
            call dest_SparseMatrix(this%SubStates(iSS)%Hamiltonian)
            deallocate(this%SubStates(iSS)%Hamiltonian)
         endif
         do iB=1,this%NBands
            do iS=1,2
               do iCA=1,2
                  if(associated(this%SubStates(iSS)%Psis))&
                     call dest_SparseMatrix(this%SubStates(iSS)%Psis(iB,iS,iCA))
                  if (associated(this%SubStates(iSS)%Psis_EB)) then
                     if(allocated(this%SubStates(iSS)%Psis_EB(iB,iS,iCA)%Op))&
                        deallocate(this%SubStates(iSS)%Psis_EB(iB,iS,iCA)%Op)
                  end if
               enddo
            enddo
         enddo
         if(associated(this%SubStates(iSS)%Psis))&
            deallocate(this%SubStates(iSS)%Psis)
         if(associated(this%SubStates(iSS)%Psis_EB))&
            deallocate(this%SubStates(iSS)%Psis_EB)
      enddo
   endif
   if(associated(this%StatesSStates))&
      deallocate(this%StatesSStates)
   if(associated(this%SubStates))&
      deallocate(this%SubStates)
end subroutine dest_States

!===============================================================================
end module MStates
!===============================================================================


#ifdef States_Test

!===============================================================================
program States_Prog
!===============================================================================
use MParameters
use MStates
!local
type(TStates)                          :: DStates
integer                             :: nsubstates, is
integer, allocatable                :: states2substates(:)

call read_ParameterFile("Parameters.in")

write(stdout,*)"================================================================"
write(stdout,*)"            initialising the module states"
write(stdout,*)"================================================================"
call init_States(DStates)
write(stdout,*)"================================================================"
write(stdout,*)"                  number of bands"
write(stdout,*)"================================================================"
write(stdout,'("                      "I3)'),DStates%NBands
write(stdout,*)"================================================================"
write(stdout,*)"              printing initialised states"
write(stdout,*)"================================================================"
call qns2substates(DStates, nsubstates, states2substates)
call init_SubStates(DStates, nsubstates, states2substates)
call print_States(DStates)
!stop
write(stdout,*)"================================================================"
write(stdout,*)"              printing superstates"
write(stdout,*)"================================================================"
call print_SubStates(DStates)
write(stdout,*)"================================================================"
write(stdout,*)"              printing StatesSStates"
write(stdout,*)"================================================================"
call ausgabe(DStates%StatesSStates(:,1),shape(DStates%StatesSStates(:,1)),"DStates%StatesSStates(:,1)",3,.false.)
call ausgabe(DStates%StatesSStates(:,2),shape(DStates%StatesSStates(:,2)),"DStates%StatesSStates(:,2)",3,.false.)
call ausgabe(DStates%StatesSStates(:,3),shape(DStates%StatesSStates(:,3)),"DStates%StatesSStates(:,3)",3,.false.)
write(stdout,*)"================================================================"
write(stdout,*)"           destructing the module states"
write(stdout,*)"================================================================"
call dest_States(DStates)
call dest_Parameters()
!===============================================================================
end program States_Prog
!===============================================================================


#endif
