!===============================================================================
!> This module generates the Hamiltonian, creation and annihilation Operators.
!! They can then be stored in the type specified in MStates.
module MOperator
!===============================================================================
use MParameters
use MStates
use MSparseMatrix
use MAngularMomentum
use MAusgabe
   
   type :: TPsis
!! creation and annihilation operators,last index 1 creation,2 annihilation
      type(TOperator),pointer                :: Psis(:,:,:)=>null()
   end type   

   type :: TOperator
      type(TSubOperator),pointer             :: SubOps(:)=>null()
   end type TOperator
   
   type :: TSubOperator
      real(KINDR),pointer                    :: Op(:,:)=>null()
   end type TSubOperator

contains


!===============================================================================
subroutine u_readin(unit,u_matrix)
!===============================================================================
     !< Sets values of the U matrix based on a file
     integer, intent(in) :: unit

     integer :: line
     integer :: nbands
     character(len=200) :: msg, linestr
     character(len=2) :: c1, c2, c3, c4
     integer :: b1, b2, b3, b4
     real(KINDR) :: u_value
     logical :: header
     real(KINDR),allocatable :: u_matrix(:,:,:,:)


     line = 0
     linestr = '[no line read yet]'
     header = .true.
     do
         line = line + 1
         read(unit, '(A)', err=98, end=99) linestr
         if(linestr(1:1) == '#') cycle
         if(header) then
             ! Read number of bands and compare for safety
             read(linestr, *) nbands
             if(2*nbands /= size(u_matrix,1)) stop '[u_readin] Wrong #bands'
             header = .false.
         else
             read(linestr, '(4(BN, I1, A1, 1X), F20.10)', err=98) &
                 b1, c1, b2, c2, b3, c3, b4, c4, u_value
             if (b1 < 1 .or. b1 > nbands) call raise('Invalid band 1')
             if (b2 < 1 .or. b2 > nbands) call raise('Invalid band 2')
             if (b3 < 1 .or. b3 > nbands) call raise('Invalid band 3')
             if (b4 < 1 .or. b4 > nbands) call raise('Invalid band 4')
             u_matrix(2*(b1-1)+char_spin(c1), 2*(b2-1)+char_spin(c2), &
                      2*(b3-1)+char_spin(c3), 2*(b4-1)+char_spin(c4)) = u_value
         endif
     enddo
  98 call raise(msg)
  99 return
 contains
     integer function char_spin(c) result(sp)
         character(len=1), intent(in) :: c
         select case(c)
         case('u','U')
             sp = 1
         case('d','D')
             sp = 2
         case default
             call raise('Invalid spin, expecting u or d, got `' // c // '''')
         end select
     end function

     subroutine raise(msg)
         character(len=*), intent(in) :: msg
         character(len=200) :: fname

         inquire(unit, name=fname)
         write (0,"('Error reading U matrix file ',A,' on line',I7,':')") trim(fname), line
         write (0,"(A,/,'>>> ',A)") trim(msg), trim(linestr)
         stop
     end subroutine
end subroutine u_readin


!===============================================================================
subroutine init_TOperator(this,DStates)
!===============================================================================
   type(TOperator)                           :: this
!input
   type(TStates)                             :: DStates
!local
   integer                                   :: iSSt,NStates
   allocate(this%SubOps(0:DStates%NSStates-1))
   do iSSt=0,DStates%NSStates-1
      NStates=DStates%SubStates(iSSt)%NStates
      allocate(this%SubOps(iSSt)%Op(0:NStates-1,0:NStates-1))
      this%SubOps(iSSt)%Op(:,:)=0d0
   enddo
end subroutine init_TOperator

!===============================================================================
!> This subroutine allocates the arrays necessary to store the creation and 
!! annihilation operators.
!! @params this a pointer to type TPsis which is filled by this routine
!! @params DStates holds information about the states of the system
subroutine init_Psi(this,DStates)
!===============================================================================
   type(TPsis)                               :: this
!input
   type(TStates)                             :: DStates
!local
   integer                                   :: iB,iS,iI,iSSt1,iSSt2,iBand,iSpin
   integer                                   :: iSt1,iSt2,State1,State2
   real(KINDR)                               :: sig
   real(KINDR),allocatable                   :: temp(:,:)
   real(KINDR),pointer                       :: temp1(:,:)

   ! Fortran 90 does not guarantee nullified pointers, esp. over multiple calls
   nullify(temp1)

! allocate the necessary space for the dense psi matrices
   allocate(this%Psis(DStates%NBands,2,2)) 
   do iB=1,DStates%NBands
   do iS=1,2
   do iI=1,2
      allocate(this%Psis(iB,iS,iI)%SubOps(0:DStates%NSStates-1))
      do iSSt1=0,DStates%NSStates-1
         this%Psis(iB,iS,iI)%SubOps(iSSt1)%Op=>Null()
      enddo
   enddo
   enddo
   enddo
   allocate(temp(0:DStates%NStates,0:DStates%NStates))
!! generate the marices by counting hoppings over electrons 
!! (read phase factor *-1 from every hop)
   do iSSt1=0,DStates%NSStates-1
   if (associated(DStates%SubStates(iSSt1)%Connect)) deallocate(DStates%SubStates(iSSt1)%Connect)
   allocate(DStates%SubStates(iSSt1)%Connect(DStates%NBands,2,2))
   DStates%SubStates(iSSt1)%Connect(:,:,:)=-1
   do iSSt2=0,DStates%NSStates-1
   if(iSSt1.eq.iSSt2)cycle
   do iSpin=1,2
   do iBand=1,DStates%NBands
      do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
      do iSt2=0,DStates%SubStates(iSSt2)%NStates-1
         State1=DStates%SubStates(iSSt1)%States(iSt1)
         State2=DStates%SubStates(iSSt2)%States(iSt2)
         iI=(iBand-1)+DStates%NBands*(iSpin-1)
         if(ieor(State1,State2).eq.ibset(0,iI))then
            sig=1d0
!            do iB=1,iBand-1
!               iI=(iB-1)+DStates%NBands*(iSpin-1)
 
            spin: do iS=1,2
            do iB=1,DStates%NBands
               iI=(iB-1)+DStates%NBands*(iS-1)
               if((iBand-1)+DStates%NBands*(iSpin-1).eq.iI)exit spin
               if(btest(State1,iI))then 
                  sig=sig*(-1d0)
               endif
            enddo
            enddo spin
            if(State2.gt.State1)then
               DStates%SubStates(iSSt1)%Connect(iBand,iSpin,1)=iSSt2
            elseif(State2.lt.State1)then
               DStates%SubStates(iSSt1)%Connect(iBand,iSpin,2)=iSSt2
            endif
            temp(iSt2,iSt1)=sig
         else 
            temp(iSt2,iSt1)=0d0
         endif
      enddo
      enddo
      if(DStates%SubStates(iSSt1)%Connect(iBand,iSpin,1).eq.iSSt2)then
         allocate(this%Psis(iBand,iSpin,1)%SubOps(iSSt1)%&
            Op(0:DStates%SubStates(iSSt2)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1))
         this%Psis(iBand,iSpin,1)%SubOps(iSSt1)%Op(:,:)=&
            temp(0:DStates%SubStates(iSSt2)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1)
      elseif(DStates%SubStates(iSSt1)%Connect(iBand,iSpin,2).eq.iSSt2)then
         allocate(this%Psis(iBand,iSpin,2)%SubOps(iSSt1)%&
            Op(0:DStates%SubStates(iSSt2)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1))
         this%Psis(iBand,iSpin,2)%SubOps(iSSt1)%Op(:,:)=&
            temp(0:DStates%SubStates(iSSt2)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1)
      endif
   enddo
   enddo
   enddo
   enddo


!! check anti commutation relation
!! @todo probably put this stuff into a seperate routine
!   do iSSt1=0,DStates%NSStates-1
!   do iBand=1,DStates%NBands
!   do iSpin=1,2
!   do iB=1,DStates%NBands
!   do iS=1,2
!      temp1=>Null()
!      SSt1=DStates%SubStates(iSSt1)%Connect(iB,iS,2) 
!      if(SSt1.ne.-1)then
!      if(DStates%SubStates(SSt1)%Connect(iBand,iSpin,1).eq.iSSt1)then 
!         allocate(temp1(0:DStates%SubStates(iSSt1)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1))
!! create space for first matrix in temp1
!         temp1=matmul(this%Psis(iBand,iSpin,1)%Subops(SSt1)%Op,this%Psis(iB,iS,2)%SubOps(iSSt1)%Op)
!      endif
!      endif
!      SSt1=DStates%SubStates(iSSt1)%Connect(iBand,iSpin,1)
!      if(SSt1.ne.-1)then
!      if(DStates%SubStates(SSt1)%Connect(iB,iS,2).eq.iSSt1)then
!         if(associated(temp1))then
!         temp1=temp1+matmul(this%Psis(iB,iS,2)%Subops(SSt1)%Op,this%Psis(iBand,iSpin,1)%SubOps(iSSt1)%Op)
!         else
!         allocate(temp1(0:DStates%SubStates(iSSt1)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1))
!         temp1=matmul(this%Psis(iB,iS,2)%Subops(SSt1)%Op,this%Psis(iBand,iSpin,1)%SubOps(iSSt1)%Op)
!         endif
!      endif
!      endif
!      if(.not.associated(temp1))cycle
!      do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
!         if(iB.eq.iBand.and.iS.eq.iSpin)then
!         if(.not.eq(temp1(iSt1,iSt1),1d0))write(stdout,'("Error in diag operator::init_Psi")')
!         else
!         if(.not.eq(temp1(iSt1,iSt1),0d0))write(stdout,'("Error in diag operator::init_Psi")')
!         endif
!         do iSt2=0,DStates%SubStates(iSSt1)%NStates-1
!            if(.not.eq(temp1(iSt1,iSt2),0d0).and.iSt1.ne.iSt2)write(stdout,*)"Error in offdiag operator::init_Psi"
!         enddo
!      enddo 
!      if(associated(temp1))deallocate(temp1)
!   enddo
!   enddo
!   enddo
!   enddo
!   enddo
   deallocate(temp)
end subroutine init_Psi

!> we analyze the hamiltonian for the minimal block size and build
!! up the array states2substates for later building the right
!! substates.

!! This proceeds as follows:
!!- In each old superstate, 'nets' [subsets] of contained states that
!!  are connected by application of the Hamiltonian (and thus mix under
!!  time evolution) are constructed. A subset x is represented by a
!!  slice nets(x, :), whose entries are the indices of the states it
!!  contains in the States array of the old superstate or -1 at
!!  'unoccupied' places (the nets array is preallocated to be able to
!!  hold the maximum number of subsets with the maximum number of
!!  subset members each). For initialization, the state with index 0 is
!!  put into the first subset. Then, one loops over the subsets (index
!!  net), over the contained states (index node) and over all states
!!  (States-index edge) to find new states (edge) with non-zero
!!  Hamiltonian matrix element (dh%subops(sst)%op(nets(net,node),edge)
!!  between them and an already contained state (States-index nets(net,
!!  node)), which are then appended to the subset (net) to be iterated
!!  over in the "node"-loop. At the same time, the assignments of
!!  states (bin. rep.) to newly subdivided superstates are kept track
!!  of in states2substates, and when one has fully completed the
!!  "node"-loop over a subset, another subset needs to be created if a
!!  state in the old superstate that is not contained in a subset yet,
!!  which is initialized with that state.
!!
!!- To fulfill the second condition for choosing superstates, that
!!  application of a creator or annihilator may lead to at most one
!!  other superstate, the mappings from states to superstates after
!!  application of an operator are created (for each flavour & ca),
!!  and if any two states from the same newly subdivided superstate
!!  map to different superstates after the operator application, the
!!  two target superstates are merged together (by assigning all
!!  states contained in the one with higher index to the one with the
!!  lower index and removing the former).
!!
!! we are doing the following for each substate in pseudo code:
!!
!! nets.append(net(nodes.pop))
!!
!! for net in nets
!!   for node in net
!!     for edge in node
!!       if edge.node not in net
!!           net.append(edge.node)
!!           nodes.remove(edge.node)
!!       end if
!!     end for edge
!!   end for node
!!   if nodes are empty
!!     done
!!   else
!!     nets.append(net(nodes.pop))
!!   end if
!! end for net
!!
subroutine analyze_hamiltonian(dh, dstates, nsubstates, states2substates, threshold, log)
   type(TOperator)                        :: dh
   type(TStates),intent(in)               :: dstates
   integer                                :: nsubstates
   integer, allocatable                   :: states2substates(:), cstates2substates(:)
   real(KINDR), intent(in)                :: threshold
   logical, intent(in)                    :: log
   real(KINDR)                            :: maxignored
   integer, allocatable                   :: nets(:,:)
   integer                                :: sst, net, node, edge, netnode
   integer                                :: cstate,nstates
   integer                                :: idx, idx1, idx2
   integer                                :: empty_node, ib, is , ca
   integer                                :: bigger_substate, smaller_substate
   logical                                :: member, change
   
   if (.not.allocated(states2substates))&
      allocate(states2substates(0:dstates%nstates-1))
   if (.not.allocated(cstates2substates))&
      allocate(cstates2substates(0:dstates%nstates-1))
   states2substates(:) = -1
   nsubstates = -1
   nstates = dstates%nstates
   maxignored = 0.0d0

   do sst = 0, dstates%nsstates-1
     allocate(nets(0:dstates%substates(sst)%nstates-1,0:dstates%substates(sst)%nstates-1))
     nets(:,:) = -1
     nets(0,0) = 0
     nsubstates = nsubstates + 1
     states2substates(dstates%substates(sst)%states(0)) = nsubstates
!     write(*,*)dh%subops(sst)%op
     do net = 0, dstates%substates(sst)%nstates-1
       if (nets(net,0).eq.-1) exit
       do node = 0, dstates%substates(sst)%nstates-1
         if (nets(net,node).eq.-1) exit
         do edge = 0, dstates%substates(sst)%nstates-1 
           if (dh%subops(sst)%op(nets(net,node),edge).ne.0d0) then
             if (abs(dh%subops(sst)%op(nets(net, node), edge)) > threshold) then
!             write(*,*)nets(net,node),edge,dh%subops(sst)%op(node,edge)
               member = .false.
               do netnode = 0, dstates%substates(sst)%nstates-1
                 if (nets(net,netnode).eq.edge) then
                   member = .true.
                   exit
                 else if (nets(net,netnode).eq.-1) then
                   empty_node = netnode
                   exit
                 endif
               enddo
               if (.not.member) then
                 nets(net,empty_node) = edge
                 states2substates(dstates%substates(sst)%states(edge)) = nsubstates
               endif
             else
                maxignored = max(maxignored, abs(dh%subops(sst)%op(nets(net, node), edge)))
             endif
           endif
         enddo
       enddo
       do node = 0, dstates%substates(sst)%nstates-1
         if (states2substates(dstates%substates(sst)%states(node)).eq.-1) then
           nets(net+1,0) = node
           nsubstates = nsubstates + 1
           states2substates(dstates%substates(sst)%states(node)) = nsubstates
           exit
         endif
       enddo
     enddo
     deallocate(nets)
   enddo
   nsubstates = nsubstates + 1
!   write(*,*)states2substates
! joing blocks which are doubly connected
   change = .true.
outer:   do while(change)
      change = .false.
      do ib = 1, dstates%nbands
      do is = 1, 2
      do ca = 1, 2
        cstates2substates(:) = -1
        do idx = 0, dstates%nstates-1
          if (.not.btest(dstates%statessstates(idx,1),(ib-1)+dstates%nbands*(is-1)).and.&
              ca == 1)then
             cstate = ibset(dstates%statessstates(idx,1),(ib-1)+dstates%nbands*(is-1))
          else if (btest(dstates%statessstates(idx,1),(ib-1)+dstates%nbands*(is-1)).and.&
                   ca == 2)then
             cstate = ibclr(dstates%statessstates(idx,1),(ib-1)+dstates%nbands*(is-1))
          else
            cstate = -1
          endif
          if (cstate.ne.-1)then
            cstates2substates(dstates%statessstates(idx,1)) = states2substates(cstate)
          endif
        enddo
!        write(*,*)"operator: ", ib,is,ca,cstates2substates
        do idx1 = 0, nstates-1
        do idx2 = idx1+1, nstates-1
          if(states2substates(dstates%statessstates(idx1,1)) == &
             states2substates(dstates%statessstates(idx2,1)) .and.&
             cstates2substates(dstates%statessstates(idx1,1)) /= &
             cstates2substates(dstates%statessstates(idx2,1)) .and.&
             cstates2substates(dstates%statessstates(idx1,1)) /= -1 .and.&
             cstates2substates(dstates%statessstates(idx2,1)) /= -1) then
             change = .true.
             bigger_substate = max(cstates2substates(dstates%statessstates(idx1,1)),&
                                   cstates2substates(dstates%statessstates(idx2,1)))
             smaller_substate = min(cstates2substates(dstates%statessstates(idx1,1)),&
                                    cstates2substates(dstates%statessstates(idx2,1)))
!             write(*,*)"substate: ",bigger_substate," becomes ",smaller_substate," since "
!             write(*,*)"substate: ",states2substates(dstates%statessstates(idx1,1))," are the same"
!             write(*,*)"but ",cstates2substates(dstates%statessstates(idx1,1)) &
!              &,cstates2substates(dstates%statessstates(idx2,1)), " are not"
             do idx = 0, nstates-1
               if (states2substates(dstates%statessstates(idx,1)) == bigger_substate) then
                  states2substates(dstates%statessstates(idx,1)) = smaller_substate
               endif
               if (states2substates(dstates%statessstates(idx,1)) .gt. bigger_substate) then
                  states2substates(dstates%statessstates(idx,1)) = &
                    states2substates(dstates%statessstates(idx,1))-1
               endif
             enddo
             cycle outer
          end if
        enddo
        enddo
!        write(*,*)"states2substates:", states2substates
      enddo
      enddo
      enddo
   enddo outer
   nsubstates =  maxval(states2substates)+1
!   write(*,*) states2substates
   if (log) then
      write (*,"('Analyzed Hamiltonian, found',I5,' blocks')") nsubstates
      if (maxignored > 0.0d0) then
         write (*, "('Maximum size of entries smaller than EPSBLOCK potentially ignored by automatic partitioning: ', E12.3)") maxignored
      endif
   endif
end subroutine analyze_hamiltonian

!===============================================================================
!> This subroutine initialises the Hamiltonian container and the one particle
!! part of the Hamiltonian.
!! @params this A pointer to TOperator which after this routine stores the
!! block diagonal Hamiltonian
!! @params DStates A pointer to TStates
!! @params DPsis A pointer to TPsis
!! @params muimp an array containing the impurity chemical potential
!! @params phonon_shift an array containing the shifts due to the phonons
subroutine init_Hamiltonian(this,DStates,muimp)
!===============================================================================
   type(TOperator)                        :: this
!input
   type(TStates),intent(in)               :: DStates
   real(KINDR),intent(in)                 :: muimp(DStates%NBands,2)
!local
   integer                                :: iSSt1,iSt1,iB1,iS1,State1
   
!! allocate the necessary space
   call init_TOperator(this,DStates)

! set the one-particle quantities of the hamiltonian
   do iSSt1=0,DStates%NSStates-1
   do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
     State1=DStates%SubStates(iSSt1)%States(iSt1)
     do iB1=1,DStates%NBands
       do iS1=1,2
         if(get_Occ(DStates,State1,iB1,iS1).ne.1)cycle
         this%SubOps(iSSt1)%Op(iSt1,iSt1) = this%SubOps(iSSt1)%Op(iSt1, iSt1) +& 
                                            muimp(iB1, iS1) 
       enddo
     enddo
   enddo
   enddo

end subroutine init_Hamiltonian

!===============================================================================
subroutine oneparticle_hamiltonian(h, DStates, DPsis, muimp)
!===============================================================================
      !< \brief Add interaction term corresponding to U to the Hamiltonian

      type(toperator), intent(inout) :: h
      type(tstates), intent(in)      :: dstates
      type(tpsis), intent(in)        :: dpsis
      real(KINDR)                    :: muimp(:,:,:,:)

      integer :: sst, sst1, sst2
      integer :: b1, b2, sp1, sp2
      real(KINDR) :: muimp_value

      if(.not.associated(dpsis%psis)) stop '[oneparticle_hamiltonian] DAFUQ?'

      !! allocate the necessary space
      call init_TOperator(h,DStates)

      do b1 = 1,dstates%NBands
      do sp1 = 1,2
        do b2 = 1,dstates%NBands
        do sp2 = 1,2
            ! Check that U_1234 is non-zero (there is an interaction)
            muimp_value = muimp(b1,sp1,b2,sp2)
            if(muimp_value /= 0) then
               do sst = 0,dstates%NSStates-1
                  ! Make the connections: check if annihilator corresponding
                  ! to c[b_4,s_4] connects to any block (sst4 >= 0), then
                  ! continue  with that superstate and check for the next
                  ! operators until a chain c^+ c^+ c c is built.
                  !sst4 = DStates%SubStates(sst)%Connect(b3,sp3,2)
                  !if (sst4 < 0) cycle
                  !sst3 = DStates%SubStates(sst4)%Connect(b4,sp4,2)
                  !if (sst3 < 0) cycle

                  sst2 = DStates%SubStates(sst)%Connect(b2,sp2,2)
                  if (sst2 < 0) cycle
                  sst1 = DStates%SubStates(sst2)%Connect(b1,sp1,1)
                  if (sst1 < 0) cycle
                  if (sst /= sst1) then
                    cycle
!                    write (0,"('U(',3(2I1,','),2I1,')=',F10.2)") &
!                    b1, sp1, b2, sp2, b3, sp2, b4, sp1, u_value
!                    write (0,"('sst chain:',5(I5,'<'))") sst1, sst2, sst3, sst4, sst
!                    stop 'Hamiltonian not block-diagnonal'

                  endif

                  ! Now, compute the product c^+ c acting on the
                  ! corresponding subspace
                  h%SubOps(sst)%Op = h%SubOps(sst)%Op + muimp_value * &
                           matmul(DPsis%Psis(b1,sp1,1)%SubOps(sst2)%Op, &
                                  DPsis%Psis(b2,sp2,2)%SubOps(sst)%Op)
               enddo
            endif
               ! End U check
        enddo
        enddo
      enddo
      enddo
end subroutine oneparticle_hamiltonian


!===============================================================================
!> This subroutine generates the S^2 operator using 1/2(S+S- + S-S+)+S_z^2.
!! It uses S=\sum_{n,\sigma,sigma^\prime}c^\dagger_{n,\sigma}\frac{\hbar \sigma_{\sigma,\sigma^\prime}}{2}c_{n,\sigma^\prime}
subroutine init_S2alt2(this,DStates,DPsis)
!===============================================================================
   type(TOperator)                        :: this
!input
   type(TStates),intent(in)               :: DStates
   type(TPsis),intent(in)                 :: DPsis
!local
   integer                                :: iSSt1,iSt1,iSt2,iB1,iS1,iS2,iB3,iS3,iS4
   integer                                :: SSt1,SSt2,SSt3,SSt4
   integer                                :: xyz
   real(KINDR)                            :: sm1(2,2),sm2(2,2)
   real(KINDR),pointer                    :: temp1(:,:),temp2(:,:),temp3(:,:),temp4(:,:)

   nullify(temp1)
   nullify(temp2)
   nullify(temp3)
   nullify(temp4)

   call init_TOperator(this,DStates)
! generate the operator using creation/annihilation operators
   do xyz=1,3
      if (xyz.eq.1) then
         sm1=Spm(1,:,:)
         sm2=Spm(2,:,:)/2.
      elseif (xyz.eq.2) then
         sm1=Spm(2,:,:)
         sm2=Spm(1,:,:)/2.
      elseif (xyz.eq.3) then
         sm1=real(PauliSM(3,:,:))
         sm2=real(PauliSM(3,:,:)/4.)
      endif
      do iSSt1=0,DStates%NSStates-1
         do iB1=1,DStates%NBands
         do iS1=1,2
         do iS2=1,2
            if(associated(temp1))deallocate(temp1)
            if(associated(temp2))deallocate(temp2)
            if(abs(sm1(iS2,iS1)).eq.0d0)then
               cycle
            endif
!    start from SuperState iSSt1 --> reach SSt1
            SSt1=DStates%SubStates(iSSt1)%Connect(iB1,iS1,2)
            if(SSt1.eq.-1)cycle    !Skip the not connected SuperStates
            allocate(temp1(0:DStates%SubStates(SSt1)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1))
            temp1=DPsis%Psis(iB1,iS1,2)%SubOps(iSSt1)%Op
!    start from SuperState SSt1 --> reach SSt2
            SSt2=DStates%SubStates(SSt1)%Connect(iB1,iS2,1)
            if(SSt2.eq.-1)cycle    !Skip the not connected SuperStates
            allocate(temp2(0:DStates%SubStates(SSt2)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1))
            temp2=matmul(DPsis%Psis(iB1,iS2,1)%SubOps(SSt1)%Op,temp1)
            do iB3=1,DStates%NBands
            do iS3=1,2
            do iS4=1,2
               if(associated(temp3))deallocate(temp3)
               if(associated(temp4))deallocate(temp4)
               if(abs(sm2(iS4,iS3)).eq.0d0)then
                  cycle
               endif
   !    start from SuperState SSt2 --> reach SSt3
               SSt3=DStates%SubStates(SSt2)%Connect(iB3,iS3,2)
               if(SSt3.eq.-1)cycle    !Skip the not connected SuperStates
!               write(*,*)SSt3
               allocate(temp3(0:DStates%SubStates(SSt3)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1))
               temp3=matmul(DPsis%Psis(iB3,iS3,2)%SubOps(SSt2)%Op,temp2)
   !    start from SuperState SSt3 --> reach SSt4
               SSt4=DStates%SubStates(SSt3)%Connect(iB3,iS4,1)
               if(SSt4.eq.-1)cycle    !Skip the not connected SuperStates
               allocate(temp4(0:DStates%SubStates(SSt4)%NStates-1,0:DStates%SubStates(iSSt1)%NStates-1))
               temp4=matmul(DPsis%Psis(iB3,iS4,1)%SubOps(SSt3)%Op,temp3)
               do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
               do iSt2=0,DStates%SubStates(SSt4)%NStates-1
                  if(temp4(iSt2,iSt1).ne.0d0)then
                     this%SubOps(iSSt1)%Op(iSt2,iSt1)=this%SubOps(iSSt1)%Op(iSt2,iSt1)+&
                        sm1(iS2,iS1)*sm2(iS4,iS3)*temp4(iSt2,iSt1)
                  endif
               enddo
               enddo
            enddo
            enddo
            enddo
         enddo
         enddo
         enddo
      enddo
   enddo
end subroutine init_S2alt2

!===============================================================================
! higher precision can be reached by using dsyevr?
subroutine diag_Operator(this,DStates,EVec,EVal)
!===============================================================================
   type(TOperator)                        :: this
!input
   type(TStates)                          :: DStates
   real(KINDR), allocatable               :: EVal(:)
   type(TOperator)                        :: EVec
!output
!local
   real(KINDR), allocatable               :: WORK(:)
   integer                                :: iSSt,INFO,NElements

   do iSSt=0,DStates%NSStates-1
      NElements=DStates%SubStates(iSSt)%NStates
      EVec%SubOps(iSSt)%Op=this%SubOps(iSSt)%Op
#ifdef LAPACK77_Interface
      allocate(WORK(3*NElements))
      call DSYEV('V','U',NElements,EVec%SubOps(iSSt)%Op(0, 0),NElements,&
         EVal(DStates%SubStates(iSSt)%Offset), WORK,3*NElements,INFO)
      deallocate(WORK)
#endif
#ifdef LAPACK95_Interface
      call SYEV(EVec%SubOps(iSSt)%Op,&
         EVal(DStates%SubStates(iSSt)%Offset:DStates%SubStates(iSSt)%NStates-1),&
         'V','U',INFO)
#endif
      if(INFO.ne.0)then
         write(stderr,*) "Error in MOperator::Diag_Operator", INFO
         stop
      endif
   enddo


end subroutine Diag_Operator

!===============================================================================
subroutine set_Psis(this,DStates)
!===============================================================================
   type(TPsis)                               :: this
!input
   type(TStates)                             :: DStates
!local
   integer                                   :: iSSt,SSt,iB,iS,CA
   do iSSt=0,DStates%NSStates-1
   allocate(DStates%SubStates(iSSt)%Psis(DStates%NBands,2,2))
   do iB=1,DStates%NBands
   do iS=1,2
   do CA=1,2
      SSt=DStates%SubStates(iSSt)%Connect(iB,iS,CA)
      if(SSt.ne.-1)then
         call init_SparseMatrix(DStates%SubStates(iSSt)%Psis(iB,iS,CA),&
                              this%Psis(iB,iS,Ca)%SubOps(iSST)%Op,&
                              DStates%SubStates(SSt)%NStates,&
                              DStates%SubStates(iSSt)%NStates)
      endif
   enddo
   enddo
   enddo
   enddo
end subroutine set_Psis


!===============================================================================
subroutine set_EB_Psis(DTransformed_Psis,DStates)
!===============================================================================
type(TPsis)                               :: DTransformed_Psis
!input
type(TStates), intent(inout)              :: DStates
!local
integer                                   :: iSSt,SSt,iB,iS,CA
integer                                   :: form(2)

do iSSt=0,DStates%NSStates-1
allocate(DStates%SubStates(iSSt)%Psis_EB(DStates%NBands,2,2))
do iB=1,DStates%NBands
do iS=1,2
do CA=1,2
   SSt=DStates%SubStates(iSSt)%Connect(iB,iS,CA)
   if(SSt.ne.-1) then
      form=shape(DTransformed_Psis%Psis(iB,iS,CA)%SubOps(isst)%Op)
      allocate(DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op(form(1),form(2)))
!     call ausgabe(form,shape(form),"form",3,.false.)
!     call ausgabe(shape(DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op),shape(shape(DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op)),&
!     &"shape(DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op)",3,.false.)

      DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op=DTransformed_Psis%Psis(iB,iS,CA)%SubOps(isst)%Op
!     call ausgabe(DTransformed_Psis%Psis(iB,iS,CA)%SubOps(isst)%Op,shape(DTransformed_Psis%Psis(iB,iS,CA)%SubOps(isst)%Op),"DTransformed_Psis%Psis(iB,iS,CA)%SubOps(isst)%Op",3,.false.)
!     call ausgabe(DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op,shape(DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op),"DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op",3,.false.)
!     write(unit=*,fmt=*) "-------------"
   endif
enddo
enddo
enddo
enddo

   !call print_ev(DStates)
end subroutine set_EB_Psis


!===============================================================================
subroutine set_Hamiltonian(this,DStates)
!===============================================================================
   type(TOperator)                           :: this
!input
   type(TStates)                             :: DStates
!local
   integer                                   :: iSSt
   do iSSt=0,DStates%NSStates-1
      allocate(DStates%SubStates(iSSt)%Hamiltonian)
      call init_SparseMatrix(DStates%SubStates(iSSt)%Hamiltonian,&
                           this%Subops(iSSt)%Op,&
                           DStates%SubStates(iSSt)%NStates,&
                           DStates%SubStates(iSSt)%NStates)
   enddo
end subroutine set_Hamiltonian

!===============================================================================
subroutine set_HEVal(this,DStates)
!===============================================================================
   real(KINDR),allocatable                   :: this(:)
!input
   type(TStates)                             :: DStates
!local
   integer                                   :: iSSt
   do iSSt=0,DStates%NSStates-1
      allocate(DStates%SubStates(iSSt)%EVal(0:DStates%SubStates(iSST)%NStates-1))
      DStates%SubStates(iSSt)%EVal(:)&
         =this(DStates%SubStates(iSSt)%Offset:&
               DStates%SubStates(iSSt)%Offset+DStates%SubStates(iSSt)%NStates-1)
   enddo
end subroutine set_HEval

! setting the eigenvectors, check for numerical zeros and renormalize them
!===============================================================================
subroutine set_HEvec(this,DStates)
!===============================================================================
   type(TOperator)                           :: this
!input
   type(TStates)                             :: DStates
!local
   integer                                   :: iSSt,i,j,NStates
   do iSSt=0,DStates%NSStates-1
      NStates=DStates%SubStates(iSSt)%NStates
      allocate(DStates%SubStates(iSSt)%EVec(0:NStates-1,0:NStates-1))
      DStates%SubStates(iSSt)%EVec=this%SubOps(iSSt)%Op
      do i=0, nstates-1
        do j=0, nstates-1
          if(abs(dstates%substates(isst)%evec(j,i)).lt.get_Real_Parameter("EPSEVEC"))then
            dstates%substates(isst)%evec(j,i)=0d0
          endif
        enddo
        dstates%substates(isst)%evec(:,i)=dstates%substates(isst)%evec(:,i)/&
                                       sum(dstates%substates(isst)%evec(:,i)**2)
      enddo
!! check orthogonality of states
!      do i=0, nstates-1
!        do j=0, nstates-1
!          write(*,*) dot_product(this%SubOps(iSSt)%Op(:,i),&
!                                 this%SubOps(iSSt)%Op(:,j)),&
!                     dot_product(dstates%substates(isst)%evec(:,i),&
!                                 dstates%substates(isst)%evec(:,j))
!        enddo
!      enddo
   enddo
end subroutine set_HEvec

!===============================================================================
subroutine transform_S2(DStates,DS2)
!===============================================================================
!input
   type(TStates)              :: DStates
   type(TOperator)            :: DS2
   integer                    :: isst

   do isst=0, DStates%NSStates-1

      !apply unitary basis transformation to S2:
      !          S2'=U^\dagger S2 U

      DS2%SubOps(isst)%Op=&
      matmul(transpose(DStates%SubStates(isst)%Evec),&
      matmul(DS2%SubOps(isst)%Op,DStates%Substates(isst)%EVec))

   enddo

end subroutine transform_S2


! setting the S2 operator
!===============================================================================
subroutine set_S2(this,DStates)
!===============================================================================
   type(TOperator)                           :: this
!input
   type(TStates)                             :: DStates
!local
   integer                                   :: iSSt,i,j,NStates
   do iSSt=0,DStates%NSStates-1
      NStates=DStates%SubStates(iSSt)%NStates
      if(.not.associated(DStates%SubStates(iSSt)%S2)) then
          allocate(DStates%SubStates(iSSt)%S2(0:NStates-1,0:NStates-1))
      endif
      DStates%SubStates(iSSt)%S2=this%SubOps(iSSt)%Op
   enddo
end subroutine set_S2

!===============================================================================
subroutine print_Psi(this,DStates)
!===============================================================================
   type(TPsis)                               :: this
!input
   type(TStates),intent(in)                  :: DStates
!local
   integer                                   :: iSSt1,iSSt2,iSt1,iSt2,iB,iS,iDag
   character(5)                              :: formatstr

   do iDag=1,2
   write(stdout,*)"================================================================"
      if(iDag.eq.1)write(stdout,'("Creation Operators: ")')
      if(iDag.eq.2)write(stdout,'("Annihilation Operators: ")')
   write(stdout,*)"================================================================"
   do iSSt1=0,DStates%NSStates-1 
   do iSSt2=0,DStates%NSStates-1
   do iB=1,DStates%NBands
   do iS=1,2
      if(iSSt1.eq.iSSt2)cycle
      if(DStates%SubStates(iSSt1)%Connect(iB,iS,iDag).ne.iSSt2)cycle
      write(stdout,'("Super State ",I3," to ",I3," with Orbital,Spin",2I3)') iSSt1,iSSt2,iB,iS
      do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
         formatstr="(A  )"
         write(formatstr(3:4),'(I2)')len_trim(get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt1)))+4
         write(stdout,formatstr,advance='no')get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt1))
      enddo
      write(stdout,*)
      do iSt2=0,DStates%SubStates(iSSt2)%NStates-1
         formatstr="(A  )"
         write(formatstr(3:4),'(I2)')len_trim(get_OccStr(DStates,DStates%SubStates(iSSt2)%States(iSt2)))
         write(stdout,formatstr,advance="no")get_OccStr(DStates,DStates%SubStates(iSSt2)%States(iSt2))
      do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
         write(stdout,'("    ",f4.1)',advance="no")this%Psis(iB,iS,iDag)%SubOps(iSSt1)%Op(iSt2,iSt1)

      enddo
      write(stdout,*)
      enddo
   enddo
   enddo
   enddo
   write(unit=*,fmt=*) " "
   enddo
   enddo

end subroutine print_Psi

!===============================================================================
subroutine print_Operator(this,DStates)
!===============================================================================
   type(TOperator)                           :: this
!input
   type(TStates)                             :: DStates
!local
   integer                                   :: iSSt1,iSt1,iSt2
   character(5)                              :: formatstr

   write(stdout,*)"Operator: "
   do iSSt1=0,DStates%NSStates-1
!      if(get_Nt(DStates,DStates%SubStates(iSSt1)%States(0)).gt.get_Real_Parameter("NElMin")&
!         .or.get_Nt(DStates,DStates%SubStates(iSSt1)%States(0)).lt.get_Real_Parameter("NElMax"))&
!         cycle
      write(*,*)"SubState: ",iSSt1
   do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
      formatstr="(A  )"
      write(formatstr(3:4),'(I2)')len_trim(get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt1)))+4
      write(stdout,formatstr,advance='no')get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt1))
   enddo
   write(stdout,*)
   do iSt2=0,DStates%SubStates(iSSt1)%NStates-1
         formatstr="(A  )"
         write(formatstr(3:4),'(I2)')len_trim(get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt2)))
         write(stdout,formatstr,advance="no")get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt2))
   do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
      if (abs(this%SubOps(iSSt1)%Op(iSt2,iSt1))>1d-4.or.abs(this%SubOps(iSSt1)%Op(iSt2,iSt1))==0)then
        write(stdout,'("   ",f6.2)',advance="no")this%SubOps(iSSt1)%Op(iSt2,iSt1)
      else
        write(stdout,'(" ",e8.1)',advance="no")this%SubOps(iSSt1)%Op(iSt2,iSt1)
      endif
   enddo
   write(stdout,*)
   enddo

   enddo
end subroutine print_Operator


!===============================================================================
subroutine print_ev(DStates)
!===============================================================================
!input
   type(TStates)                             :: DStates
!local
   integer                                   :: iSSt1,iSt1
   character(5)                              :: formatstr
!   real(kindr) :: mean

   do iSSt1=0,DStates%NSStates-1
   write(*,*) "-------isst--------", isst1
!   mean = sum(DStates%SubStates(isst1)%Eval)/max(1,size(DStates%SubStates(isst1)%Eval))
   !write (*, *) "mean ", mean
   !write (*, *) "variance", sqrt(sum((mean-DStates%SubStates(isst1)%Eval)**2)/max(1,size(DStates%SubStates(isst1)%Eval)))
!      if(get_Nt(DStates,DStates%SubStates(iSSt1)%States(0)).gt.get_Real_Parameter("NElMin")&
!         .or.get_Nt(DStates,DStates%SubStates(iSSt1)%States(0)).lt.get_Real_Parameter("NElMax"))&
!         cycle
   do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
      !if (get_nt(dstates,DStates%SubStates(iSSt1)%States(iSt1)) .ne. 2) cycle
      formatstr="(A  )"
      write(formatstr(3:4),'(I2)')len_trim(get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt1)))+4
      write(stdout,formatstr,advance='no')get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt1))
      write(*,*) DStates%SubStates(isst1)%Eval(ist1)
   enddo
      write(*,*) ""
   enddo

end subroutine print_ev

! !===============================================================================
! subroutine print_Operator(this,DStates)
! !===============================================================================
!    type(TOperator)                           :: this
! !input
!    type(TStates)                             :: DStates
! !local
!    integer                                   :: iSSt1,iSt1,iSt2
!    character(5)                              :: formatstr
! 
!    !if(.not.present(fileid)) fileid = stdout
!    do iSSt1=0,DStates%NSStates-1
! !      if(get_Nt(DStates,DStates%SubStates(iSSt1)%States(0)).gt.get_Real_Parameter("NElMin")&
! !         .or.get_Nt(DStates,DStates%SubStates(iSSt1)%States(0)).lt.get_Real_Parameter("NElMax"))&
! !         cycle
!       write(*,"(A10,I2:)")"SUBSTATE: ", iSSt1
!       write(*,"(A14)",advance='no') ""
!    do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
!       formatstr="(A  )"
! !       write(formatstr(3:4),'(I2)')len_trim(get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt1)))+3
! !       write(*,formatstr,advance='no')get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt1))
!         write(*,"(A9)",advance='no')get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt1))
!    enddo
!    write(*,*)
!    do iSt2=0,DStates%SubStates(iSSt1)%NStates-1
!          formatstr="(A  )"
!          write(formatstr(3:4),'(I2)')len_trim(get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt2)))
!          write(*,formatstr,advance="no")get_OccStr(DStates,DStates%SubStates(iSSt1)%States(iSt2))
!    do iSt1=0,DStates%SubStates(iSSt1)%NStates-1
!       if (abs(this%SubOps(iSSt1)%Op(iSt2,iSt1))>1d-4.or.abs(this%SubOps(iSSt1)%Op(iSt2,iSt1))==0)then
!         write(*,'("   ",f6.2)',advance="no")this%SubOps(iSSt1)%Op(iSt2,iSt1)
!       else
!         write(*,'(" ",e8.1)',advance="no")this%SubOps(iSSt1)%Op(iSt2,iSt1)
!       endif
!    enddo
!    write(*,*)
!    enddo
! 
!    enddo
! end subroutine print_Operator

!===============================================================================
subroutine dest_Psis(this)
!===============================================================================
   type(TPsis)                               :: this
!local
   integer                                   :: iB,iS,CA,iSSt

   do iB=1,size(this%Psis(:,1,1))
   do iS=1,2
   do CA=1,2
      if(associated(this%Psis(iB,iS,CA)%SubOps))then
         do iSSt=0,size(this%Psis(iB,iS,CA)%SubOps)-1
            if(associated(this%Psis(iB,iS,CA)%SubOps(iSSt)%Op))&
               deallocate(this%Psis(iB,iS,CA)%SubOps(iSSt)%Op)
         enddo
         deallocate(this%Psis(iB,iS,CA)%SubOps)
      endif
   enddo
   enddo
   enddo
   deallocate(this%Psis)
   

end subroutine

!===============================================================================
subroutine dest_TOperator(this)
!===============================================================================
   type(TOperator)                        :: this
!local
   integer                                :: iSSt

   if(.not.associated(this%Subops))return
   do iSSt=0,size(this%Subops)-1
      if(.not.associated(this%Subops(iSSt)%Op))cycle
      deallocate(this%SubOps(iSSt)%Op)
   enddo
   deallocate(this%Subops)
end subroutine dest_TOperator

!===============================================================================
subroutine u_hamiltonian(u_matrix, h, dstates, dpsis)
!===============================================================================
      !< \brief Add interaction term corresponding to U to the Hamiltonian

      type(toperator), intent(inout) :: h
      type(tstates), intent(in)      :: dstates
      type(tpsis), intent(in)        :: dpsis
      real(KINDR)                    :: u_matrix(:,:,:,:)

      integer :: sst, sst1, sst2, sst3, sst4
      integer :: b1, b2, b3, b4, sp1, sp2, sp3, sp4
      real(KINDR) :: u_value

      if(.not.associated(dpsis%psis)) stop '[u_hamiltonian] DAFUQ?'

      do b1 = 1,dstates%NBands
      do sp1 = 1,2
        do b2 = 1,dstates%NBands
        do sp2 = 1,2
          do b3 = 1,dstates%NBands
          do sp3 = 1,2
            do b4 = 1,dstates%NBands
            do sp4 = 1,2
               ! Check that U_1234 is non-zero (there is an interaction)
               u_value = 0.5*u_matrix(2*(b1-1)+sp1,2*(b2-1)+sp2,2*(b3-1)+sp3,2*(b4-1)+sp4)
               if(u_value /= 0) then
               do sst = 0,dstates%NSStates-1
                  ! Make the connections: check if annihilator corresponding
                  ! to c[b_4,s_4] connects to any block (sst4 >= 0), then
                  ! continue  with that superstate and check for the next
                  ! operators until a chain c^+ c^+ c c is built.
                  sst4 = DStates%SubStates(sst)%Connect(b3,sp3,2)
                  if (sst4 < 0) cycle
                  sst3 = DStates%SubStates(sst4)%Connect(b4,sp4,2)
                  if (sst3 < 0) cycle
                  sst2 = DStates%SubStates(sst3)%Connect(b2,sp2,1)
                  if (sst2 < 0) cycle
                  sst1 = DStates%SubStates(sst2)%Connect(b1,sp1,1)
                  if (sst1 < 0) cycle
                  if (sst /= sst1) then
                    cycle
!                 write (0,"('U(',3(2I1,','),2I1,')=',F10.2)") &
!                 b1, sp1, b2, sp2, b3, sp2, b4, sp1, u_value
!                 write (0,"('sst chain:',5(I5,'<'))") sst1, sst2, sst3, sst4, sst
!                 stop 'Hamiltonian not block-diagnonal'
                  endif

                  ! Now, compute the product c^+ c^+ c c acting on the
                  ! corresponding subspace
                  h%SubOps(sst)%Op = h%SubOps(sst)%Op + u_value * &
                           matmul(DPsis%Psis(b1,sp1,1)%SubOps(sst2)%Op, &
                           matmul(DPsis%Psis(b2,sp2,1)%SubOps(sst3)%Op, &
                           matmul(DPsis%Psis(b4,sp4,2)%SubOps(sst4)%Op, &
                                 DPsis%Psis(b3,sp3,2)%SubOps(sst)%Op)))
               enddo
               endif
               ! End U check
            enddo
          enddo
        enddo
        enddo
        enddo
        enddo
      enddo
      enddo
end subroutine

!===============================================================================
pure subroutine force_hermitian(H, DStates)
!===============================================================================
   type(toperator), intent(inout) :: H
   type(tstates), intent(in)      :: DStates

   integer :: sst, i, j

   do sst = 0, DStates%NSStates - 1
      do i = 0, size(H%SubOps(sst)%Op, 1) - 2
         do j = i + 1, size(H%SubOps(sst)%Op, 2) - 1
            H%SubOps(sst)%Op(i, j) = (H%SubOps(sst)%Op(i, j) + H%SubOps(sst)%Op(j, i))/2.0d0
            H%SubOps(sst)%Op(j, i) = H%SubOps(sst)%Op(i, j)
         end do
      end do
   end do
end subroutine force_hermitian

!===============================================================================
!! This subroutine transforms the PSIS to the basis that is defined by the 
!! transformation-matrix HEVectors (in this case the eigenbasis of local Hamiltonian).
subroutine transform_Psis(DStates,transformed_Psis)
!===============================================================================
!input
   type(TStates)              :: DStates
   type(TPsis)                :: transformed_Psis
   integer                    :: iB, iS, CA, isst1, isst2      

   !generate psis in occupation-number-basis
   call init_Psi(transformed_Psis,DStates)
   ! call print_Psi(transformed_Psis,DStates)
   ! this%Psis(iB,iS,iDag)%SubOps(iSSt1)%Op(iSt2,iSt1)

   do iB=1,size(transformed_Psis%Psis(:,1,1))      
   do iS=1,2
   do CA=1,2
      !write(unit=*,fmt=*) "ib is ca", ib, is, ca
      do isst1=0, DStates%NSStates-1
      do isst2=0, DStates%NSStates-1
         if (isst1.eq.isst2) cycle
         if (DStates%Substates(isst1)%Connect(iB,iS,CA).eq.isst2) then
         !write(unit=*,fmt=*) isst1, isst2
         
         !apply unitary basis transformation to PSIS:
         !          PSI'=U^+(new superstate) PSI U(old superstate)
!        transformed_Psis%Psis(ib,is,CA)%SubOps(isst1)%Op=&
!        matmul(transpose(HEVectors%SubOps(isst2)%op),&
!        matmul(transformed_Psis%Psis(ib,is,CA)%SubOps(isst1)%Op,HEVectors%SubOps(isst1)%op))

         transformed_Psis%Psis(ib,is,CA)%SubOps(isst1)%Op=&
         matmul(transpose(DStates%SubStates(isst2)%Evec),&
         matmul(transformed_Psis%Psis(ib,is,CA)%SubOps(isst1)%Op,DStates%Substates(isst1)%EVec))
         endif
      enddo
      enddo
   enddo
   enddo
   enddo

end subroutine transform_Psis


!===============================================================================
!! Helper function for debugging and testing; it transforms a given state_vector
!! in eigenbasis back to occupation-number-basis.
subroutine transform_state_eb_to_occbasis(DStates,state_vector,debug)
!===============================================================================

   type(TStates)              :: DStates
   real(KINDR), intent(inout) :: state_vector(0:DStates%NStates-1)
   integer :: isst

   ! integer, dimension(1) :: sh
   ! integer :: iii
   logical :: debug


   if (debug) write(unit=*,fmt=*) "==============transform_state_eb_to_occbasis=============="
   if (debug) then
      call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
   endif

   do isst=0,DStates%NSStates-1
      !if (debug) write(unit=*,fmt=*) "isst ", isst
      !if (debug)write(unit=*,fmt=*) "offset", DStates%SubStates(isst)%Offset
      !if (debug)write(unit=*,fmt=*) "nstates in substate", DStates%Substates(isst)%NStates

      ! apply unitary basis transformation:
      !       state_vector' = U state_vector
      state_vector(DStates%SubStates(isst)%Offset:DStates%SubStates(isst)%Offset+DStates%Substates(isst)%NStates-1)=&
      &matmul(DStates%substates(isst)%Evec,state_vector(DStates%SubStates(isst)%Offset:&
      DStates%SubStates(isst)%Offset+DStates%Substates(isst)%NStates-1))

   enddo

   if (debug) then  
      call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
   endif

end subroutine transform_state_eb_to_occbasis



!===============================================================================
!! Helper function for debugging and testing; it transforms a given state_vector
!! in occupation-number-basis to eigenbasis of local Hamiltonian.
subroutine transform_state_occbasis_to_eb(DStates,state_vector,debug)
!===============================================================================

   type(TStates)              :: DStates
   real(KINDR), intent(inout) :: state_vector(0:DStates%NStates-1)
   integer :: isst
   logical :: debug


   if (debug) write(unit=*,fmt=*) "==============transform_state_occbasis_to_eb=============="
   if (debug) then
      call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
   endif

   do isst=0,DStates%NSStates-1
   
!  if (debug) write(unit=*,fmt=*) "isst ", isst
!  if (debug)write(unit=*,fmt=*) "offset", DStates%SubStates(isst)%Offset
!  if (debug)write(unit=*,fmt=*) "nstates in substate", DStates%Substates(isst)%NStates
!  call ausgabe(transpose(HEVectors%SubOps(isst)%op),shape(transpose(HEVectors%SubOps(isst)%op)),"HEVectors%SubOps(isst)%op",3,.false.)

   ! apply unitary basis transformation:
   !       state_vector' = U^+ state_vector
   state_vector(DStates%SubStates(isst)%Offset:DStates%SubStates(isst)%Offset+DStates%Substates(isst)%NStates-1)=&
   &matmul(transpose(DStates%substates(isst)%Evec),state_vector(DStates%SubStates(isst)%Offset:&
   DStates%SubStates(isst)%Offset+DStates%Substates(isst)%NStates-1))

   enddo

   if (debug) then  
      call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
   endif

end subroutine transform_state_occbasis_to_eb


!===============================================================================
!! Helper function that applies a given operator Psis of kind (iB,iS,CA) to state_vector;
!! state_vector and Psis can be in whatever basis. 
!! Be careful: actual_superstate has to be given to function but it will be changed to
!! the label of target-superstate!
subroutine apply_psi_to_ket_EB(state_vector,actual_superstate,iB,IS,CA,DStates,debug)
!===============================================================================
   type(TStates) , intent(in)                  :: DStates
!  type(TPsis), intent(in)                     :: Psis
   real(KINDR), intent(inout) :: state_vector(0:DStates%NStates-1)
   integer :: iB, iS, CA
   integer :: connect
   integer, intent(inout) :: actual_superstate
   logical :: debug

   if (debug) write(unit=*,fmt=*) "==============apply_psi_to_ket=============="
   if (debug) then
      call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
   endif

   !check if operator of kind (iB,iS,CA) can lead to another superstate;
   !if not, set state_vector to zero
   connect=DStates%Substates(actual_superstate)%Connect(iB,iS,CA)
   if (connect .eq. -1) then
      state_vector=0
      if (debug) write(unit=*,fmt=*) "====>>>> not connected!"
      return
   endif

   if (debug) then
      write(unit=*,fmt=*) "ib, is, ca", ib, is, ca
      write(unit=*,fmt=*) "actual_superstate", actual_superstate
      write(unit=*,fmt=*) "connect ", connect
      write(unit=*,fmt=*) "offset", DStates%SubStates(connect)%Offset
      write(unit=*,fmt=*) "number of states in substate", DStates%SubStates(connect)%NStates
      write(unit=*,fmt=*) " "
      write(unit=*,fmt=*) " "
      write(unit=*,fmt=*) "actual superstate", state_vector(DStates%SubStates(actual_superstate)%Offset:&
      DStates%SubStates(actual_superstate)%Offset+DStates%SubStates(actual_superstate)%NStates-1)

      write(unit=*,fmt=*) "goal superstate", & 
      & state_vector(DStates%SubStates(connect)%Offset:DStates%SubStates(connect)%Offset+&
      & DStates%SubStates(connect)%NStates-1)
   endif

   !    state_vector' = Psi state_vector
!  call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
   state_vector(DStates%SubStates(connect)%Offset:DStates%SubStates(connect)%Offset+DStates%SubStates(connect)%NStates-1)=&
      &matmul(DStates%SubStates(actual_superstate)%Psis_EB(iB,iS,CA)%Op,&
      &state_vector(DStates%SubStates(actual_superstate)%Offset:DStates%SubStates(actual_superstate)%Offset+& 
      & DStates%SubStates(actual_superstate)%NStates-1))
!  call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)

! !      DStates%SubStates(actual_superstate)%Psis_EB(iB,iS,CA)%Op   

!  call ausgabe(DStates%SubStates(actual_superstate)%Psis_EB(iB,iS,CA)%Op,shape(DStates%SubStates(actual_superstate)%Psis_EB(iB,iS,CA)%Op),&
!  &"DStates%SubStates(actual_superstate)%Psis_EB(iB,iS,CA)%Op",3,.false.)
!  
!  call ausgabe(Psis%Psis(ib,is,ca)%SubOps(actual_superstate)%Op,shape(Psis%Psis(ib,is,ca)%SubOps(actual_superstate)%Op),"Psis%Psis(ib,is,ca)%SubOps(actual_superstate)%Op",3,.false.)
   ! set old entries of state_vector to zero
   state_vector(DStates%SubStates(actual_superstate)%Offset:DStates%SubStates(actual_superstate)%Offset+& 
   & DStates%SubStates(actual_superstate)%NStates-1)=0

   if (debug) write(unit=*,fmt=*) "modified goal superstate",&
   & state_vector(DStates%SubStates(connect)%Offset:DStates%SubStates(connect)%Offset+DStates%SubStates(connect)%NStates-1)
   
   ! save in which superstate we are now
   actual_superstate=connect

   if (debug) then
      call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
      write(unit=*,fmt=*) "-------------------------------------------------------------------------"
   endif

end subroutine apply_psi_to_ket_EB


!===============================================================================
!! Helper function that applies a given operator Psis of kind (iB,iS,CA) to state_vector;
!! state_vector and Psis can be in whatever basis. 
!! Be careful: actual_superstate has to be given to function but it will be changed to
!! the label of target-superstate!
subroutine apply_psi_to_ket_OCCB(state_vector,actual_superstate,iB,IS,CA,DStates,Psis,debug)
!===============================================================================
   type(TStates) , intent(in)                  :: DStates
   type(TPsis), intent(in)                     :: Psis
   real(KINDR), intent(inout) :: state_vector(0:DStates%NStates-1)
!   real(KINDR) :: state_vector_evolved(0:DStates%NStates-1)
   integer :: iB, iS, CA
   integer :: connect
   integer, intent(inout) :: actual_superstate
   logical :: debug

   if (debug) write(unit=*,fmt=*) "==============apply_psi_to_ket=============="
   if (debug) then
      call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
   endif

   !check if operator of kind (iB,iS,CA) can lead to another superstate;
   !if not, set state_vector to zero
   connect=DStates%Substates(actual_superstate)%Connect(iB,iS,CA)
   if (connect .eq. -1) then
      state_vector=0
      if (debug) write(unit=*,fmt=*) "====>>>> not connected!"
      return
   endif

   if (debug) then
      write(unit=*,fmt=*) "ib, is, ca", ib, is, ca
      write(unit=*,fmt=*) "actual_superstate", actual_superstate
      write(unit=*,fmt=*) "connect ", connect
      write(unit=*,fmt=*) "offset", DStates%SubStates(connect)%Offset
      write(unit=*,fmt=*) "number of states in substate", DStates%SubStates(connect)%NStates
      write(unit=*,fmt=*) " "
      write(unit=*,fmt=*) " "
      write(unit=*,fmt=*) "actual superstate", state_vector(DStates%SubStates(actual_superstate)%Offset:&
      DStates%SubStates(actual_superstate)%Offset+DStates%SubStates(actual_superstate)%NStates-1)

      write(unit=*,fmt=*) "goal superstate", state_vector(DStates%SubStates(connect)%Offset:DStates%SubStates(connect)%Offset+&
      & DStates%SubStates(connect)%NStates-1)
   endif

   !    state_vector' = Psi state_vector
   state_vector(DStates%SubStates(connect)%Offset:DStates%SubStates(connect)%Offset+DStates%SubStates(connect)%NStates-1)=&
      &matmul(Psis%Psis(ib,is,ca)%SubOps(actual_superstate)%Op,&
      &state_vector(DStates%SubStates(actual_superstate)%Offset:DStates%SubStates(actual_superstate)%Offset+&
      &DStates%SubStates(actual_superstate)%NStates-1))

!  if (DStates%SubStates(actual_superstate)%Connect(ib,is,ca)==-1) then
! !      write(unit=*,fmt=*) "not connected#!!!!!"
!     state_vector_evolved=0
!  else
   
!  if (DStates%SubStates(actual_superstate)%Connect(ib,is,ca)==-1) then
! !      write(unit=*,fmt=*) "not connected#!!!!!"
!  else
!      call SMatVecProd(DStates%SubStates(actual_superstate)%Psis(ib,is,ca),&
!     state_vector_evolved(DStates%SubStates(actual_superstate)%Offset:DStates%SubStates(actual_superstate)%Offset+DStates%SubStates(actual_superstate)%NStates-1),&
!     state_vector(DStates%SubStates(DStates%SubStates(actual_superstate)%Connect(ib,is,ca))%Offset:&
!     DStates%SubStates(DStates%SubStates(actual_superstate)%Connect(ib,is,ca))%Offset+DStates%SubStates(DStates%SubStates(actual_superstate)%Connect(ib,is,ca))%NStates-1))
!  endif
   
! !      DStates%SubStates(actual_superstate)%Psis_EB(iB,iS,CA)%Op   

!     call ausgabe(DStates%SubStates(actual_superstate)%Psis_EB(iB,iS,CA)%Op,shape(DStates%SubStates(actual_superstate)%Psis_EB(iB,iS,CA)%Op),&
!     &"DStates%SubStates(actual_superstate)%Psis_EB(iB,iS,CA)%Op",3,.false.)
!     
!     call ausgabe(Psis%Psis(ib,is,ca)%SubOps(actual_superstate)%Op,shape(Psis%Psis(ib,is,ca)%SubOps(actual_superstate)%Op),"Psis%Psis(ib,is,ca)%SubOps(actual_superstate)%Op",3,.false.)
   ! set old entries of state_vector to zero
   state_vector(DStates%SubStates(actual_superstate)%Offset:DStates%SubStates(actual_superstate)%Offset+&
   & DStates%SubStates(actual_superstate)%NStates-1)=0
   
   
   if (debug) write(unit=*,fmt=*) "modified goal superstate", &
   & state_vector(DStates%SubStates(connect)%Offset:DStates%SubStates(connect)%Offset+DStates%SubStates(connect)%NStates-1)
   
   ! save in which superstate we are now
   actual_superstate=connect

   if (debug) then
      call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
      write(unit=*,fmt=*) "-------------------------------------------------------------------------"
   endif

end subroutine apply_psi_to_ket_OCCB


!===============================================================================
!! Small helper function to look at Hamiltonian under unitary transformation.
subroutine test_hamiltonian(DStates,HEVectors,DHamiltonian)
!===============================================================================
   type(TStates) , intent(in)        :: DStates
   type(TOperator),intent(in)        :: HEVectors
   type(TOperator),intent(in)          ::  DHamiltonian
   integer :: isst

   do isst=0,DStates%NSStates-1
      write(unit=*,fmt=*) "isst", isst
      !print old hamiltonian
      call ausgabe(DHamiltonian%SubOps(isst)%Op,shape(DHamiltonian%SubOps(isst)%Op),"DHamiltonian%SubOps(isst)%Op",3,.false.)
      
      !transform Hamiltonian
      !      H' = EV^+ H EV
      DHamiltonian%SubOps(isst)%Op=&
            matmul(transpose(HEVectors%SubOps(isst)%op),&
            matmul(DHamiltonian%SubOps(isst)%Op,HEVectors%SubOps(isst)%op))
      
      !print new hamiltonian
      call ausgabe(DHamiltonian%SubOps(isst)%Op,shape(DHamiltonian%SubOps(isst)%Op),"DHamiltonian%SubOps(isst)%Op",3,.false.)
   enddo

end subroutine test_hamiltonian


!===============================================================================
end module MOperator
!===============================================================================

#ifdef Operator_Test

!===============================================================================
program Prog_Operator
!===============================================================================
use MParameters
use MStates
use MOperator
use MAusgabe

! integer, parameter :: nbands= 3
integer                                :: iB,iS,CA,iSSt,Offset
type(TStates)                          :: DStates
type(TPsis)                            :: DPsis, DTransformed_Psis
type(TOperator)                        :: DHamiltonian,DS2,HEVectors
real(KINDR),allocatable                :: HEValues(:),HS2Values(:)
!> the number of blocks in the hamiltonian = number of substates
integer                             :: nsubstates
!> saves the states belonging to one substate
integer, allocatable                :: states2substates(:)
real(KINDR), allocatable :: u_matrix(:,:,:,:)

real(KINDR), allocatable :: state_vector(:), state_vector_eb(:), original_state_vector(:)
integer :: iii, i, actual_superstate
integer, dimension(1) :: sh
integer :: i_anfangszustand, state_number_in_array_of_states
integer :: i_substate, SSt,isst1,isst2
logical :: debug =.false.

real(KINDR),allocatable                :: muimp(:,:)
real(KINDR),allocatable                :: temp(:)

call read_ParameterFile("Parameters.in")
call init_States(DStates)
call qns2substates(DStates, nsubstates, states2substates)
call init_SubStates(DStates, nsubstates, states2substates)

call u_allocate(DStates%NBands,u_matrix)
write(*,*) "DStates%NBands",DStates%NBands

!write(unit=*,fmt=*) "size u m", shape(u_matrix)

!Read in a U-Matrix from file, since it now is generated in pyhton-module
!u_matrix_2_band_kana.dat
!u_matrix_5_band_coulomb.dat

open(4711, file='u_matrix_2_band_kana.dat', status='old', action='read', access='sequential')
call u_readin(4711,u_matrix)
close(4711)

allocate(state_vector(0:DStates%NStates-1))
allocate(state_vector_eb(0:DStates%NStates-1))
allocate(original_state_vector(0:DStates%NStates-1))

!write(stdout,*)"================================================================"
!write(stdout,*)"             initialising psis"
!write(stdout,*)"================================================================"
call init_Psi(DPsis,DStates)
!write(stdout,*)"================================================================"
!write(stdout,*)"             printing psis"
!write(stdout,*)"================================================================"
! call print_Psi(DPsis,DStates)
!write(stdout,*)"================================================================"
!write(stdout,*)"             initialising Hamiltonian"
!write(stdout,*)"================================================================"
allocate(muimp(DStates%NBands,2))
muimp(:,:) = get_Real_Parameter("mu") 
!write(*,*) "muimp", muimp
call init_Hamiltonian(DHamiltonian,DStates,muimp)
write(stdout,*)"================================================================"
write(stdout,*)"             print S^2"
write(stdout,*)"================================================================"
call init_S2alt2(DS2,DStates,DPsis)
call print_Operator(ds2,DStates) 
!write(stdout,*)"================================================================"
!write(stdout,*)"             printing Hamiltonian"
!write(stdout,*)"================================================================"
!call print_Operator(DHamiltonian,DStates) 
!call print_Operator(DS2,DStates) 
call init_TOperator(HEVectors,DStates)
allocate(HEValues(0:DStates%NStates-1))
!allocate(HS2Values(0:DStates%NStates-1))
!allocate(temp(DStates%NStates))
!write(stdout,*)"================================================================"
!write(stdout,*)"             add interaction part to hamiltonian"
!write(stdout,*)"================================================================"
call u_hamiltonian(u_matrix,DHamiltonian,DStates,DPsis)
write(stdout,*)"================================================================"
write(stdout,*)"             printing Hamiltonian with interaction part"
write(stdout,*)"================================================================"
call print_Operator(DHamiltonian,DStates) 
!stop
write(stdout,*)"================================================================"
write(stdout,*)"             diagonalising Hamiltonian"
write(stdout,*)"================================================================"
call diag_Operator(DHamiltonian,DStates,HEVectors,HEValues)
!call diag_Operator(DS2,DStates,HEVectors,HS2Values)
!write(stdout,*)"================================================================"
!write(stdout,*)"             printing eigenvalues of Hamiltonian"
!write(stdout,*)"================================================================"

!do iSSt=0,DStates%NSStates-1
!write(stdout,*)"Superstate:",iSSt,"Nt=", int(get_Nt(DStates,DStates%SubStates(iSSt)%States(0))),&
!& HEValues(DStates%SubStates(iSSt)%offset:DStates%SubStates(iSSt)%offset+&
   !DStates%SubStates(iSSt)%NStates-1)
!enddo

write(stdout,*)"================================================================"
write(stdout,*)"             printing eigenvectors of Hamiltonian"
write(stdout,*)"================================================================"
call print_Operator(HEVectors,DStates) 
call set_HEvec(HEVectors,DStates)
write(stdout,*)"================================================================"
write(stdout,*)"             transforming psis"
write(stdout,*)"================================================================"
call transform_Psis(DStates,DTransformed_Psis)
write(stdout,*)"================================================================"
write(stdout,*)"             printing transformed psis"
write(stdout,*)"================================================================"
!call print_Psi(DTransformed_Psis,DStates)
!write(stdout,*)"================================================================"
!write(stdout,*)"            look at transformed hamiltonian"
!write(stdout,*)"================================================================"
!call test_hamiltonian(DStates,HEVectors,DHamiltonian)
!write(stdout,*)"================================================================"
!write(stdout,*)"            print psis as matrix"
!write(stdout,*)"================================================================"

! do iB=1,size(DTransformed_Psis%Psis(:,1,1))      
! do iS=1,2
! do CA=1,2
!  !write(unit=*,fmt=*) "ib is ca", ib, is, ca
!  do isst1=0, DStates%NSStates-1
!  do isst2=0, DStates%NSStates-1
!     if (isst1.eq.isst2) cycle
!     if (DStates%Substates(isst1)%Connect(iB,iS,CA).eq.isst2) then
!     write(unit=*,fmt=*) "sst ", isst1, " --> ", isst2
!     call ausgabe(DTransformed_Psis%Psis(iB,iS,CA)%SubOps(isst1)%Op,shape(DTransformed_Psis%Psis(iB,iS,CA)%SubOps(isst1)%Op),"DTransformed_Psis%Psis(iB,iS,CA)%SubOps(isst1)%Op",3,.false.)
!     endif
!  enddo
!  enddo
! enddo
! enddo
! enddo

write(stdout,*)"================================================================"
write(stdout,*)"             write transformed psis in DStates"
write(stdout,*)"================================================================"

call set_EB_Psis(DTransformed_Psis,DStates)

write(stdout,*)"================================================================"
write(stdout,*)"             test one transformed psis"
write(stdout,*)"================================================================"



!!!!! Samller test-function that applies one PSI of kind (iB,iS,CA) to a predefined ket
!debug=.false.
!state_number_in_array_of_states=6
!state_vector=0
!state_vector(state_number_in_array_of_states)=1
!actual_superstate=DStates%StatesSStates(state_number_in_array_of_states,2)
!iB=2
!iS=1
!CA=1

!write(*,*) "state_number", state_number_in_array_of_states
!write(unit=*,fmt=*) "connect", DStates%SubStates(actual_superstate)%Connect(iB,iS,CA)

!write(unit=*,fmt=*) "original state"
!call ausgabe(state_vector,shape(state_vector),"state_vector",10,.false.)

!call transform_state_occbasis_to_eb(DStates,state_vector,debug)

!write(unit=*,fmt=*) "transformed state in EB"
!call ausgabe(state_vector,shape(state_vector),"state_vector",10,.false.)

!call apply_psi_to_ket_EB(state_vector,actual_superstate,iB,iS,CA,DStates,debug)

!write(unit=*,fmt=*) "applied psi"
!call ausgabe(state_vector,shape(state_vector),"state_vector",10,.false.)

!call transform_state_eb_to_occbasis(DStates,state_vector,debug)

!write(unit=*,fmt=*) "transformed back to OCCB, psi applied"
!call ausgabe(state_vector,shape(state_vector),"state_vector",10,.false.)

!write(unit=*,fmt=*) "=========================="

!state_vector=0
!state_vector=0
!state_vector(state_number_in_array_of_states)=1
!actual_superstate=DStates%StatesSStates(state_number_in_array_of_states,2)

!write(unit=*,fmt=*) "original state" 
!call ausgabe(state_vector,shape(state_vector),"state_vector",10,.false.)

!call apply_psi_to_ket_OCCB(state_vector,actual_superstate,iB,iS,CA,DStates,DPsis,debug)

!write(unit=*,fmt=*) "transformed state in OCCB"
!call ausgabe(state_vector,shape(state_vector),"state_vector",10,.false.)


write(stdout,*)"================================================================"
write(stdout,*)"             test all transformed psis"
write(stdout,*)"================================================================"

!!!!! Test-function that goes through all possible states, iB, iS, CA.
!!!!! It compares the results obtained in occupation-number-basis and eigenbasis of local Hamiltonian
do i_anfangszustand=0,DStates%NStates-1
   do iB=1,DStates%NBands
   do iS=1,2
   do CA=1,2

      !define superstate one starts with
      actual_superstate=DStates%StatesSStates(i_anfangszustand,2)

      !define state_vector in occupation-number-basis
      state_vector=0
      state_vector_eb=0
      state_vector(i_anfangszustand)=1
      state_vector_eb(i_anfangszustand)=1
      original_state_vector=state_vector

!     write(unit=*,fmt=*) "i_anfangszustand", i_anfangszustand

      !apply psi to ket in occupation-number-basis
      call apply_psi_to_ket_OCCB(state_vector,actual_superstate,iB,iS,CA,DStates,dpsis,debug)
!     write(unit=*,fmt=*) "hier!!!!"

      !reset actual_superstate, because it is changed by subroutine apply_psi_to_ket
      actual_superstate=DStates%StatesSStates(i_anfangszustand,2)
      !apply transformed psi to transformed state_vector and transform state_vector back to occupation-number-basis
      call transform_state_occbasis_to_eb(DStates,state_vector_eb,debug)
      call apply_psi_to_ket_EB(state_vector_eb,actual_superstate,iB,iS,CA,DStates,debug)
      call transform_state_eb_to_occbasis(DStates,state_vector_eb,debug)

      !compare the two results
      if (sum(abs(state_vector-state_vector_eb)).ge.1.0e-10) then
         write(unit=*,fmt=*) "FAIL!!!!!!!!!!"
         write(unit=*,fmt=*) "ib is ca", ib, is, ca      
         call ausgabe(original_state_vector,shape(original_state_vector),"original_state_vector",30,.false.)
         call ausgabe(state_vector,shape(state_vector),"state_vector",30,.false.)
         call ausgabe(state_vector_eb,shape(state_vector_eb),"state_vector_eb",30,.false.)
         stop     
      else
!        write(unit=*,fmt=*) "OK!"
      endif

   enddo
   enddo
   enddo

enddo

write(unit=*,fmt=*) "All Psis OK!"


write(stdout,*)"================================================================"
write(stdout,*)"             test writing of psis in DStates"
write(stdout,*)"================================================================"



! do iSSt=0,DStates%NSStates-1
! do iB=1,DStates%NBands
! do iS=1,2
! do CA=1,2
!    SSt=DStates%SubStates(iSSt)%Connect(iB,iS,CA)
!    if(SSt.ne.-1)then
!     write(unit=*,fmt=*) "sst ", iSSt, " --> ", SSt
!     write(unit=*,fmt=*) "ib", ib
!     write(unit=*,fmt=*) "is", is
!     write(unit=*,fmt=*) "CA", CA
!     call ausgabe(DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op,shape(DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op),"DStates%SubStates(iSSt)%Psis_EB(iB,iS,CA)%Op",3,.false.)
! 
!    endif
! enddo
! enddo
! enddo
! enddo


! stop

! ! ! write(stdout,*)"================================================================"
! ! ! write(stdout,*)"             setting sparse Psis"
! ! ! write(stdout,*)"================================================================"
! ! ! call set_Psis(DPsis,DStates)
! ! ! write(stdout,*)"================================================================"
! ! ! write(stdout,*)"             destructing Psis"
! ! ! write(stdout,*)"================================================================"
call dest_Psis(DPsis)
call dest_Psis(DTransformed_Psis)
! ! ! write(stdout,*)"================================================================"
! ! ! write(stdout,*)"             setting sparse Hamiltonian"
! ! ! write(stdout,*)"================================================================"
! ! ! call set_Hamiltonian(DHamiltonian,DStates)
call dest_TOperator(DHamiltonian)
! ! ! write(stdout,*)"================================================================"
! ! ! write(stdout,*)"             setting Hamiltonian eigenvalues"
! ! ! write(stdout,*)"================================================================"
! ! ! call set_HEVal(HEValues,DStates)
! ! ! deallocate(HEValues)
! ! ! write(stdout,*)"================================================================"
! ! ! write(stdout,*)"             setting Hamiltonian eigenvectors"
! ! ! write(stdout,*)"================================================================"
! ! ! call set_HEVec(HEVectors,DStates)
! ! ! write(stdout,*)"================================================================"
! ! ! write(stdout,*)"           destructing Hamiltonian"
! ! ! write(stdout,*)"================================================================"
call dest_TOperator(HEVectors)
call dest_States(DStates)
call dest_Parameters()
end program Prog_Operator
#endif
