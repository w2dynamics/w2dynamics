!> This module is responsible for generating operators at random times and
!! calculating the weight of this configuration, i.e. calulating the trace over
!! local operators.
!!
!! @ author Nicolaus Parragh
!!
!! Module for calculating the trace in the CTQMC. The operators are stored in
!! a doubly-linked list to be able to compute the trace from left to right or
!! the other way round. 
!! We need an array to flag states which need to be considered in the trace due
!! to QNs. iSSt and EVec can be looked up in DStates.
!===============================================================================
module MTrace
!===============================================================================
use MParameters
use MStates
use MRandom
use MSparseMatrix
use MLanczos
use MMatrixUpdate
use MAusgabe
use MOperator
use MCompoundIndex

!> Type which stores informations about an operator, e.g. time and orbital
!! character. It is an element of a doubly linked list. 
   type :: TOper
!> Pointers to the next and previous operators. The operators are insterted into
!! the list in a time ordered fashion so next and prev point to operators at
!! earlier or later times.
      type(TOper),pointer              :: next=>null(),prev=>null()
      type(TOper),pointer              :: fnext=>null(),fprev=>null()
!> tau stores the time at which the hybridization event occurs.
      real(KINDR)                      :: tau
!> The three variables store the orbital and spin character of the hybridization
!! and specifies if an electron hops onto this site (Creation = 1) or hops into
!! the bath (Annihiliation = 2).
      integer                          :: orbital, spin, ca
!!> These arrays store the current state of the trace at time tau
!! (after the last accepted move) and possibly the state that will
!! replace it if the current move is accepted. Whether the state is a
!! bra or a ket depends on the position of the operator and is
!! indicated by the following logical values. The states are always
!! stored after application of the operator, i.e. at tau+ for kets and
!! tau- for bras.  cache_written indicates that the cache is to be
!! substituted for the state if the current move is accepted.
      real(KINDR), allocatable         :: state(:, :), cache(:, :)
      real(KINDR), allocatable         :: slogmax(:), clogmax(:)
      logical                          :: ket_state=.true., ket_cache=.true., cache_written=.false.
!> We need to renomaralize the states at each step since the time-evolution in
!! imaginary time is not unitarian. 
!     real(KINDR),allocatable          :: norml(:), normr(:)
!> Stores the information in which super state we are currently in 
!     integer,allocatable              :: sstl(:), sstr(:); sstl -> preceding_sst, sstr -> following_sst
!> Store/cache the superstate after and before this operator (from tau+ on/ up to tau-)
      integer                          :: following_sst=-2, preceding_sst=-2
      integer                          :: cache_post_sst=-2, cache_pre_sst=-2
!> Flag: from this operator on the trace has to be recalculated
      logical                          :: calc_new
!> Flag: stores if operator is connected to hyb or a worm operator
      logical                          :: has_hyb
   end type TOper

!> A helper type to be able to store vectors of different length in one array.
   type :: TSubArray
      real(KINDR),pointer              :: Vec(:)=>null()
   end type TSubArray

   type :: alloc_vector
      ! good, because it manages its memory itself
      real(KINDR), allocatable :: v(:)
   end type

!> A helper type to be able to store matrices of different length in one array.
   type :: TSubMatrix
!     real(KINDR)                      :: Det = 1.0
      real(KINDR),pointer              :: Mat(:,:)=>null()
   end type TSubMatrix

   type :: TLogTr
      real(KINDR) :: log, sign
   end type TLogTr

   interface operator(*)
      module procedure multiplyLogTr
   end interface operator(*)

   interface operator(/)
      module procedure divideLogTr
   end interface operator(/)

   interface operator(+)
      module procedure addLogTr
   end interface operator(+)

!> A helper type to be able to create arrays of operators which can then be
!! inserted and removed from the list.
   type :: TOperPointer
      type(TOper),pointer              :: p=>null()
   end type TOperPointer 

!> The main type of the trace module as in we use an object based approach. 
   type :: TTrace
!> The ket, bra, and superstates at beta half and whether they need to be recalculated.
      real(KINDR),pointer              :: ket_b2(:, :) => null(), ket_b2_logmax(:) => null()
      real(KINDR),pointer              :: bra_b2(:, :) => null(), bra_b2_logmax(:) => null()
      integer                          :: sst_ket_b2
      integer                          :: sst_bra_b2
      logical                          :: ket_b2_invalid=.true., bra_b2_invalid=.true.

!> The pointers to the beginning and end of the doubly linked list
      type(TOper),pointer              :: first=>null(),last=>null()
!> Pointers to the last operator before and first after the current/previous window
      type(TOper),pointer              :: prewin=>null(), postwin=>null()
      type(TOperPointer), allocatable  :: fprewin(:, :)
!> Pointer to the last operator before the window when the states
!  stored in the operators were last updated and pointers for the
!  range of operators that was potentially recalculated during the
!  last call to get_trace_eb
      type(TOper),pointer              :: last_prewin=>null()
      type(TOper),pointer              :: UpdateR=>null(), UpdateL=>null(), LastEndR=>null()

!> array of pointers to worm operators
!  we store in the order c cdag c cdag ...
      type(TOperPointer), allocatable  :: wormContainer(:)
      
!> The inverse temperature beta and the current value of the trace and the
!! bosonic trace initialised to one.
      real(KINDR)                      :: beta=1d0, cfgsign = -1.0_KINDR, BosonicTrace=1d0
      type(TLogTr)                     :: Trace = TLogTr(log=0.0_KINDR, sign=1.0_KINDR)
      type(TLogDet)                    :: Det = TLogDet(log=0.0_KINDR, sign=1.0_KINDR)

!!equal time offset to assure time-order coincides with linked lists
      real(KINDR)                      :: equal_time_offset=1d0

!> The number of operators in the trace.
      integer                          :: NOper=0
!> The maximum tau difference between operators for which a pair insertion or removal
!  may be proposed. This value is read in from the configuration file.
      real(KINDR)                      :: taudiff_max=5._KINDR
!> The lower and upper edge of the sliding window in the current position and in
!  the position to which the currently stored edge states belong.
      real(KINDR)                      :: tauwin_min, tauwin_max
!> The number of positions for the sliding window (not the divisor
!  NWin for the window width as in width = beta/NWin) and the current
!  window position expressed as number of half window widths from tau
!  = 0 to the beginning of the window
      integer                          :: num_windows, window_position=0
!> Number of operator pairs inside the window that are eligible for
!> removal. ensure_valid_NPairs must be called before using this value
!> if it is supposed to be correct.
      integer                          :: NPairs
!> Number of operators per flavour before the window, validity of the entries and the
!> last operator included in the count. ensure_valid_PreWinFCount must be called before
!> using the entries if they are supposed to be correct.
      integer,allocatable              :: PreWinFCount(:,:,:)
      logical                          :: PreWinFCountValid=.false.
      type(TOper), pointer             :: FCountLastIncluded => null()
!> The number of creation/annihilation operators in the trace. Orbital and spin
!! resolved.
!! this now saves the operators individually, and not in pairs as originally
      integer,pointer                  :: NOSOper(:,:)=>null()
      integer,pointer                  :: NOSCAOper(:,:,:)=>null()
!> The hybridization function orbital and spin resolved. Last index is tau.
      real(KINDR),pointer              :: HybrF(:,:,:)=>null()
      real(KINDR),pointer              :: HybrF_full(:,:,:,:,:)=>null()
!> The screening function orbital and spin resolved. Last index is tau.
      real(KINDR),pointer              :: screening_function(:,:,:,:,:)=>null()
!> The inverted matrix of the hybridization lines.
      type(TSubMatrix),pointer         :: MInv(:,:)=>null()
      type(TSubMatrix),pointer         :: MInv_full=>null()
!> The global update.
      integer,pointer                  :: gu(:)=>null()
!> The States array contains the necessary information to know which states need
!! to be considered when calculating the local trace. The format is:
!! first index in States is: 
!! 1: SuperState which contributes to trace
!! 2: number of states in this superstate which contribute
      integer,allocatable              :: States(:,:)
!> Lookup array for the indices of superstates in States.
!! (sst_to_statesindex(States(x, 1)) = x)
      integer, allocatable             :: sst_to_statesindex(:)
!> Store the current/target and source outer state, superstate and its size
!  - outer_sst may take on values in 0..DStates%NSStates-1 for which
!    at least one state is included in the outer truncation (if
!    enabled) and the superstate it represents is
!    DStates%SubStates(outer_sst)
!  - If state sampling is enabled, outer_state may take on values in
!    1..DTrace%States(DTrace%sst_to_statesindex(outer_sst), 2), which
!    is equal to 1..DStates%SubStates(outer_sst)%NStates when no outer
!    truncation is applied, and represents the state
!    DStates%SubStates(outer_sst)%States(outer_state - 1).
!    If state sampling is not enabled, outer_state and outer_state_old
!    MUST always be set to 1 for convenience
!  - outer_sst_size is the (potentially truncated) size of the
!    superstate represented by outer_sst (between moves incl. at
!    measurement time guaranteed, otherwise watch out for update
!    order).
!    For convenience, it MUST be set to 1 if state sampling
!    is enabled
      integer                          :: outer_sst, outer_sst_old, outer_sst_size
      integer                          :: outer_state, outer_state_old
      logical                          :: b_statesampling, b_offdiag, b_fix_prealloc_stvec

!> The number of lowest energy states which we consider in the trace.
      integer                          :: NTruncStates, NTruncStatesMax
!> hack: make NStatesMax available without DStates for b_fix_prealloc_stvec
      integer                          :: NStatesMax
!> The ground state energy of the system in the atomic limit.
      real(KINDR)                      :: Egs
! in this array the state at beta half is stored every time the local trace is calculated
!! This information is necessary to e.g. calculate the density matrix.
      type(TLogTr), allocatable        :: parttrace(:),tparttrace(:)
!     real(KINDR),allocatable          :: parttrace_oper_resolved(:,:)
!     real(KINDR),allocatable          :: parttrace_oper_resolved_summed(:,:)

!> The operator that is nearest to, but left of (tau greater than) beta/2.
!! This definition "wraps around" the trace, so if all operators are to the
!! right of beta/2, it instead takes the rightmost (smallest tau) operator.
!! This is used to be able to measure time-independent quantities to the right
!! of this operator.
      type(TOper),pointer              :: beta_half_operator

!> The Lanczos histogram. Currently not stored or handed back anywhere.
      integer                          :: lhisto(PMAX)
!> The non-reversibility histogram
      integer                          :: rhisto(PMAX)
!> epsilon criterion for amplitudes in the states to be numerically zero
!! EPSTRACEVEC in the parameter file
      real(KINDR)                      :: epstracevec
!> epsilon criterion for the sandwich of the trace to be numerically zero
!! EPSSANDWICH in the parameter file
      real(KINDR)                      :: epssandwich
!> Flag set if we are doing a paramagnetic calculation.
      integer                          :: ParaMag
!> The number of symmtry moves specified in the paramter file.
      integer                          :: NSymMove
!> The symmetry moves specified in the paramter file.
      integer,allocatable              :: SymMoves(:,:)
!> Flag set to one if we perform a calculation with phonons.
      integer                          :: Phonon
!> The number of operators in a pool we can recycle.
      integer                          :: iOperPool
!> This array stores old operators for recycling.
      type(TOperPointer)               :: OperPool(10000)

      type(alloc_vector), allocatable  :: tau_c(:,:)  !< minv creation op taus
      type(alloc_vector), allocatable  :: tau_a(:,:)  !< minv annihilation op taus
!      type(alloc_vector), allocatable  :: urho(:,:)  !< stores U x rho^(1)
!> types for segment
      !!! index: band,spin
      logical,allocatable :: outerstate_tmp(:,:)
      logical,allocatable :: initial_outerstate(:,:)
      logical,allocatable :: initial_outerstate_old(:,:)
      real(kindr),allocatable :: mu_accumulator(:,:)
      real(kindr),allocatable :: mu_accumulator_tmp(:,:)
      real(kindr),allocatable :: int_accumulator(:,:,:,:)
      real(kindr),allocatable :: int_accumulator_tmp(:,:,:,:)
      !!! u-matrix and mumip
      real(kindr),allocatable :: u_matrix(:,:,:,:)
      real(kindr),allocatable :: u_matrix_dd(:,:,:,:)
      real(kindr),allocatable :: muimp(:,:)

      type(TOperPointer), dimension(:,:), allocatable :: ffirst
      type(TOperPointer), dimension(:,:), allocatable :: flast

      logical :: outerstate_full
      logical :: outerstate_empty

!> special types for calculation of the full trace in segment
   logical,allocatable :: actual_state_full(:,:)
   logical,allocatable :: has_operator(:,:)
   real(kindr),allocatable :: mu_accumulator_full_tmp(:,:)
   real(kindr),allocatable :: mu_accumulator_full(:,:)
   real(kindr),allocatable :: int_accumulator_full_tmp(:,:,:,:)
   real(kindr),allocatable :: int_accumulator_full(:,:,:,:)

   end type TTrace

   type :: TTrace_pointer
      type(TTrace),pointer :: ptr
   end type TTrace_pointer

contains

pure elemental function multiplyLogTr(tr1, tr2) result(tr3)
   type(TLogTr), intent(in)  :: tr1, tr2
   type(TLogTr)              :: tr3
   tr3 = TLogTr(log = tr1%log + tr2%log, sign = tr1%sign * tr2%sign)
end function multiplyLogTr

pure elemental function divideLogTr(tr1, tr2) result(tr3)
   type(TLogTr), intent(in)  :: tr1, tr2
   type(TLogTr)              :: tr3
   tr3 = TLogTr(log = tr1%log - tr2%log, sign = tr1%sign * tr2%sign)
end function divideLogTr

pure elemental function addLogTr(tr1, tr2) result(tr3)
   type(TLogTr), intent(in), target  :: tr1, tr2
   type(TLogTr)                      :: tr3

   if (tr1%log >= tr2%log) then
      tr3 = TLogTr(log = tr1%log + log(1.0_KINDR + tr1%sign * tr2%sign * exp(tr2%log - tr1%log)),&
                   sign = tr1%sign)
   else
      tr3 = TLogTr(log = tr2%log + log(1.0_KINDR + tr1%sign * tr2%sign * exp(tr1%log - tr2%log)),&
                   sign = tr2%sign)
   end if
end function addLogTr

pure elemental real(KINDR) function trval(tr)
   type(TLogTr), intent(in) :: tr
   trval = tr%sign * exp(tr%log)
end function trval

pure elemental real(KINDR) function wrat(trrat, detrat)
! UNSIGNED weight ratio
   type(TLogTr), intent(in)  :: trrat
   type(TLogDet), intent(in) :: detrat
   wrat = exp(trrat%log + detrat%log)
end function wrat

! flavor comparison function
! op1: pointer to operator
! op2: pointer to operator
!
! returns: true if flavors are equal, false otherwise
pure logical function flaveq(op1, op2)
   type(TOper), pointer, intent(in) :: op1, op2

   flaveq = op1%Orbital == op2%Orbital .and. op1%Spin == op2%Spin
end function flaveq


!===============================================================================
!> The following subroutine initialises the type trace. It has to do a lot more
!! than the equivalent functions in other modules since we need to compute some
!! more information about the system before we can start the simulation, e.g.
!! find all the states we want to consider in the outer trace.
subroutine init_Trace(this,DStates,FTau,FTau_full,screening_function,Nftau,muimp,u_matrix)
!===============================================================================
   type(TTrace)                        :: this
!input
!> The initialised states type already containing the initialised superstates,
!! the local hamiltonian and also  its eigenvalues and eigenvectors.
   type(TStates),intent(in)            :: DStates
!> The number of tau points in the hybridization function array.
   integer                             :: Nftau
!> The discretized hybridization function diagonal in both orbit and spin.
   real(KINDR)                         :: FTau(DStates%NBands,2,Nftau)
   real(KINDR)                         :: FTau_full(DStates%NBands,2,DStates%NBands,2,Nftau)
   real(KINDR)                         :: screening_function(DStates%NBands,2,&
                                            DStates%NBands,2,Nftau)
!local
   integer                             :: iB,iS,i,iSt,j,iSSt,TNSSt
   integer                             :: Nd
   real(KINDR),allocatable             :: energies(:), emin(:), emax(:)
   real(KINDR),allocatable             :: empty_trace_values(:)
   real(KINDR)                         :: accumulator, rand_outer
   integer,allocatable                 :: degen(:)
   real(KINDR)                         :: dummy
   character(len=9)                    :: IDSymMove
   character(len=2)                    :: ISymMove
   real(KINDR)                :: u_matrix(2*DStates%NBands,2*DStates%NBands,2*DStates%NBands,2*DStates%NBands)
   real(KINDR)                :: muimp(DStates%NBands,2)
   integer :: b1,sp1,b2,sp2

! In here we cycle over the EVal array and look for the lowest energy states.
! We can then determine the degeneracy and allocate enough space in the States
! array. This array will be as big as the total amount of states considered in 
! the trace. The elements will be the number of the state considered and in the
! second column its iSSt value for each cycle (-1 if the trace yields 0 due to QNs
! or the superstate the state belongs to.).
   this%beta=get_Real_Parameter("beta")
   this%epstracevec=get_Real_Parameter("EPSTRACEVEC")
   this%epssandwich=get_Real_Parameter("EPSSANDWICH")
   if (get_Integer_Parameter("segment")==0) then
      if (get_Integer_Parameter("statesampling") /= 0) then
         this%b_statesampling = .true.
      else
         this%b_statesampling = .false.
      end if
   else
      this%b_statesampling = .false.
   end if
   if (get_Integer_Parameter("offdiag") == 0) then
      this%b_offdiag = .false.
   else
      this%b_offdiag = .true.
   endif

   if (get_Integer_Parameter("FixedPreallocatedVecs") == 0) then
      this%b_fix_prealloc_stvec = .false.
   else
      this%b_fix_prealloc_stvec = .true.
   endif

   this%equal_time_offset=this%beta*1d-14

   allocate(this%HybrF(DStates%NBands,2,Nftau))
   allocate(this%HybrF_full(DStates%NBands,2,DStates%NBands,2,Nftau))
   allocate(this%screening_function(DStates%NBands,2,DStates%NBands,2,Nftau))
   this%NSymMove=get_Integer_Parameter("NSymMove")
   allocate(this%SymMoves(2 * this%NSymMove, 2 * DStates%NBands))
   !!! allocate segement stuff
   Nd = get_Integer_Parameter("Nd")
   allocate(this%outerstate_tmp(Nd,2))
   allocate(this%initial_outerstate(Nd,2))
   this%initial_outerstate(:,:)=.false.
   allocate(this%initial_outerstate_old(Nd,2))
   this%initial_outerstate_old(:,:)=.false.
   allocate(this%mu_accumulator(Nd,2))
   allocate(this%mu_accumulator_tmp(Nd,2))
   allocate(this%int_accumulator(Nd,2,Nd,2))
   allocate(this%u_matrix(2*Nd,2*Nd,2*Nd,2*Nd))
   allocate(this%u_matrix_dd(Nd,2,Nd,2))
   allocate(this%muimp(Nd,2))
   allocate(this%int_accumulator_tmp(Nd,2,Nd,2))
   !!! segment stuff for full calc of trace
   allocate(this%has_operator(Nd,2))
   allocate(this%mu_accumulator_full(Nd,2))
   allocate(this%mu_accumulator_full_tmp(Nd,2))
   allocate(this%int_accumulator_full(Nd,2,Nd,2))
   allocate(this%int_accumulator_full_tmp(Nd,2,Nd,2))
   allocate(this%actual_state_full(Nd,2))
   !allocate(this%fehlende_flavour(Nd,2))

   !673 !    Density-density part of matrix
   !674         do i = 1, nband
   !675          do j = 1, nband
   !676            uumat(i,  j+nband,iat) = real(uu(i,j,i,j))
   !677            uumat(i+nband,j,iat)   = real(uu(i,j,i,j))
   !678            uumat(i,j,iat)     = real(uu(i,j,i,j)) - real(uu(i,j,j,i))
   !680            uumat(i+nband,j+nband,iat) = real(uu(i,j,i,j)) - real(uu(i,j,j,i))
   !682          end do
   !683         end do

   this%u_matrix(:,:,:,:)=u_matrix(:,:,:,:)

   this%u_matrix_dd(:,:,:,:)=0d0
   do b1 = 1,dstates%NBands
   do sp1 = 1,2
   do b2 = 1,dstates%NBands
   do sp2 = 1,2
      !write(*,*) "b1,sp1,b2,sp2,b3,sp3,b4,sp4 ", b1,sp1,b2,sp2,b3,sp3,b4,sp4 
      !write(*,*) "u_matrix(2*(b1-1)+sp1,2*(b2-1)+sp2,2*(b3-1)+sp3,2*(b4-1)+sp4) ", &
      !& u_matrix(2*(b1-1)+sp1,2*(b2-1)+sp2,2*(b3-1)+sp3,2*(b4-1)+sp4)
      if(b1.eq.b2)then
         this%u_matrix_dd(b1,sp1,b2,sp2)=u_matrix(2*(b1-1)+sp1,2*(b2-1)+sp2,2*(b1-1)+sp1,2*(b2-1)+sp2)
      endif
      if(b1.ne.b2)then
         if(sp1.eq.sp2)then
            this%u_matrix_dd(b1,sp1,b2,sp2)=u_matrix(2*(b1-1)+sp1,2*(b2-1)+sp2,2*(b1-1)+sp1,2*(b2-1)+sp2)-&
            &                               u_matrix(2*(b1-1)+sp1,2*(b2-1)+sp2,2*(b2-1)+sp2,2*(b1-1)+sp1)
         else
            this%u_matrix_dd(b1,sp1,b2,sp2)=u_matrix(2*(b1-1)+sp1,2*(b2-1)+sp2,2*(b1-1)+sp1,2*(b2-1)+sp2)
         endif
      endif
   enddo
   enddo
   enddo
   enddo

   this%muimp(:,:)=muimp(:,:)
   allocate(this%ffirst(dstates%nbands,2))
   allocate(this%flast(dstates%nbands,2))

   this%outerstate_full=.false.
   this%outerstate_empty=.false.
   
   do i=1,this%NSymMove
      write(ISymMove,'(I2.2)')i
      IDSymMove="SymMove"//ISymMove
      call get_Integer_List(IDSymMove,this%SymMoves(i,:))
      do j = 1, 2 * DStates%NBands
         this%SymMoves(this%NSymMove + i, this%SymMoves(i, j)) = j
      end do
   enddo

   this%NStatesMax = DStates%NStatesMax

   allocate(this%sst_to_statesindex(0:DStates%NSStates-1))
   this%sst_to_statesindex = -1
   if (get_String_Parameter("EnableOuterTruncation") == "YES") then
      allocate(energies(get_Integer_Parameter("truncation")))
      allocate(emin(get_Integer_Parameter("truncation")))
      allocate(emax(get_Integer_Parameter("truncation")))
      allocate(degen(get_Integer_Parameter("truncation")))
      energies(:)=huge(dummy)
      emin(:) = huge(emin(1))
      emax(:) = -huge(dummy)
      degen(:)=0

      do iSSt=0,DStates%NSStates-1
         ! generate list of eigenenergies included in outer truncation
         do iSt=0,DStates%SubStates(iSSt)%NStates-1
            do j=1,size(energies(:))
               if (abs(DStates%SubStates(iSSt)%EVal(iSt)-emin(j)).le.EPSDEG.or.&
                   abs(DStates%SubStates(iSSt)%EVal(iSt)-emax(j)).le.EPSDEG.or.&
                   DStates%SubStates(iSSt)%EVal(iSt).lt.emax(j).and.&
                   DStates%SubStates(iSSt)%EVal(iSt).gt.emin(j)) then
                  degen(j)=degen(j)+1
                  emin(j) = min(emin(j),DStates%SubStates(iSSt)%EVal(iSt))
                  emax(j) = max(emax(j),DStates%SubStates(iSSt)%EVal(iSt))
                  exit
               else if (DStates%SubStates(iSSt)%EVal(iSt).lt.emin(j)-EPSDEG) then
                  do i=size(energies(:))-1,j,-1
                     energies(i+1)=energies(i)
                     emin(i+1) = emin(i)
                     emax(i+1) = emax(i)
                     degen(i+1)=degen(i)
                  enddo
                  energies(j)=DStates%SubStates(iSSt)%EVal(iSt)
                  emin(j) = DStates%SubStates(iSSt)%EVal(iSt)
                  emax(j) = DStates%SubStates(iSSt)%EVal(iSt)
                  degen(j)=1
                  exit
               endif
            enddo
         enddo
      enddo
      this%Egs=energies(1)

      ! count superstates contributing to outer truncation
      TNSSt=0
      do iSSt=0,DStates%NSStates-1
         do j=1,size(energies)
            if (abs(DStates%SubStates(iSSt)%EVal(0)-emin(j)).le.EPSDEG.or.&
                abs(DStates%SubStates(iSSt)%EVal(0)-emax(j)).le.EPSDEG.or.&
                DStates%SubStates(iSSt)%EVal(0).lt.emax(j).and.&
                DStates%SubStates(iSSt)%EVal(0).gt.emin(j)) then
               TNSSt=TNSSt+1
               exit
            endif
         enddo
      enddo

      ! generate list of contributing superstates and the numbers of
      ! contributing states in each one
      allocate(this%States(TNSSt,2))
      this%States=0
      TNSSt=0
      do iSSt=0,DStates%NSStates-1
         do j=1,size(energies)
            if (abs(DStates%SubStates(iSSt)%EVal(0)-emin(j)).le.EPSDEG.or.&
                abs(DStates%SubStates(iSSt)%EVal(0)-emax(j)).le.EPSDEG.or.&
                DStates%SubStates(iSSt)%EVal(0).lt.emax(j).and.&
                DStates%SubStates(iSSt)%EVal(0).gt.emin(j)) then
               TNSSt=TNSSt+1
               exit
            endif
         enddo

         do j=1,size(energies)
            do i=0,size(DStates%SubStates(iSSt)%EVal)-1
               if (abs(DStates%SubStates(iSSt)%EVal(i)-emin(j)).le.EPSDEG.or.&
                   abs(DStates%SubStates(iSSt)%EVal(i)-emax(j)).le.EPSDEG.or.&
                   DStates%SubStates(iSSt)%EVal(i).lt.emax(j).and.&
                   DStates%SubStates(iSSt)%EVal(i).gt.emin(j)) then
                  this%States(TNSSt,2)=this%States(TNSSt,2)+1
                  this%States(TNSSt,1)=iSSt
                  this%sst_to_statesindex(iSSt) = TNSSt
               endif
            enddo
         enddo
      enddo
      this%NTruncStates=sum(this%States(:,2))
      this%NTruncStatesMax=maxval(this%States(:,2))

      deallocate(energies)
      deallocate(emin)
      deallocate(emax)
      deallocate(degen)
   else
      allocate(this%States(DStates%NSStates, 2))
      this%Egs = DStates%SubStates(0)%EVal(0)
      do iSSt = 0, DStates%NSStates - 1
         this%States(iSSt + 1, 1) = iSSt
         this%States(iSSt + 1, 2) = DStates%SubStates(iSSt)%NStates
         this%sst_to_statesindex(iSSt) = iSSt + 1
         if (DStates%SubStates(iSSt)%EVal(0) < this%Egs)&
            this%Egs = DStates%SubStates(iSSt)%EVal(0)
      end do
      this%NTruncStates = DStates%NStates
      this%NTruncStatesMax = DStates%NStatesMax
   end if


   allocate(this%gu(0:DStates%NBands*2))
   this%gu=0
   allocate(this%NOSOper(DStates%NBands,2))
   this%NOSOper(:,:)=0
   allocate(this%NOSCAOper(DStates%NBands,2,2))
   this%NOSCAOper(:,:,:)=0
   this%NOper=0
   this%lhisto = 0
   this%rhisto = 0
   allocate(this%parttrace(this%NTruncStatesMax))
   allocate(this%tparttrace(this%NTruncStatesMax))

   allocate(this%PreWinFCount(DStates%NBands, 2, 2))
   allocate(this%fprewin(DStates%NBands, 2))
   this%PreWinFCount(:,:,:) = 0
   this%PreWinFCountValid = .true.
   call init_window(this, get_Real_Parameter("TaudiffMax"))
   this%NPairs = 0

! calling get_trace to set all values to the right values of an empty trace
   this%iOperPool=0
   !!! ??? call get_trace_seg the first time ???
   if(.not.get_Integer_Parameter("segment")==1)then
      if (this%b_statesampling) then
         this%outer_sst_size = 1
         ! choose initial outer state with probability proportional to the empty trace value
         allocate(empty_trace_values(this%NTruncStates))
         ! force full trace calculation
         this%outer_sst_old = -1
         this%outer_state_old = -1
         iSt = 1
         do i = 1, size(this%States(:, 1))
            do j = 1, this%States(i, 2)
               this%outer_sst = this%States(i, 1)
               this%outer_state = j
               empty_trace_values(iSt) = trval(get_Trace_EB(this, DStates, global=.true.))
               iSt = iSt + 1
            end do
         enddo
         empty_trace_values = empty_trace_values / sum(empty_trace_values)

         rand_outer = grnd()
         accumulator = 0._KINDR
         iSt = 1
         sstloop: do i = 1, size(this%States(:, 1))
            do j = 1, this%States(i, 2)
               this%outer_sst = this%States(i, 1)
               this%outer_state = j
               accumulator = accumulator + empty_trace_values(iSt)
               if (accumulator >= rand_outer .and. (empty_trace_values(iSt) /= 0._KINDR))&
                  exit sstloop
               iSt = iSt + 1
            end do
         enddo sstloop
         deallocate(empty_trace_values)
      else
         this%outer_state = 1
         this%outer_state_old = 1
         ! choose initial outer superstate with probability proportional to the empty trace value
         allocate(empty_trace_values(size(this%States(:, 1))))
         ! force full trace calculation
         this%outer_sst_old = -1
         do i = 1, size(this%States(:, 1))
            this%outer_sst = this%States(i, 1)
            empty_trace_values(i) = trval(get_Trace_EB(this, DStates, global=.true.))
         enddo
         empty_trace_values = empty_trace_values / sum(empty_trace_values)

         rand_outer = grnd()
         accumulator = 0._KINDR
         do i = 1, size(this%States(:, 1))
            this%outer_sst = this%States(i, 1)
            accumulator = accumulator + empty_trace_values(i)
            if (accumulator >= rand_outer .and. (empty_trace_values(i) /= 0._KINDR)) exit
         enddo
         deallocate(empty_trace_values)
      end if

      this%Trace=get_Trace_EB(this,DStates, global=.true.)
      if (this%Trace%log < -huge(0.0_KINDR)) stop "Zero trace initial state"
      call update_trace_EB(this)
   else
      this%Trace = TLogTr(log = 0.0_KINDR, sign = 1.0_KINDR)
   endif

   allocate(this%tau_c(DStates%NBands,2))
   allocate(this%tau_a(DStates%NBands,2))
!   allocate(this%urho(DStates%NBands,2))
   allocate(this%MInv(DStates%NBands,2))
   if(associated(this%MInv))then
      do iB=1,size(this%MInv(:,1))
      do iS=1,2
         allocate(this%MInv(iB,iS)%Mat(0,0))
      enddo
      enddo
   endif

   allocate(this%MInv_full)
   allocate(this%MInv_full%Mat(0,0))

   this%HybrF=FTau      
   this%HybrF_full=FTau_full
   this%screening_function = screening_function
   this%ParaMag=get_Integer_Parameter("ParaMag")

!    if (get_Integer_Parameter("Phonon")==1 .or. get_Integer_Parameter("Screening")==1) then
   if (get_Integer_Parameter("Uw")==1) then
      this%Phonon = 1
   else
      this%Phonon = 0
   endif

end subroutine init_Trace


!===============================================================================
subroutine dest_Trace(this)
!===============================================================================
   type(TTrace)                        :: this
!local
   type(TOper),pointer                 :: Element,element_temp
   integer                             :: iB,iS,i
   Element=>this%first
   if(associated(Element))then
      do while(associated(Element))
         if (allocated(Element%state)) deallocate(Element%state)
         if (allocated(Element%cache)) deallocate(Element%cache)
         if (allocated(Element%slogmax)) deallocate(Element%slogmax)
         if (allocated(Element%clogmax)) deallocate(Element%clogmax)
!          deallocate(Element%norml)
!          deallocate(Element%normr)
         element_temp=>Element
         Element=>Element%next
         deallocate(element_temp)
      enddo
   endif
   this%first=>null()
   this%last=>null()
   deallocate(this%tau_c)
   deallocate(this%tau_a)

   if(associated(this%MInv))then
      do iB=1,size(this%MInv(:,1))
         do iS=1,2
            if (associated(this%MInv(iB,iS)%Mat)) deallocate(this%MInv(iB,iS)%Mat)
         enddo
      enddo
      deallocate(this%Minv)
   endif
   if(associated(this%MInv_full))then
      if (associated(this%MInv_full%Mat)) deallocate(this%MInv_full%Mat)
      deallocate(this%Minv_full)
   endif

   if(associated(this%gu))deallocate(this%gu)
   if(associated(this%NOSOper))deallocate(this%NOSOper)
   if(associated(this%NOSCAOper))deallocate(this%NOSCAOper)
   if(associated(this%HybrF))deallocate(this%HybrF)
   if(associated(this%HybrF_full))deallocate(this%HybrF_full)
   if(associated(this%screening_function))deallocate(this%screening_function)
   if(allocated(this%parttrace))deallocate(this%parttrace)
!    write(unit=*,fmt=*) "hier!!!!"

   !write(*,*) "DESTROY IT !!!"
   if (associated(this%ket_b2)) deallocate(this%ket_b2)
   if (associated(this%bra_b2)) deallocate(this%bra_b2)
   if (associated(this%ket_b2_logmax)) deallocate(this%ket_b2_logmax)
   if (associated(this%bra_b2_logmax)) deallocate(this%bra_b2_logmax)

   this%iOperPool=0
   i=1
   do while(associated(this%OperPool(i)%p))
     if (allocated(this%OperPool(i)%p%state)) deallocate(this%OperPool(i)%p%state)
     if (allocated(this%OperPool(i)%p%cache)) deallocate(this%OperPool(i)%p%cache)
     if (allocated(this%OperPool(i)%p%slogmax)) deallocate(this%OperPool(i)%p%slogmax)
     if (allocated(this%OperPool(i)%p%clogmax)) deallocate(this%OperPool(i)%p%clogmax)
!    deallocate(this%OperPool(i)%p%norml)
!    deallocate(this%OperPool(i)%p%normr)
     deallocate(this%OperPool(i)%p) 
     i=i+1
   enddo
   
   if (allocated(this%SymMoves)) deallocate(this%SymMoves)
   !!! deallocate segement stuff
   if (allocated(this%outerstate_tmp)) deallocate(this%outerstate_tmp)
   if (allocated(this%initial_outerstate)) deallocate(this%initial_outerstate)
   if (allocated(this%initial_outerstate_old)) deallocate(this%initial_outerstate_old)
   if (allocated(this%mu_accumulator)) deallocate(this%mu_accumulator)
   if (allocated(this%mu_accumulator_tmp)) deallocate(this%mu_accumulator_tmp)
   if (allocated(this%int_accumulator_tmp)) deallocate(this%int_accumulator_tmp)
   if (allocated(this%int_accumulator)) deallocate(this%int_accumulator)
   if (allocated(this%u_matrix)) deallocate(this%u_matrix)
   if (allocated(this%u_matrix_dd)) deallocate(this%u_matrix_dd)
   if (allocated(this%muimp)) deallocate(this%muimp)

   deallocate(this%ffirst)
   deallocate(this%flast)
   deallocate(this%fprewin)

   deallocate(this%PreWinFCount)
   if (allocated(this%States)) deallocate(this%States)   
   if (allocated(this%sst_to_statesindex)) deallocate(this%sst_to_statesindex)
   if (allocated(this%tparttrace)) deallocate(this%tparttrace)
!   if (allocated(this%urho)) deallocate(this%urho)
   !if (allocated(this%omega0)) deallocate(this%omega0)
   !if (allocated(this%g_phonon)) deallocate(this%g_phonon) 
end subroutine dest_Trace

!==============================================================================
subroutine minv_matching_taus(this,MInv)
   !FIXME minv_full
    !TODO update that every time we insert something in the trace.
    type(TSubMatrix),pointer            :: MInv(:,:)
    type(TTrace), intent(inout) :: this
    logical, parameter :: debug = .false.
    character(1), parameter :: CHAR_CA(2) = (/'C','A'/), CHAR_UD(2) = (/'^','v'/)

    integer, dimension(size(MInv,1),2) :: where_c, where_a
    type(TOper), pointer :: op
    integer :: iorb, ispin, order

    if(debug) write(0,'(A)',advance='no') '[minv_matching_taus] sizes: '
    do iorb = 1, size(MInv,1)
        do ispin = 1, 2
            if(allocated(this%tau_c(iorb,ispin)%v)) then
                deallocate(this%tau_c(iorb,ispin)%v)
                deallocate(this%tau_a(iorb,ispin)%v)
            endif
            order = size(MInv(iorb,ispin)%mat,1)  !should be square since pairs
            allocate(this%tau_c(iorb,ispin)%v(order))
            allocate(this%tau_a(iorb,ispin)%v(order))
            if(debug) write(0,"(I3)",advance='no') order
        enddo
        if(debug) write(0,"(A)",advance='no') ','
    enddo
    if(debug) write(0,*)

    where_c = 1
    where_a = 1
    op => this%first
    if(debug) write(0,'(A)',advance='no') '[minv_matching_taus] opers: '
    do while(associated(op))
        if(debug) write(0,"(A1,I1,A1,'(',F5.1,')',1X)",advance='no') &
                      CHAR_CA(op%CA), op%orbital, CHAR_UD(op%spin), op%tau
        if(op%CA == 1) then  !creation
            this%tau_c(op%orbital, op%spin)%v(where_c(op%orbital,op%spin)) = op%tau
            where_c(op%orbital, op%spin) = where_c(op%orbital, op%spin) + 1
        else !annihilation
            this%tau_a(op%orbital, op%spin)%v(where_a(op%orbital,op%spin)) = op%tau
            where_a(op%orbital, op%spin) = where_a(op%orbital, op%spin) + 1
        endif
        op => op%next
    enddo
    if(debug) then
        write(0,*)
        do iorb = 1, size(MInv,1)
          do ispin = 1, 2
            write(0,"('C',I1,A1,': ',100F6.1)") iorb, CHAR_UD(ispin), this%tau_c(iorb,ispin)%v
            write(0,"('A',I1,A1,': ',100F6.1)") iorb, CHAR_UD(ispin), this%tau_a(iorb,ispin)%v
          enddo
        enddo
    endif
end subroutine

! subroutine urho_contracted(u_matrix, this, dstates)
!     !< \brief Computes the U times density, required for improved estimators
!     !!
!     !! Computes the one-particle reduced density operator, fully contracted with
!     !! the U-matrix:
!     !! \f[
!     !!        \sum_{kl} c^+_k(\tau) U_{t k l t} c_l(\tau)
!     !! \f]
!     !! for every creation operator in the local trace (type t and at tau coordinate
!     !! \f$ \tau \f$. This is used in the computation of the improved estimators.
!     use mstates
!     use moperator
!     use msparsematrix

!     type(TTrace), intent(inout) :: this
!     type(TStates), intent(in)   :: dstates
!     real(KINDR)                 :: u_matrix(:,:,:,:)

!     logical, parameter :: debug = .false.
!     character(1), parameter :: CHAR_CA(2) = (/'C','A'/), CHAR_UD(2) = (/'^','v'/)

    
!     integer, dimension(size(this%minv,1),2) :: where_op
!     type(TOper), pointer :: op
!     integer :: iorb, ispin, order, imat, korb, kspin, kmat, lorb, lspin, lmat, i
!     integer :: index_sst, sst, norb, lconn
!     real(KINDR) :: urho_value, u_value, trace, norm

!     !TODO: improved estimator for density-density
!     stop "tried to call urho_contracted"

!     norb = size(this%minv,1)

!     !TODO: correct comment
!     !if(.not.allocated(u_matrix)) stop 'u matrix uninitialized'

!     do iorb = 1, norb
!         do ispin = 1, 2
!             if(allocated(this%urho(iorb,ispin)%v)) then
!                 deallocate(this%urho(iorb,ispin)%v)
!             endif
!             order = size(this%minv(iorb,ispin)%mat,1)  !should be square since pairs
!             allocate(this%urho(iorb,ispin)%v(order))
!         enddo
!     enddo

!     !    trace = sum(this%parttrace)
!     trace = trval(this%Trace)
!     where_op = 1
!     op => this%first
!     if(debug) write(0,'(A)') '[urho_contracted] opers: '
!     do while(associated(op))
!       if(op%CA == 1) then  ! only creation operators!
!         ! cache some improtant quantities
!         iorb = op%orbital
!         ispin = op%spin
!         if(debug) then
!             write(0,"(4X,A1,I1,A1,'(',F5.1,') ')",advance='no') &
!                       CHAR_CA(op%CA), iorb, CHAR_UD(ispin), op%tau
! !           write(0,"(100(I3,L1,:,','))",advance='no') &
! !                     (op%preceding_sst, associated(op%bral(i)%vec), i=1,size(op%bral))
!             write(0,"('>')",advance='no')
!         endif

!         ! compute the corresponding value of urho for that operator by doing
!         ! the sum_kl c+_k U_ikli c_l
!         ! TODO: This is potentially slow and could be sped up by caching the
!         !       operator family R_i := \sum_kl c_k^+ U_ikli c_l
!         ! TODO: one could try and at least pre-compute c_l |ketl> and
!         !       (c_k |bral>)^+
!         urho_value = 0
!         imat = 2*(iorb-1)+ispin
!         do lorb = 1, norb
!           do lspin = 1, 2
!             lmat = 2*(lorb-1)+lspin
!             do korb = 1, norb
!               do kspin = 1, 2
!                 kmat = 2*(korb-1)+kspin
!                 ! Get the corresponding U value and check if they mix
!                 ! FIXME: Here we assumed that the U-matrix is only filled for one
!                 !        particular configuration, i.e., ignoring the "reflected"
!                 !        terms that arise when you go from density-density to
!                 !        general interactions. While this is usually a good assumpution
!                 !        and speeds up the creation of the Hamiltonian, care must
!                 !        be taken when using read-in U matrices. Eventually, this
!                 !        should be handled in interaction.f90 for performance and clarity
!                 u_value = (u_matrix(imat, kmat, imat, lmat) + &
!                            u_matrix(imat, kmat, lmat, imat) + &
!                            u_matrix(kmat, imat, imat, lmat) + &
!                            u_matrix(kmat, imat, lmat, imat))
!                 if(u_value == 0) cycle
                
!                 if(debug) write(0,"(' U',4(I1,A1),':')",advance="no") &
!                       iorb,CHAR_UD(ispin),korb,CHAR_UD(kspin),lorb,CHAR_UD(lspin),&
!                       iorb,CHAR_UD(ispin)
                
!                 ! Now iterate over all superstates (=substates) that 
!                 ! contribute to the current trace at this point.
!                 sst = op%preceding_sst
! !               do index_sst = 1,size(op%bral)
!                     ! First, we need to check if the bral and ketl states are
!                     ! associated, i.e., they have a non-zero contribution to 
!                     ! the trace. This needs to be done first, because the sst
!                     ! array might contain bullshit in this case
! !                   if(.not.associated(op%bral(index_sst)%vec)) cycle
! !                   if(.not.associated(op%ketl(index_sst)%vec)) cycle
                    
!                     ! Get c_korb,kspin, which is an annihilator (2), and see 
!                     ! if it connects to the same block as the c_lorb,lspin. 
!                     ! If it does not connect (< 0), try next (must be seperate)
!                     lconn = dstates%SubStates(sst)%Connect(lorb,lspin,2)
!                     if(lconn < 0) cycle
!                     if(lconn /= dstates%SubStates(sst)%Connect(korb,kspin,2))&
!                         cycle
                        
!                     if(debug) write(0,"(I3,':',I3)",advance="no") index_sst,sst
                    
!                     ! Since the states in the trace are not normed because of
!                     ! the non-Unitary "time" evolution, we need to divide the
!                     ! expectation value by that and also take into account the
!                     ! weight of the current state within the trace, i.e., how
!                     ! strongly the state contributes to this configuration.
! !                   norm = this%parttrace(index_sst) / trace / &
! !                           dot_product(op%bral(index_sst)%vec(:), op%ketl(index_sst)%vec(:))
!                     if(norm == 0) cycle  ! TODO: do some fuzzier check
                    
!                     if(debug) then
!                         write(0,"('(',F5.3,')')",advance="no") abs(norm)
!                         if(.not.associated(dstates%SubStates(sst)%Psis)) then
!                             write (0,*) 'Substate psi is null'
!                         elseif(.not.allocated(dstates%SubStates(sst)%Psis(korb,kspin,2)%dat)) then
!                             write (0,"(/,'c',I1,A1,' invalid for sst')") korb, CHAR_CA(kspin)
!                         elseif(.not.allocated(dstates%SubStates(sst)%Psis(lorb,lspin,2)%dat)) then
!                             write (0,"(/,'c',I1,A1,' invalid for sst')") lorb, CHAR_CA(lspin)
!                         endif
!                     endif
! !                   urho_value = urho_value + u_value*norm*dot_product(&
! !                           matvecprod(dstates%SubStates(sst)%Psis(korb,kspin,2),&
! !                                      op%bral(index_sst)%vec),&
! !                           matvecprod(dstates%SubStates(sst)%Psis(lorb,lspin,2),&
! !                                      op%ketl(index_sst)%vec))
!                     if(debug) write(0,"('=',F5.3,')')",advance="no") urho_value
! !               enddo   ! I
!                 if(debug) write(0,*)
!               enddo   ! HATE
!             enddo   ! THIS
!           enddo   ! SO
!         enddo   ! VERY

!         ! add it to the array
!         this%urho(iorb, ispin)%v(where_op(iorb,ispin)) = urho_value
!         where_op(iorb, ispin) = where_op(iorb, ispin) + 1
!       endif
!       op => op%next
!     enddo  ! MUCH
!     if(debug) write(0,*) 'Done.'
! end subroutine

!===============================================================================
!> Calculate the value of the hybridization function for a certain delta tau
!! using a linear interpolation.
real(KINDR) function get_HybrF(this,tau,Orb,Spin)
!===============================================================================
   type(TTrace)                        :: this
!input
   real(KINDR)                         :: tau
   integer                             :: Orb,Spin
!local
   real(KINDR)                         :: sig
   real(KINDR)                         :: n

   if(tau.lt.0d0)then
      tau=tau+this%beta
      sig=-1d0 
   else
      sig=1d0
   endif
   n=tau/this%beta*dble(size(this%HybrF(Orb,Spin,:))-1)
   get_HybrF=sig*(this%HybrF(Orb,Spin,int(n)+1)+(n-dble(int(n)))*&
      (this%HybrF(Orb,Spin,int(n)+2)-this%HybrF(Orb,Spin,int(n)+1)))
end function get_HybrF

!===============================================================================
!> Calculate the value of the full hybridization function for a certain delta tau
!! using a linear interpolation.
real(KINDR) function get_HybrF_full(this,tau,b1,s1,b2,s2)
!===============================================================================
   type(TTrace)                        :: this
!input
   real(KINDR)                         :: tau
   integer                             :: b1,s1,b2,s2
!local
   real(KINDR)                         :: sig
   real(KINDR)                         :: n

   if(tau.lt.0d0)then
      tau=tau+this%beta
      sig=-1d0 
   else
      sig=1d0
   endif
   n=tau/this%beta*dble(size(this%HybrF_full(b1,s1,b2,s2,:))-1)
   get_HybrF_full=sig*(this%HybrF_full(b1,s1,b2,s2,int(n)+1)+(n-dble(int(n)))*&
      (this%HybrF_full(b1,s1,b2,s2,int(n)+2)-this%HybrF_full(b1,s1,b2,s2,int(n)+1)))
end function get_HybrF_full

!===============================================================================
!> Calculate the value of the bosonic screening function for a certain delta tau
!! using a linear interpolation.
real(KINDR) function get_lin_screening_function(this,tau,Orb1,Spin1,Orb2,Spin2)
!===============================================================================
   type(TTrace)                        :: this
!input
   real(KINDR)                         :: tau
   integer                             :: Orb1,Spin1,Orb2,Spin2
!local
   real(KINDR)                         :: n

   if(tau.lt.0d0) tau = tau + this%beta
   n=tau/this%beta*dble(size(this%screening_function(Orb1,Spin1,Orb2,Spin2,:))-1)
   get_lin_screening_function=&
      (this%screening_function(Orb1,Spin1,Orb2,Spin2,int(n)+1)+&
       (n-dble(int(n)))*&
       (this%screening_function(Orb1,Spin1,Orb2,Spin2,int(n)+2)-&
        this%screening_function(Orb1,Spin1,Orb2,Spin2,int(n)+1)))

end function get_lin_screening_function

!===============================================================================
!> Generate a random permutation of an array.
subroutine rand_perm(arr,n)
!===============================================================================
!local
   integer :: i,j,n,t
   integer :: arr(n)

   !!! algorithm "Random permutation of k out of n objects", Green 1963
   do j=1,n
      i=randint(j, n)

      !!! exchange
      t=arr(j)
      arr(j)=arr(i)
      arr(i)=t
   enddo

end subroutine rand_perm

!===============================================================================
!> Generate a random global update
subroutine gen_GlobalUpdate_new(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   integer                             :: i, gu(DStates%NBands*2)

   ! initialize gu array to set of flavor numbers
   do i=1,DStates%NBands*2
      gu(i)=i
   enddo

   ! permute the array randomly
   call rand_perm(gu,DStates%NBands*2)
   this%gu(1:)=gu

   ! choose randomly whether to flip creators and annihilators
   this%gu(0) = randint(0, 1)

end subroutine gen_GlobalUpdate_new

! !===============================================================================
! !> Generate a random global update. (obsolete, apparently permutations not uniform)
! logical function gen_GlobalUpdate(this,DStates)
! !===============================================================================
!    type(TTrace)                        :: this
! !input
!    type(TStates)                       :: DStates
! !local
!    integer                             :: i,gucheck(-DStates%NBands*2:DStates%NBands*2)
 
!    do i=-DStates%NBands*2,DStates%NBands*2
!       gucheck(i)=i
!    enddo
!    this%gu=gucheck
! !   do while(all(this%gu.eq.gucheck))
!       do i=1,DStates%NBands*2-1
!          this%gu(i:DStates%NBands*2)=cshift(this%gu(i:DStates%NBands*2),int(grnd()*dble(DStates%NBands-i)))
!       enddo
!       this%gu(0)=nint(grnd())
!       gen_GlobalUpdate=.not.all(this%gu.eq.gucheck)
! !   enddo
! end function gen_GlobalUpdate

!===============================================================================
!> Generate a spin flip update.
subroutine gen_SpinFlipUpdate(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   integer                             :: i
 
   this%gu=0

   do i=1,DStates%NBands
      this%gu(i)=i+DStates%NBands
   enddo
   do i=1+DStates%NBands,2*DStates%NBands
      this%gu(i)=i-DStates%NBands
   enddo
end subroutine gen_SpinFlipUpdate

!===============================================================================
!> Choose a random symmetry update specified in the paramter file.
subroutine gen_SymUpdate(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   integer                             :: isym
 
   this%gu=0

   isym = randint(1, 2 * this%NSymMove)
   this%gu(1:DStates%NBands*2)=this%SymMoves(isym,:)
end subroutine gen_SymUpdate

!===============================================================================
!> Prepare global update array for creator/annihilator flip.
subroutine gen_CAFlip(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   integer                             :: i

   this%gu(0) = 1
   do i = 1, 2 * DStates%NBands
      this%gu(i) = i
   end do
end subroutine gen_CAFlip

!===============================================================================
!> Generate an inverse update to reconstruct the original state. Used when the
!! move is not accepted.
subroutine gen_InverseGlobalUpdate(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   integer                             :: i,tgu(0:DStates%NBands*2)

   tgu=this%gu
   do i=1,DStates%NBands*2
     tgu(this%gu(i))=i
   enddo
   this%gu=tgu
end subroutine gen_InverseGlobalUpdate

!===============================================================================
!> Get superstates currently in trace right after begin and right before end.
subroutine get_current_ssts(this, begin, end, begin_sst, end_sst)
!===============================================================================
   type(TTrace)         :: this
!input
   type(TOperPointer)   :: begin, end
!output
   integer, intent(out) :: begin_sst, end_sst

   if (.not. associated(begin%p)) then
      begin_sst = this%outer_sst
   else
      begin_sst = begin%p%following_sst
   endif

   if (.not. associated(end%p)) then
      end_sst = this%outer_sst
   else
      end_sst = end%p%preceding_sst
   endif
end subroutine get_current_ssts

!===============================================================================
!> Check if a valid superstate sequence between the given operators can be found.
logical function check_sst_sequence(this, DStates, begin, end, begin_sst, end_sst)
!===============================================================================
   type(TTrace)                           :: this
!input
   type(TStates)                          :: DStates
   type(TOperPointer)                     :: begin, end
   integer                                :: begin_sst, end_sst
!local
   integer                                :: current_sst
   type(TOper), pointer                   :: Element

   ! begin should point to the last (max tau) operator before the
   ! section to be checked, end to the first (min tau) after
   ! it. begin_sst is the superstate going out of begin and end_sst
   ! the one coming into end. A null pointer represents tau = 0+ as
   ! begin and tau = beta- as end. If end is at lower tau than begin,
   ! the outer superstate is changed if necessary.
   current_sst = begin_sst
   if (associated(begin%p)) then
      Element => begin%p%next
   else
      Element => this%first
   endif

   do
      if (current_sst == -1) then
         check_sst_sequence = .false.
         return
      else if (associated(Element, end%p)) then
         if (current_sst == end_sst) then
            check_sst_sequence = .true.
         else
            check_sst_sequence = .false.
         endif
         return
      else if (.not. associated(Element)) then
         if (this%sst_to_statesindex(current_sst) == -1) then
            check_sst_sequence = .false.
            return
         end if

         this%outer_sst = current_sst

         if (.not. associated(end%p)) then
            if (current_sst == end_sst) then
               check_sst_sequence = .true.
            else
               check_sst_sequence = .false.
            endif
            return
         else
            Element => this%first
         endif
      else
         current_sst = DStates%SubStates(current_sst)&
                              %Connect(Element%Orbital, Element%Spin, Element%CA)
         Element => Element%next
      endif
   enddo
end function check_sst_sequence

!===============================================================================
!> Perform spin-flip update and return probability factor (0 if qn violation, else 1).
real(KINDR) function perform_spinflip_check_qn(this, DStates)
!===============================================================================
   type(TTrace)        :: this
!input
   type(TStates)       :: DStates
!local
   integer             :: old_state, spin_up, spin_down, new_state
   type (TOperPointer) :: nullptr

   old_state = DStates%SubStates(this%outer_sst)%States(0)
   spin_down = ibits(old_state, 0, DStates%NBands)
   spin_up = ibits(old_state, DStates%NBands, DStates%NBands)
   new_state = ior(spin_up, ishft(spin_down, DStates%NBands))

   this%outer_sst = DStates%StatesSStates(new_state, 2)
   call globalUpdate(this, DStates)
   if (this%sst_to_statesindex(this%outer_sst) == -1) then
      perform_spinflip_check_qn = 0._KINDR
   else if (check_sst_sequence(this, DStates, nullptr, nullptr, this%outer_sst, this%outer_sst)) then
      if (this%b_statesampling) then
         perform_spinflip_check_qn = WeightedOuterStateChoice(this, DStates,&
                                                   this%beta - this%last%tau + this%first%tau,&
                                                   this%beta - this%last%tau + this%first%tau,&
                                                   this%outer_sst_old,&
                                                   this%outer_state,&
                                                   (/ this%outer_sst_old /),&
                                                   (/ this%outer_sst /))
      else
         perform_spinflip_check_qn = 1._KINDR
      end if
   else
      perform_spinflip_check_qn = 0._KINDR
   endif
end function perform_spinflip_check_qn

!===============================================================================
!> Select a random new outer superstate of those compatible with all
!  quantum numbers after performing the global update and return the
!  acceptance probability factor.
real(KINDR) function perform_global_check_qn(this, DStates)
!===============================================================================
   type(TTrace)             :: this
!input
   type(TStates)            :: DStates
!local
   integer                  :: possible_ssts_old, possible_ssts
   integer                  :: ssts_old(1:size(this%States(:, 1)))
   integer                  :: ssts_new(1:size(this%States(:, 1))), i, k
   type(TOper),pointer      :: Element

   ! quantum number checking for current configuration (needed for
   ! acceptance probability)
   possible_ssts_old = check_qn_all_ssts(this, DStates, ssts_old)
   if (possible_ssts_old <= 0) then
      stop "unreachable: impossible old configuration in perform_global_check_qn"
   endif

   call globalUpdate(this, DStates)

   ! quantum number checking for potential new configurations
   possible_ssts = check_qn_all_ssts(this, DStates, ssts_new)
   if (possible_ssts <= 0) then
      perform_global_check_qn = 0._KINDR
      return
   endif

   if (this%b_statesampling) then
      perform_global_check_qn = WeightedOuterStateChoice(this, DStates,&
                           this%beta - this%last%tau + this%first%tau,&
                           this%beta - this%last%tau + this%first%tau,&
                           this%outer_sst,&
                           this%outer_state,&
                           ssts_old,&
                           ssts_new)
   else
      perform_global_check_qn = real(possible_ssts, KINDR)/possible_ssts_old

      k = randint(1, possible_ssts)
      do i = 1, size(ssts_new)
         if (ssts_new(i) /= -1) then
            k = k - 1
            if (k <= 0) then
               this%outer_sst = ssts_new(i)
               exit
            endif
         endif
      enddo
   end if
end function perform_global_check_qn

!===============================================================================
!> Perform quantum number checking for all superstates from 0 to beta.
!
!  Returns the number of superstates for which the check passes and
!  writes -1 into array elements of ssts corresponding to superstates
!  for which it does not pass.
function check_qn_all_ssts(this, DStates, ssts) result(possible_ssts)
!===============================================================================
   type(TTrace), intent(in)  :: this
   type(TStates), intent(in) :: DStates
   integer, intent(out)      :: ssts(1:size(this%States(:, 1)))
   integer                   :: possible_ssts, i, k
   type(TOper), pointer      :: Element

   ssts = this%States(:, 1)
   possible_ssts = size(ssts)
   Element => this%first
   do i = 1, this%NOper
      do k = 1, size(ssts)
         if (ssts(k) == -1) cycle
         ssts(k) = DStates%SubStates(ssts(k))&
                          %Connect(Element%Orbital, Element%Spin, Element%CA)
         if (ssts(k) == -1) possible_ssts = possible_ssts - 1
      end do
      if (possible_ssts <= 0) then
         possible_ssts = 0
         return
      end if
      Element => Element%next
   end do
   do k = 1, size(ssts)
      if (ssts(k) == -1) cycle
      if (ssts(k) /= this%States(k, 1)) then
         ssts(k) = -1
         possible_ssts = possible_ssts - 1
      end if
   end do
end function check_qn_all_ssts

!===============================================================================
!> Save current outer superstate (to be used every time before it might be restored)
subroutine save_outer_sst(this)
!===============================================================================
   type(TTrace) :: this

   this%outer_sst_old = this%outer_sst
   this%outer_state_old = this%outer_state
end subroutine save_outer_sst

!===============================================================================
!> Restore old outer superstate (to be used when a move is rejected)
subroutine restore_outer_sst(this)
!===============================================================================
   type(TTrace) :: this

   this%outer_sst = this%outer_sst_old
   this%outer_state = this%outer_state_old
end subroutine restore_outer_sst

!===============================================================================
!> Choose a new outer superstate and state and return the
!  proposal-probability-induced factor for the Metropolis acceptance
!  probability.
!
!  Prerequisites:
!  - save_outer_sst should be called before calling this subroutine
!
!  Parameters:
!  - outer_tau_old: sum of the lengths between 0/beta and the
!    first/last operator in the current ("old") configuration
!  - outer_tau_new: sum of the lengths between 0/beta and the
!    first/last operator in the ("new") configuration to be proposed
!  - outer_sst_old: outer superstate in the old configuration
!  - outer_state_old: outer state in the old configuration
!  - old_weight_ssts: superstates which may be selected as outer
!    superstates when moving from the new to the old configuration
!  - new_weight_ssts: superstates which may be selected as outer
!    superstates when moving from the old to the new configuration
!
!  old_weight_ssts and new_weight_ssts may contain negative entries
!  which are ignored by this function.
real(KINDR) function WeightedOuterStateChoice(this, DStates,&
                                              outer_tau_old,&
                                              outer_tau_new,&
                                              outer_sst_old,&
                                              outer_state_old,&
                                              old_weight_ssts,&
                                              new_weight_ssts)
!===============================================================================
!input
   type(TTrace)             :: this
   type(TStates)            :: DStates
   real(KINDR), intent(in)  :: outer_tau_old, outer_tau_new
   integer, intent(in)      :: outer_sst_old, outer_state_old
   integer, intent(in)      :: old_weight_ssts(:), new_weight_ssts(:)
!local
   integer                  :: i, j, total_size, offset, sstindex
   real(KINDR)              :: total_weight, accumulator, rand
   real(KINDR), allocatable :: new_weights(:)

   ! get proposal probability for the current configuration
   total_weight = 0
   do i = 1, size(old_weight_ssts)
      if (old_weight_ssts(i) >= 0) then
         sstindex = this%sst_to_statesindex(old_weight_ssts(i))
         if (sstindex /= -1) then
            total_weight = total_weight&
            + sum(exp(- outer_tau_old&
                      * (DStates%SubStates(old_weight_ssts(i))&
                                %Eval(0:this%States(sstindex, 2) - 1)&
                         - this%Egs)))
         end if
      end if
   end do

   WeightedOuterStateChoice = exp(- outer_tau_old&
                                    * (DStates%SubStates(outer_sst_old)&
                                              %Eval(outer_state_old - 1)&
                                       - this%Egs))&
                              / total_weight

   ! get number of potential new outer states
   total_size = 0
   do i = 1, size(new_weight_ssts)
      if (new_weight_ssts(i) >= 0) then
         sstindex = this%sst_to_statesindex(new_weight_ssts(i))
         if (sstindex /= -1) then
            total_size = total_size + this%States(sstindex, 2)
         end if
      end if
   end do

   if (total_size > 0) then
      ! get proposal probabilities for all potential new outer states
      allocate(new_weights(total_size))
      offset = 1
      do i = 1, size(new_weight_ssts)
         if (new_weight_ssts(i) >= 0) then
            sstindex = this%sst_to_statesindex(new_weight_ssts(i))
            if (sstindex /= -1) then
               new_weights(offset:offset + this%States(sstindex, 2) - 1)&
                  = exp(- outer_tau_new&
                        * (DStates%SubStates(new_weight_ssts(i))&
                                  %Eval(0:this%States(sstindex, 2) - 1)&
                           - this%Egs))
               offset = offset + this%States(sstindex, 2)
            end if
         end if
      end do
      total_weight = sum(new_weights)

      if (total_weight <= 0_KINDR) then
         this%outer_sst = outer_sst_old
         if (this%b_statesampling) then
            this%outer_state = outer_state_old
         else
            this%outer_state = 1
         end if
         WeightedOuterStateChoice = 0_KINDR
         return
      end if

      rand = grnd() * total_weight

      ! find state and superstate corr. to random number in [0:total_weight]
      accumulator = 0_KINDR
      offset = 0
      new_ssts: do i = 1, size(new_weight_ssts)
         if (new_weight_ssts(i) >= 0) then
            sstindex = this%sst_to_statesindex(new_weight_ssts(i))
            if (sstindex /= -1) then
               do j = 1, this%States(sstindex, 2)
                  accumulator = accumulator + new_weights(offset + j)
                  if (rand <= accumulator .and. (new_weights(offset + j) /= 0_KINDR)) then
                     this%outer_sst = new_weight_ssts(i)
                     if (this%b_statesampling) then
                        this%outer_state = j
                     else
                        this%outer_state = 1
                     end if
                     exit new_ssts
                  end if
               end do
               offset = offset + this%States(sstindex, 2)
            end if
         end if
      end do new_ssts

      if (offset == total_size) then
         write (*, *) "warning: failed to choose new outer state in spite of non-zero weights"
         this%outer_sst = outer_sst_old
         if (this%b_statesampling) then
            this%outer_state = outer_state_old
         else
            this%outer_state = 1
         end if
         WeightedOuterStateChoice = 0_KINDR
         return
      end if

      ! Metropolis: p_acc(i -> j) = w(j)/w(i) * p_prop(j -> i)/p_prop(i -> j)
      ! WeightedOuterStateChoice currently holds p_prop(j -> i)
      WeightedOuterStateChoice =  WeightedOuterStateChoice * (total_weight / new_weights(offset + j))
      return
   else
      WeightedOuterStateChoice = 0_KINDR
      return
   end if
end function WeightedOuterStateChoice

!===============================================================================
!> (Re-)initialize the sliding window.
subroutine init_window(this, taudiff_max)
!===============================================================================
   type(TTrace) :: this
   real(KINDR)  :: taudiff_max
   integer      :: i, j

   this%window_position = 0

   if (get_Integer_Parameter("segment") == 0 .and. taudiff_max > 0_KINDR) then
      this%taudiff_max = taudiff_max
      this%num_windows = max(floor(this%beta / this%taudiff_max) - 1, 1)
      this%tauwin_min = (dble(this%window_position)/(dble(this%num_windows + 1))) * this%beta
      this%tauwin_max = (dble(this%window_position + 2)/(dble(this%num_windows + 1))) * this%beta

      write (*, "('t_max = ', F6.1, ', ', I3, ' window(s), NWin = ', F6.1, ', t_win = ', F6.1)")&
         this%taudiff_max,&
         this%num_windows,&
         dble(this%num_windows + 1)/2d0,&
         this%beta / (dble(this%num_windows + 1)/2d0)
   else
      this%taudiff_max = this%beta
      this%num_windows = 1
      this%tauwin_min = 0
      this%tauwin_max = this%beta
   end if

   this%NPairs = -1
   this%PreWinFCount = 0
   this%PreWinFCountValid = .true.

   this%prewin => null()
   this%postwin => null()
   this%last_prewin => null()
   do i = 1, size(this%fprewin, 1)
      do j = 1, size(this%fprewin, 2)
         this%fprewin(i, j)%p => null()
      end do
   end do

   call set_window(this)
end subroutine init_window

!===============================================================================
!> Shift the sliding window to the next position.
subroutine shift_window(this)
!===============================================================================
   type(TTrace)         :: this
   integer              :: i, j

   if (this%num_windows > 1) then
      this%window_position = this%window_position + 1
      if (this%window_position == this%num_windows) then
         this%window_position = 0
         do j = 1, size(this%fprewin, 2)
            do i = 1, size(this%fprewin, 1)
               this%fprewin(i, j)%p => null()
            end do
         end do
         this%prewin => null()
         this%postwin => null()
      end if
      call set_window(this)
      this%NPairs = -1
      if (this%window_position == 0) then
         call update_PreWinFCount(this)
         this%last_prewin => null()
      end if
   else
      this%window_position = 0
   end if
end subroutine shift_window

!===============================================================================
!> Adjust related window state to a newly set position.
!
!  Preconditions:
!  - this%window_position must be set to a valid position
!  - every this%fprewin(i, j)%p must point to null or an operator before the window
!  - this%postwin must not point to an operator after the first one after the window
!
!  Note that this subroutine does not touch other quantities whose
!  validity might change on moving the window such as NPairs or
!  PreWinFCount.
subroutine set_window(this)
!===============================================================================
   type(TTrace)         :: this
   integer              :: i, j
   real(KINDR)          :: maxtau

   if (this%num_windows > 1) then
      this%tauwin_min = (dble(this%window_position)/(dble(this%num_windows + 1))) * this%beta
      this%tauwin_max = (dble(this%window_position + 2)/(dble(this%num_windows + 1))) * this%beta

      maxtau = -1d0
      this%prewin => null()
      do i = 1, size(this%fprewin, 1)
         do j = 1, size(this%fprewin, 2)

            if (.not. associated(this%fprewin(i, j)%p))&
               this%fprewin(i, j)%p => this%ffirst(i, j)%p

            if (associated(this%fprewin(i, j)%p)) then
               if (this%fprewin(i, j)%p%tau < this%tauwin_min) then
                  do
                     if (.not. associated(this%fprewin(i, j)%p%fnext)) exit
                     if (this%fprewin(i, j)%p%fnext%tau >= this%tauwin_min) exit
                     this%fprewin(i, j)%p => this%fprewin(i, j)%p%fnext
                  end do

                  if (this%fprewin(i, j)%p%tau > maxtau) then
                     maxtau = this%fprewin(i, j)%p%tau
                     this%prewin => this%fprewin(i, j)%p
                  end if
               else
                  this%fprewin(i, j)%p => null()
               end if
            end if

         end do
      end do

      if (.not. associated(this%postwin) .or. this%window_position == 0) then
         if (associated(this%prewin)) then
            this%postwin => this%prewin
         else
            this%postwin => this%first
         end if
      end if
      do
         if (.not. associated(this%postwin)) exit
         if (this%postwin%tau > this%tauwin_max) exit
         this%postwin => this%postwin%next
      end do
   else
      this%prewin => null()
      this%last_prewin => null()
      this%postwin => null()
      do i = 1, size(this%fprewin, 1)
         do j = 1, size(this%fprewin, 2)
            this%fprewin(i, j)%p => null()
         end do
      end do
   end if
end subroutine set_window

!===============================================================================
!> Recalculate the stored number of pairs eligible for removal in the window
!  if necessary.
subroutine ensure_valid_NPairs(this)
!===============================================================================
   type(TTrace)         :: this

   if (this%b_offdiag) then
      call ensure_valid_NPairs_offdiag(this)
   else
      call ensure_valid_NPairs_diag(this)
   end if
end subroutine ensure_valid_NPairs

!===============================================================================
!> Recalculate the stored number of pairs eligible for removal in the window
!  if necessary.
subroutine ensure_valid_NPairs_diag(this)
!===============================================================================
   type(TTrace)         :: this
!local
   type(TOper), pointer :: El, mEl

   if (this%NPairs >= 0) return
   if (this%NPairs < -1) stop "NPairs counting logic error"

   this%NPairs = 0
   if (associated(this%prewin)) then
      El => this%prewin%next
   else
      El => this%first
   end if
   do
      if (.not. associated(El)) return
      if (El%tau >= this%tauwin_min) exit

      El => El%next
   end do

   ! count pairs of one creator and one annihilator in the current window
   ! with distance < taudiff_max and same flavours
   do
      if (.not. associated(El)) exit
      if (El%tau > this%tauwin_max) exit
      if (El%has_hyb) then

         mEl => El%fnext
         do
            if (.not. associated(mEl)) exit
            if (mEl%tau > this%tauwin_max) exit
            if (mEl%tau > El%tau + this%taudiff_max) exit

            if (mEl%has_hyb .and. mEl%CA /= El%CA)&
               this%NPairs = this%NPairs + 1

            mEl => mEl%fnext
         end do

      end if

      El => El%next
   end do
end subroutine ensure_valid_NPairs_diag

!===============================================================================
!> Recalculate the stored number of pairs eligible for removal in the window
!  if necessary.
subroutine ensure_valid_NPairs_offdiag(this)
!===============================================================================
   type(TTrace)         :: this
!local
   type(TOper), pointer :: El, mEl

   if (this%NPairs >= 0) return
   if (this%NPairs < -1) stop "NPairs counting logic error"

   this%NPairs = 0
   if (associated(this%prewin)) then
      El => this%prewin%next
   else
      El => this%first
   end if
   do
      if (.not. associated(El)) return
      if (El%tau >= this%tauwin_min) exit

      El => El%next
   end do

   ! count pairs of one creator and one annihilator in the current window
   ! with distance < taudiff_max
   do
      if (.not. associated(El)) exit
      if (El%tau > this%tauwin_max) exit
      if (El%has_hyb) then

         mEl => El%next
         do
            if (.not. associated(mEl)) exit
            if (mEl%tau > this%tauwin_max) exit
            if (mEl%tau > El%tau + this%taudiff_max) exit

            if (mEl%has_hyb .and. mEl%CA /= El%CA)&
               this%NPairs = this%NPairs + 1

            mEl => mEl%next
         end do

      end if

      El => El%next
   end do
end subroutine ensure_valid_NPairs_offdiag

!===============================================================================
!> Update the stored number of pairs eligible for removal in the window.
!  Oper: an array of two different operators in time order
!  insrem: +1 after operator insertion, -1 after operator removal
subroutine update_NPairs(this, Oper, insrem)
!===============================================================================
   type(TTrace)         :: this
   type(TOperPointer)   :: Oper(2)
   integer              :: insrem

   if (this%b_offdiag) then
      call update_NPairs_offdiag(this, Oper, insrem)
   else
      call update_NPairs_diag(this, Oper, insrem)
   end if
end subroutine update_NPairs

!===============================================================================
!> Update the stored number of pairs eligible for removal in the window.
subroutine update_NPairs_diag(this, Oper, insrem)
!===============================================================================
   type(TTrace)         :: this
   type(TOperPointer)   :: Oper(2)
   integer              :: insrem
!local
   integer              :: i
   type(TOper), pointer :: El, mEl

   if (this%NPairs <= -1) then
      call ensure_valid_NPairs_diag(this)
      return
   end if

   ! When removing a pair of non-(f-)adjacent operators, there will be no
   ! pointers to them left in the rest of the chain, so the removed pair must be
   ! subtracted separately
   if (insrem == -1 .and. .not. associated(Oper(1)%p%fnext, Oper(2)%p)) this%NPairs = this%NPairs - 1

   ! Subtract/add pairs which have one of the changed operators as member.
   ! Count the changed pair only once by excluding it on the way down and not
   ! having multiple operators with equal tau.
   do i = 1, 2
      El => Oper(i)%p

      mEl => El%fnext
      do
         if (.not. associated(mEl)) exit
         if (mEl%tau > this%tauwin_max) exit
         if (mEl%tau > El%tau + this%taudiff_max) exit

         if (mEl%has_hyb .and. mEl%CA /= El%CA)&
            this%NPairs = this%NPairs + insrem

         mEl => mEl%fnext
      end do

      mEl => El%fprev
      do
         if (.not. associated(mEl)) exit
         if (mEl%tau < this%tauwin_min) exit
         if (mEl%tau < El%tau - this%taudiff_max) exit

         if (mEl%has_hyb .and. mEl%CA /= El%CA .and. .not. associated(mEl, Oper(3 - i)%p))&
            this%NPairs = this%NPairs + insrem

         mEl => mEl%fprev
      end do
   end do
end subroutine update_NPairs_diag

!===============================================================================
!> Update the stored number of pairs eligible for removal in the window.
!  insrem: +1 after operator insertion, -1 after operator removal
subroutine update_NPairs_offdiag(this, Oper, insrem)
!===============================================================================
   type(TTrace)         :: this
   type(TOperPointer)   :: Oper(2)
   integer              :: insrem
!local
   integer              :: i
   type(TOper), pointer :: El, mEl

   if (this%NPairs <= -1) then
      call ensure_valid_NPairs_offdiag(this)
      return
   end if

   ! When removing a pair of non-adjacent operators, there will be no
   ! pointers to them left in the rest of the chain, so the removed pair must be
   ! subtracted separately
   if (insrem == -1 .and. .not. associated(Oper(1)%p%next, Oper(2)%p)) this%NPairs = this%NPairs - 1

   ! Subtract/add pairs which have one of the changed operators as member.
   ! Count the changed pair only once by excluding it on the way down and not
   ! having multiple operators with equal tau.
   do i = 1, 2
      El => Oper(i)%p

      mEl => El%next
      do
         if (.not. associated(mEl)) exit
         if (mEl%tau > this%tauwin_max) exit
         if (mEl%tau > El%tau + this%taudiff_max) exit

         if (mEl%has_hyb .and. mEl%CA /= El%CA)&
            this%NPairs = this%NPairs + insrem

         mEl => mEl%next
      end do

      mEl => El%prev
      do
         if (.not. associated(mEl)) exit
         if (mEl%tau < this%tauwin_min) exit
         if (mEl%tau < El%tau - this%taudiff_max) exit

         if (mEl%has_hyb .and. mEl%CA /= El%CA .and. .not. associated(mEl, Oper(3 - i)%p))&
            this%NPairs = this%NPairs + insrem

         mEl => mEl%prev
      end do
   end do
end subroutine update_NPairs_offdiag

!===============================================================================
!> Check if a pair of operators is a removable pair.
!  op1, op2: Operators
logical function is_rempair(this, op1, op2)
!===============================================================================
   type(TTrace), intent(in) :: this
   type(TOper), intent(in)  :: op1, op2

   if (this%b_offdiag) then
      is_rempair = is_rempair_offdiag(this, op1, op2)
   else
      is_rempair = is_rempair_diag(this, op1, op2)
   end if
end function is_rempair

logical function is_rempair_diag(this, op1, op2)
   type(TTrace), intent(in) :: this
   type(TOper), intent(in)  :: op1, op2

   if (op1%Orbital == op2%Orbital .and. op1%Spin == op2%Spin&
       .and. is_rempair_offdiag(this, op1, op2)) then
      is_rempair_diag = .true.
   else
      is_rempair_diag = .false.
   end if
end function is_rempair_diag

logical function is_rempair_offdiag(this, op1, op2)
   type(TTrace), intent(in) :: this
   type(TOper), intent(in)  :: op1, op2

   if (op1%CA /= op2%CA .and. op1%has_hyb .and. op2%has_hyb&
       .and. op1%tau >= this%tauwin_min .and. op1%tau <= this%tauwin_max&
       .and. op2%tau >= this%tauwin_min .and. op2%tau <= this%tauwin_max&
       .and. abs(op1%tau - op2%tau) <= this%taudiff_max) then
      is_rempair_offdiag = .true.
   else
      is_rempair_offdiag = .false.
   end if
end function is_rempair_offdiag

!===============================================================================
!> Ensure that the stored number of operators before the window is correct.
subroutine ensure_valid_PreWinFCount(this)
!===============================================================================
   type(TTrace)         :: this
!local
   type(TOper), pointer :: El

   if (this%PreWinFCountValid) then
      if (associated(this%FCountLastIncluded, this%prewin) .or.&
          (.not. associated(this%prewin) .and. .not. associated(this%FCountLastIncluded))) then
         return
      else
         call update_PreWinFCount(this)
      end if
   end if

   ! recount
   this%PreWinFCountValid = .true.
   this%PreWinFCount = 0
   this%FCountLastIncluded => null()
   El => this%first
   if (associated(El)) then
      if (El%tau >= this%tauwin_min) return
   else
      return
   end if

   do
      if (.not. associated(El)) exit
      if (El%has_hyb) this%PreWinFCount(El%Orbital, El%Spin, El%CA)&
                      = this%PreWinFCount(El%Orbital, El%Spin, El%CA) + 1
      if (associated(El, this%prewin)) exit
      El => El%next
   end do

   if (associated(El)) then
      this%FCountLastIncluded => El
   else
      this%FCountLastIncluded => this%last
   end if
end subroutine ensure_valid_PreWinFCount

!===============================================================================
!> Update the stored number of operators before the window after shifting it.
!
!  Preconditions:
!  - PreWinFCount must hold the correct count up to and including FCountLastIncluded
!  OR
!  - PreWinFCountValid == false
!  OR
!  - window_position == 0
subroutine update_PreWinFCount(this)
!===============================================================================
   type(TTrace)         :: this
!local
   type(TOper), pointer :: El

   if (this%window_position == 0) then
      this%PreWinFCount = 0
      this%PreWinFCountValid = .true.
      this%FCountLastIncluded => null()
      return
   end if

   ! if a full recount is necessary, defer it
   if (.not. this%PreWinFCountValid) return

   El => this%FCountLastIncluded
   if (associated(El)) then
      El => El%next
   else
      El => this%first
   end if
   if (.not. associated(El)) return
   if (El%tau >= this%tauwin_min) return

   do
      if (.not. associated(El)) exit
      if (El%has_hyb) this%PreWinFCount(El%Orbital, El%Spin, El%CA)&
                      = this%PreWinFCount(El%Orbital, El%Spin, El%CA) + 1
      if (associated(El, this%prewin)) exit
      El => El%next
   end do

   if (associated(El)) then
      this%FCountLastIncluded => El
   else
      this%FCountLastIncluded => this%last
   end if
end subroutine update_PreWinFCount

!===============================================================================
!> Perform a global update.
subroutine globalUpdate(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   type(TOper),pointer                 :: Element
   integer                             :: i, j
   integer                             :: TNOSOper(size(this%NOSOper(:,1)),2)
   integer                             :: TNOSCAOper(size(this%NOSCAOper(:,1, 1)),2, 2)
   integer                             :: PreWinFCount_temp(DStates%NBands, 2, 2)
   type(TOperPointer)                  :: ffirst_temp(DStates%NBands, 2)
   type(TOperPointer)                  :: flast_temp(DStates%NBands, 2)
   type(TOperPointer)                  :: fprewin_temp(DStates%NBands, 2)

   if(this%NOper.eq.0)return
   !write(*,*) " "
   !write(*,*) "global_update"
   
   Element=>this%first
   do while(associated(Element))
      do i = 1, DStates%NBands*2
         if(mod((i-1),DStates%NBands)+1.eq.Element%Orbital.and.&
            (i-1)/DStates%NBands+1.eq.Element%Spin)then
            Element%Orbital=abs(mod(this%gu(i)-1,DStates%NBands)+1)
            Element%Spin=abs((this%gu(i)-1)/DStates%Nbands+1)
            if(this%gu(0).ne.0)Element%CA=3-Element%CA
            exit
         endif
      enddo
      Element=>Element%next
   enddo

   do i=1,DStates%NBands*2
      TNOSOper(mod(this%gu(i)-1,DStates%NBands)+1,(this%gu(i)-1)/DStates%Nbands+1)=&
         this%NOSOper(mod((i-1),DStates%NBands)+1,(i-1)/DStates%NBands+1)
      ffirst_temp(mod(this%gu(i)-1, DStates%NBands) + 1, (this%gu(i)-1)/DStates%NBands + 1)=&
         this%ffirst(mod((i-1), DStates%NBands) + 1, (i-1)/DStates%NBands + 1)
      flast_temp(mod(this%gu(i)-1, DStates%NBands) + 1,(this%gu(i)-1)/DStates%NBands + 1)=&
         this%flast(mod((i-1), DStates%NBands) + 1, (i-1)/DStates%NBands + 1)
      fprewin_temp(mod(this%gu(i)-1, DStates%NBands) + 1,(this%gu(i)-1)/DStates%NBands + 1)=&
         this%fprewin(mod((i-1), DStates%NBands) + 1, (i-1)/DStates%NBands + 1)
      if (this%gu(0) == 0) then
         PreWinFCount_temp(mod(this%gu(i)-1, DStates%NBands) + 1,(this%gu(i)-1)/DStates%NBands + 1, :)=&
            this%PreWinFCount(mod((i-1), DStates%NBands) + 1, (i-1)/DStates%NBands + 1, :)
         TNOSCAOper(mod(this%gu(i)-1,DStates%NBands)+1,(this%gu(i)-1)/DStates%Nbands+1, :)=&
            this%NOSCAOper(mod((i-1),DStates%NBands)+1,(i-1)/DStates%NBands+1, :)
      else
         do j = 1, 2
            PreWinFCount_temp(mod(this%gu(i)-1, DStates%NBands) + 1,(this%gu(i)-1)/DStates%NBands + 1, j)=&
               this%PreWinFCount(mod((i-1), DStates%NBands) + 1, (i-1)/DStates%NBands + 1, 3 - j)
            TNOSCAOper(mod(this%gu(i)-1,DStates%NBands)+1,(this%gu(i)-1)/DStates%Nbands+1, j)=&
               this%NOSCAOper(mod((i-1),DStates%NBands)+1,(i-1)/DStates%NBands+1, 3 - j)
         end do
      end if
   enddo
   this%NOSOper = TNOSOper
   this%NOSCAOper = TNOSCAOper
   this%ffirst = ffirst_temp
   this%flast = flast_temp
   this%fprewin = fprewin_temp
   this%PreWinFCount = PreWinFCount_temp
end subroutine globalUpdate

!===============================================================================
logical function propose_flavourexchange_general(this,oper)
!===============================================================================
!input
   type(TTrace)                        :: this
!output
   type(TOperPointer)                  :: Oper(2)
!local    
   type(TOper),pointer                 :: Element
   integer                             :: pos1, pos2, N

   propose_flavourexchange_general=.true.

   ! count operators inside window
   N = 0
   if (.not. associated(this%prewin)) then
      Element => this%first
   else
      Element => this%prewin%next
   end if

   do
      if (.not. associated(Element)) exit
      if (associated(Element, this%postwin)) exit
      N = N + 1
      Element => Element%next
   end do

   if (N <= 1) then
      propose_flavourexchange_general = .false.
      return
   end if


   ! choose operator to be changed with lesser tau
   pos1 = randint(1, N - 1)
   ! with greater tau, position as pos1 + pos2
   pos2 = randint(1, N - pos1)


   ! find operators
   N = 0
   if (.not. associated(this%prewin)) then
      Element => this%first
   else
      Element => this%prewin%next
   end if

   do
      if (.not. associated(Element))&
         stop "unreachable: miscounted (1) in propose_flavourexchange_general"
      if (associated(Element, this%postwin))&
         stop "unreachable: miscounted (2) in propose_flavourexchange_general"

      pos1 = pos1 - 1
      if (pos1 <= 0) then
         Oper(1)%p => Element
         exit
      endif

      Element => Element%next
   enddo

   do
      Element => Element%next
      pos2 = pos2 - 1

      if (.not. associated(Element))&
         stop "unreachable: miscounted (3) in propose_flavourexchange_general"
      if (associated(Element, this%postwin))&
         stop "unreachable: miscounted (4) in propose_flavourexchange_general"

      if (pos2 <= 0) then
         Oper(2)%p => Element
         exit
      endif
   enddo

   ! in component sampling worm flavor needs to be kept fixed
   if( oper(1)%p%has_hyb .eqv. .false. .or. oper(2)%p%has_hyb .eqv. .false. ) then
      propose_flavourexchange_general = .false.
   endif  

 
end function propose_flavourexchange_general

!===============================================================================
!> Regenerate all arrays containing the possible hybridizations after a global
!! update and returns the inverted hybrizization matrix.
function get_InvFullHybr(this,DStates,Det)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!output
   type(TSubMatrix),pointer            :: get_InvFullHybr(:,:)
   type(TLogDet),intent(out)           :: Det
!local
   integer                             :: iB,iS,iters,Orb,Spin
   integer                             :: pos(DStates%NBands,2,2)
   type(TOper),pointer                 :: ElementS,ElementE
   type(TLogDet)                       :: tDet

   allocate(get_InvFullHybr(DStates%NBands,2))
   do iB=1,DStates%NBands
      do iS=1,2
         allocate(get_InvFullHybr(iB,iS)%Mat(this%NOSOper(iB,iS)/2,this%NOSOper(iB,iS)/2))
      enddo
   enddo
   ElementS=>this%first
   pos=1
   do iters=1,this%NOper
      if(ElementS%has_hyb) then
         Orb=ElementS%Orbital
         Spin=ElementS%Spin
         if(ElementS%CA.eq.1)then
            ElementE => this%ffirst(Orb, Spin)%p
            pos(Orb,Spin,2)=1
            do while (associated(ElementE))
               if(ElementE%has_hyb) then
                  if(ElementE%CA.eq.2)then
                     get_InvFullHybr(Orb,Spin)%Mat(pos(Orb,Spin,2),pos(Orb,Spin,1))=&
                        get_HybrF(this,ElementE%tau-ElementS%tau,Orb,Spin)
                     pos(Orb,Spin,2)=pos(Orb,Spin,2)+1
                  endif
               endif
               ElementE => ElementE%fnext
            enddo
            pos(Orb,Spin,1)=pos(Orb,Spin,1)+1
         endif
      endif
      ElementS=>ElementS%next
   enddo
   !determinant and inverse
   Det = TLogDet(log = 0.0_KINDR, sign = 1.0_KINDR)
   do iB=1,DStates%NBands
      do iS=1,2
         if(size(get_InvFullHybr(iB,iS)%Mat).ne.0)then
            call get_MatLogDetFull(this%NOSOper(iB,iS)/2,get_InvFullHybr(iB,iS)%Mat,tDet)
            Det = Det * tDet
         endif
      enddo
   enddo
end function get_InvFullHybr

!===============================================================================
!> Update all arrays containing the possible hybridizations and the
!! determinant after a tau-shift move.
subroutine update_InvFullHybrShiftTau(this, DStates, wrapped)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
   integer                             :: wrapped(DStates%NBands, 2, 2)

   if (this%b_offdiag) then
      call update_InvFullHybrShiftTau_offdiag(this, DStates, wrapped)
   else
      call update_InvFullHybrShiftTau_diag(this, DStates, wrapped)
   end if
end subroutine update_InvFullHybrShiftTau

!===============================================================================
!> Update all arrays containing the possible hybridizations and the
!! determinant after a tau-shift move.
subroutine update_InvFullHybrShiftTau_diag(this, DStates, wrapped)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
   integer                             :: wrapped(DStates%NBands, 2, 2)
!local
   integer                             :: iB, iS, i, j
   real(KINDR)                         :: DetRat
   type(TSubMatrix),pointer            :: MInvNew(:,:)

   DetRat = 1.0_KINDR
   allocate(MInvNew(DStates%NBands, 2))

   ! the effect of a tau-shift move on the hybridization matrix can
   ! always be expressed as a cyclic shift and negation of as many of
   ! the first rows and columns of the matrix as annihilation and
   ! creation operators were shifted over 0 and back in from
   ! beta. this row and column shift and negation can alternatively be
   ! expressed through the multiplication of an easily constructible
   ! matrix from the left and right each. on inversion, this product
   ! turns into the product of the inverse matrices in reverse
   ! order. the inverse of the two constructed matrices is easily
   ! found and can again be interpreted as a cyclic row and column
   ! shift and negation, which are exactly the operations performed by
   ! the code below (into a new inverse hybridization matrix, not
   ! in-place) for a fast update. it also tracks the total sign for
   ! the determinant, which depends on the number of swapped
   ! columns/rows (matrix dimension - 1 per cyclically shifted
   ! column/row) and column/row negations (1 per cyclically shifted
   ! column/row)

   ! note that for our hybridization matrix the row indices correspond
   ! to annihilator indices and thus for the inverse that is modified
   ! below to creator indices, and that the tau-shift is performed
   ! such that the operators are shifted toward lesser tau and thus
   ! "wrapping" a row means that the entries with row index 1 are
   ! cyclically shifted such that they are at the positions with the
   ! maximum row index afterward
   do iB = 1, DStates%NBands
      do iS = 1, 2

         ! minv exactly unchanged if none or all operators of one
         ! flavour wrapped around
         if ((wrapped(iB, iS, 1) == 0 .and. wrapped(iB, iS, 2) == 0) .or.&
             (wrapped(iB, iS, 1) == size(this%MInv(iB, iS)%Mat, 1) .and.&
              wrapped(iB, iS, 2) == size(this%MInv(iB, iS)%Mat, 2))) then
            MInvNew(iB, iS)%Mat => this%MInv(iB, iS)%Mat
         else
            allocate(MInvNew(iB, iS)%Mat(this%NOSOper(iB, iS)/2, this%NOSOper(iB, iS)/2))

            ! copy unwrapped columns
            do j = 1, size(this%MInv(iB, iS)%Mat, 2) - wrapped(iB, iS, 2)

               ! copy unwrapped rows
               do i = 1, size(this%MInv(iB, iS)%Mat, 1) - wrapped(iB, iS, 1)
                  MInvNew(iB, iS)%Mat(i, j) =&
                     this%MInv(iB, iS)%Mat(i + wrapped(iB, iS, 1), j + wrapped(iB, iS, 2))
               end do

               ! copy wrapped rows
               do i = 1, wrapped(iB, iS, 1)
                  MInvNew(iB, iS)%Mat(i + size(this%MInv(iB, iS)%Mat, 1) - wrapped(iB, iS, 1),&
                                      j) =&
                     -this%MInv(iB, iS)%Mat(i, j + wrapped(iB, iS, 2))
               end do

            end do

            ! copy wrapped columns
            do j = 1, wrapped(iB, iS, 2)

               do i = 1, size(this%MInv(iB, iS)%Mat, 1) - wrapped(iB, iS, 1)
                  MInvNew(iB, iS)%Mat(i,&
                                      j + size(this%MInv(iB, iS)%Mat, 2) - wrapped(iB, iS, 2)) =&
                     -this%MInv(iB, iS)%Mat(i + wrapped(iB, iS, 1), j)
               end do

               do i = 1, wrapped(iB, iS, 1)
                  MInvNew(iB, iS)%Mat(i + size(this%MInv(iB, iS)%Mat, 1) - wrapped(iB, iS, 1),&
                                      j + size(this%MInv(iB, iS)%Mat, 2) - wrapped(iB, iS, 2)) =&
                     this%MInv(iB, iS)%Mat(i, j)
               end do

            end do

            deallocate(this%MInv(iB, iS)%Mat)
         end if

         ! if the matrix dimension and the number of wrapped rows and
         ! columns are odd, the determinant of the hybridization
         ! matrix changes its sign
         if (modulo(size(MInvNew(iB, iS)%Mat, 1), 2) /= 0 .and.&
             modulo(wrapped(iB, iS, 1) + wrapped(iB, iS, 2), 2) /= 0)&
            DetRat = -DetRat

      end do
   end do

   deallocate(this%MInv)
   this%MInv => MInvNew
   this%Det%sign = this%Det%sign * DetRat
end subroutine update_InvFullHybrShiftTau_diag

!===============================================================================
!> Update all arrays containing the possible hybridizations and the
!! determinant after a tau-shift move.
subroutine update_InvFullHybrShiftTau_offdiag(this, DStates, wrapped)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
   integer                             :: wrapped(DStates%NBands, 2, 2)
!local
   ! Flavour indices of the rows and columns of the current block
   integer                             :: iBr, iSr, iBc, iSc
   ! The number of rows/columns before and row/column size of the current block
   integer                             :: iR, iC, sR, sC
   integer                             :: i, j
   real(KINDR)                         :: DetRat
   real(KINDR),pointer                 :: MInvFullNew(:,:)

   DetRat = 1.0_KINDR
   if (allocated(this%wormContainer)) then
      allocate(MInvFullNew((this%NOper-size(this%wormContainer))/2,&
                           (this%NOper-size(this%wormContainer))/2))
   else
      allocate(MInvFullNew(this%NOper/2, this%NOper/2))
   end if

   ! the effect of a tau-shift move on the hybridization matrix can
   ! always be expressed as a cyclic shift and negation of as many of
   ! the first rows and columns of the matrix as annihilation and
   ! creation operators were shifted over 0 and back in from
   ! beta. this row and column shift and negation can alternatively be
   ! expressed through the multiplication of an easily constructible
   ! matrix from the left and right each. on inversion, this product
   ! turns into the product of the inverse matrices in reverse
   ! order. the inverse of the two constructed matrices is easily
   ! found and can again be interpreted as a cyclic row and column
   ! shift and negation, which are exactly the operations performed by
   ! the code below (into a new inverse hybridization matrix, not
   ! in-place) for a fast update. it also tracks the total sign for
   ! the determinant, which depends on the number of swapped
   ! columns/rows (matrix dimension - 1 per cyclically shifted
   ! column/row) and column/row negations (1 per cyclically shifted
   ! column/row)

   ! note that for our hybridization matrix the row indices correspond
   ! to annihilator indices and thus for the inverse that is modified
   ! below to creator indices, and that the tau-shift is performed
   ! such that the operators are shifted toward lesser tau and thus
   ! "wrapping" a row means that the entries with row index 1 are
   ! cyclically shifted such that they are at the positions with the
   ! maximum row index afterward
   do iBr = 1, DStates%NBands
      do iSr = 1, 2

         ! number of rows before the current block is the number of
         ! creators of earlier flavours, where flavour ordering is
         ! ascending by orbital and then spin index
         iR = 0
         do i = 1, iBr
            do j = 1, 2
               if (i < iBr .or. j < iSr) iR = iR + this%NOSCAOper(i, j, 1)
            end do
         end do
         sR = this%NOSCAOper(iBr, iSr, 1)
         if (sR <= 0) cycle

         do iBc = 1, DStates%NBands
            do iSc = 1, 2

               ! number of columns analogously with annihilators
               iC = 0
               do i = 1, iBc
                  do j = 1, 2
                     if (i < iBc .or. j < iSc) iC = iC + this%NOSCAOper(i, j, 2)
                  end do
               end do
               sC = this%NOSCAOper(iBc, iSc, 2)
               if (sC <= 0) cycle

               if ((wrapped(iBr, iSr, 1) == 0 .and. wrapped(iBc, iSc, 2) == 0)&
                   .or. (wrapped(iBr, iSr, 1) == sR .and. wrapped(iBc, iSc, 2) == sC)) then
                  MInvFullNew(iR+1:iR+sR, iC+1:iC+sC) = this%MInv_full%Mat(iR+1:iR+sR, iC+1:iC+sC)
               else

                  ! copy unwrapped columns
                  do j = 1, sC - wrapped(iBc, iSc, 2)

                     ! copy unwrapped rows
                     do i = 1, sR - wrapped(iBr, iSr, 1)
                        MInvFullNew(iR + i, iC + j) =&
                           this%MInv_full%Mat(iR + i + wrapped(iBr, iSr, 1),&
                                              iC + j + wrapped(iBc, iSc, 2))
                     end do

                     ! copy wrapped rows
                     do i = 1, wrapped(iBr, iSr, 1)
                        MInvFullNew(iR + i + sR - wrapped(iBr, iSr, 1), iC + j) =&
                           -this%MInv_full%Mat(iR + i, iC + j + wrapped(iBc, iSc, 2))
                     end do

                  end do

                  ! copy wrapped columns
                  do j = 1, wrapped(iBc, iSc, 2)

                     do i = 1, sR - wrapped(iBr, iSr, 1)
                        MInvFullNew(iR + i, iC + j + sC - wrapped(iBc, iSc, 2)) =&
                           -this%MInv_full%Mat(iR + i + wrapped(iBr, iSr, 1), iC + j)
                     end do

                     do i = 1, wrapped(iBr, iSr, 1)
                        MInvFullNew(iR + i + sR - wrapped(iBr, iSr, 1),&
                                    iC + j + sC - wrapped(iBc, iSc, 2)) =&
                           this%MInv_full%Mat(iR + i, iC + j)
                     end do

                  end do
               end if
            end do
         end do
      end do
   end do

   ! wrapping rows and columns is done per block; the extra sign is
   ! (-1)**(size*Nwrapped) per block
   do iBr = 1, DStates%NBands
      do iSr = 1, 2
         sR = this%NOSCAOper(iBr, iSr, 1)
         if (modulo(sR * wrapped(iBr, iSr, 1), 2) /= 0) DetRat = -DetRat
      end do
   end do

   do iBc = 1, DStates%NBands
      do iSc = 1, 2
         sC = this%NOSCAOper(iBc, iSc, 2)
         if (modulo(sC * wrapped(iBc, iSc, 2), 2) /= 0) DetRat = -DetRat
      end do
   end do

   deallocate(this%MInv_full%Mat)
   this%MInv_full%Mat => MInvFullNew
   this%Det%sign = this%Det%sign * DetRat
end subroutine update_InvFullHybrShiftTau_offdiag

!===============================================================================
!> Regenerate all arrays containing the possible hybridizations after a global
!! update and returns the hybrizization matrix.
function get_FullHybr(this,DStates,Det)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!output
   type(TSubMatrix),pointer            :: get_FullHybr(:,:)
   type(TLogDet), intent(out)          :: Det
!local
   integer                             :: iB,iS,iters,Orb,Spin
   integer                             :: pos(DStates%NBands,2,2)
   type(TOper),pointer                 :: ElementS,ElementE
   type(TLogDet)                       :: tDet

   allocate(get_FullHybr(DStates%NBands,2))
   do iB=1,DStates%NBands
      do iS=1,2
         allocate(get_FullHybr(iB,iS)%Mat(this%NOSOper(iB,iS)/2,this%NOSOper(iB,iS)/2))
      enddo
   enddo
   ElementS=>this%first
   pos=1
   do iters=1,this%NOper
      if(ElementS%has_hyb) then
         Orb=ElementS%Orbital
         Spin=ElementS%Spin
         if(ElementS%CA.eq.1)then
            ElementE => this%ffirst(Orb, Spin)%p
            pos(Orb,Spin,2)=1
            do while (associated(ElementE))
               if(ElementE%has_hyb) then
                  if(ElementE%CA.eq.2)then
                     get_FullHybr(Orb,Spin)%Mat(pos(Orb,Spin,2),pos(Orb,Spin,1))=&
                        get_HybrF(this,ElementE%tau-ElementS%tau,Orb,Spin)
                     pos(Orb,Spin,2)=pos(Orb,Spin,2)+1
                  endif
               endif
               ElementE => ElementE%fnext
            enddo
            pos(Orb,Spin,1)=pos(Orb,Spin,1)+1
         endif
      endif
      ElementS=>ElementS%next
   enddo
   !determinant
   Det = TLogDet(log = 0.0_KINDR, sign = 1.0_KINDR)
   do iB=1,DStates%NBands
      do iS=1,2
         if(size(get_FullHybr(iB,iS)%Mat).ne.0)then
            call get_LogDetFull(this%NOSOper(iB,iS)/2,get_FullHybr(iB,iS)%Mat,tDet)
            Det=tDet*Det
         endif
      enddo
   enddo
end function get_FullHybr

!===============================================================================
!> Get the determinant of all matrices. Needed for global updates since we can
!! not use inversion by partitioning then.
type(TLogDet) function get_FullLogDet(this,DStates)
!===============================================================================
   type(TTrace)               :: this
!input
   type(TStates)              :: DStates
!local
   integer                    :: iB,iS
   type(TLogDet)              :: tDet

   get_FullLogDet = TLogDet(log = 0.0_KINDR, sign = 1.0_KINDR)
   do iB=1,DStates%NBands
      do iS=1,2
         if(associated(this%MinV(iB,iS)%Mat))then
         if(size(this%Minv(iB,iS)%Mat(1,:)).gt.0)then
            call get_LogDetFull(size(this%Minv(iB,iS)%Mat(1,:)),this%MInv(iB,iS)%Mat,tDet)
            get_FullLogDet = tDet * get_FullLogDet
         endif
         endif
      enddo
   enddo
end function get_FullLogDet

!===============================================================================
!> Regenerate all arrays containing the possible hybridizations after a global
!! update and returns the hybrizization matrix.
function get_FullHybr_offdiag(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   integer                             :: iters,itere
   integer                             :: pos(2),b1,s1,b2,s2
   type(TOper),pointer                 :: ElementS,ElementE
   real(KINDR),pointer                 :: get_FullHybr_offdiag(:,:)

   if (allocated(this%wormContainer)) then
      allocate(get_FullHybr_offdiag((this%NOper-size(this%wormContainer))/2,&
                           (this%NOper-size(this%wormContainer))/2))
   else
      allocate(get_FullHybr_offdiag(this%NOper/2, this%NOper/2))
   end if

   get_FullHybr_offdiag(:,:)=0d0

   pos=1
   do b1=1,dstates%nbands
   do s1=1,2
   !write(*,*) "b1,s1", b1,s1

   ElementS=>this%first
   do iters=1,this%NOper

      if(elements%orbital.eq.b1.and.elements%spin.eq.s1)then !select1

      if(ElementS%has_hyb) then
         if(ElementS%CA.eq.1)then
            !write(*,*) "elementS fits b1,s1"
            !write(*,*) "elements%orbital", elements%orbital
            !write(*,*) "elements%spin", elements%spin
            !write(*,*) "elements%tau", elements%tau

            pos(2)=1
            do b2=1,dstates%nbands
            do s2=1,2
            !write(*,*) "   b2,s2",    b2,s2

            ElementE=>this%first
            do itere=1,this%NOper
               !write(*,*) "   elemente%tau", elemente%tau

               if(elemente%orbital.eq.b2.and.elemente%spin.eq.s2)then !select2
                  !write(*,*) "   elemente fits b2,s2"

               if(ElementE%has_hyb) then
                  if(ElementE%CA.eq.2)then
                     !write(*,*) "       pos(1)", pos(1)
                     !write(*,*) "       pos(2)", pos(2)
                     !write(*,*) "       get_HybrF_full(this,ElementE%tau-ElementS%tau,b1,s1,b2,s2) ",get_HybrF_full(this,ElementE%tau-ElementS%tau,b1,s1,b2,s2)  
                     get_FullHybr_offdiag(pos(2),pos(1))=&
                        get_HybrF_full(this,ElementE%tau-ElementS%tau,ElementE%orbital,ElementE%spin,ElementS%orbital,ElementS%spin)
                     pos(2)=pos(2)+1
                  endif
               endif

               endif !select2

               ElementE=>ElementE%next
            enddo

            enddo
            enddo
            pos(1)=pos(1)+1

         endif
      endif
      
      endif !select1

      ElementS=>ElementS%next
   enddo

   enddo
   enddo

end function get_FullHybr_offdiag


!===============================================================================
!> Calculates the determinant ratio and sign for the addition of 2 operators.
type(TLogDet) function get_LogDetRatPairAdd(this,MInv,Operators,Q,R,S,FPos)
!===============================================================================
! Ok so we want to also include the sign from the hopping of the a a^+ operators
! in the determinant here to account for the block form in which we write them.
! To do this my idea is to use the NOSOper array and introduce a local
! TNOSOper array. When we iterate over all the operators we can see how many
! still have to be exchanged with the current operator (Note that the 
! NOSOper array is the total number of one operator flavour of hyb operator).
 
   type(TTrace)                        :: this
   type(TSubMatrix),pointer            :: MInv(:,:)
!input
   type(TOperPointer)                  :: Operators(2)
!output
   real(KINDR), pointer, intent(out)   :: Q(:),R(:)
   real(KINDR), intent(out)            :: S
   real(KINDR)                         :: sig
   integer                             :: iQ,iR,FPos(2)
!local
   type(TOper),pointer                 :: Element,ElementS,ElementE
! generate Q,R,S for calculating the determinant ratio
   allocate(Q(this%NOSOper(Operators(1)%p%Orbital,Operators(1)%p%Spin)/2-1))
   allocate(R(this%NOSOper(Operators(1)%p%Orbital,Operators(1)%p%Spin)/2-1))  
   if(Operators(1)%p%CA.eq.1)then
      ElementS=>Operators(1)%p
      ElementE=>Operators(2)%p
   else
      ElementS=>Operators(2)%p
      ElementE=>Operators(1)%p
   endif
   
   S=get_HybrF(this,ElementE%tau-ElementS%tau,ElementE%Orbital,ElementE%Spin)
   Element=>this%ffirst(ElementS%Orbital, ElementS%Spin)%p
   iQ=1
   iR=1
   do while (associated(Element))
      if(Element%has_hyb) then
         if(Element%tau.ne.ElementS%tau)then
            if(Element%CA.ne.ElementS%CA.and..not.associated(Element,ElementE))then
               Q(iQ)=get_HybrF(this,Element%tau-ElementS%tau,Element%Orbital,Element%Spin)
               iQ=iQ+1
            elseif(Element%CA.ne.ElementE%CA.and..not.associated(Element,ElementS))then
               R(iR)=get_HybrF(this,ElementE%tau-Element%tau,Element%Orbital,Element%Spin)
               iR=iR+1
            endif
         endif
      endif
      Element=>Element%fnext
   enddo
   S = get_DetRatAdd(MInv(Operators(1)%p%Orbital,Operators(1)%p%Spin)%Mat,&
                     this%NOSOper(Operators(1)%p%Orbital,Operators(1)%p%Spin)/2-1,Q,R,S)

   !sign from inverse by partitioning (formula for last row/column)
   sig=dble(1-2*(mod(abs(FPos(1)-FPos(2)),2)))

   get_LogDetRatPairAdd = TLogDet(log = log(abs(S)),&
                                  sign = sign(1.0_KINDR, S) * sig)
end function get_LogDetRatPairAdd

!> Calculates the determinant ratio and sign for the addition of 2 operators.
type(TLogDet) function get_LogDetRatPairAdd_full(this,Operators,DStates,Q,R,S,FPos)
!===============================================================================
! Ok so we want to also include the sign from the hopping of the a a^+ operators
! in the determinant here to account for the block form in which we write them.
! To do this my idea is to use the NOSOper array and introduce a local
! TNOSOper array. When we iterate over all the operators we can see how many
! still have to be exchanged with the current operator (Note that the 
! NOSOper array is the total number of one operator flavour of hyb operator).
 
   type(TTrace)                        :: this
!input
   type(TOperPointer)                  :: Operators(2)
   type(TStates)                       :: DStates
!output
   real(KINDR), pointer, intent(out)   :: Q(:),R(:)
   real(KINDR), intent(out)            :: S
   real(KINDR)                         :: sig
   integer                             :: iQ,iR,i,FPos(2)
   integer :: b1,s1
!local
   type(TOper),pointer                 :: Element,operc,opera
! generate Q,R,S for calculating the determinant ratio
   if (allocated(this%wormContainer)) then
      allocate(Q((this%NOper-size(this%wormContainer))/2-1))
      allocate(R((this%NOper-size(this%wormContainer))/2-1))
   else
      allocate(Q(this%NOper/2-1))
      allocate(R(this%NOper/2-1))
   end if

   if(Operators(1)%p%CA.eq.1)then
      operc=>Operators(1)%p
      opera=>Operators(2)%p
   else
      operc=>Operators(2)%p
      opera=>Operators(1)%p
   endif
   
   S=get_HybrF_full(this,opera%tau-operc%tau,opera%Orbital,opera%Spin,operc%Orbital,operc%Spin)
   iQ=1
   iR=1
   !write(*,*) "detrat_add_full"

   do b1=1,dstates%nbands
   do s1=1,2
   !write(*,*) "b1,s1", b1,s1

   Element=>this%first
   do i=1,this%NOper
   !write(*,*) "   element%tau", element%tau
      if(Element%has_hyb) then
         if(Element%Orbital.eq.b1.and.Element%Spin.eq.s1&
            .and.Element%tau.ne.operc%tau)then

            if(Element%CA.ne.1.and..not.associated(Element,opera))then
               !write(*,*) "element%orbital", element%orbital
               !write(*,*) "element%spin", element%spin
               Q(iQ)=get_HybrF_full(this,Element%tau-operc%tau,operc%orbital,operc%spin,b1,s1)
               !write(*,*) "q(iq)", q(iq)
               iQ=iQ+1
            endif
         endif
      endif
      Element=>Element%next
   enddo

   Element=>this%first
   do i=1,this%NOper
      if(Element%has_hyb) then
         if(Element%Orbital.eq.b1.and.Element%Spin.eq.s1&
            .and.Element%tau.ne.opera%tau)then
            if(Element%CA.ne.2.and..not.associated(Element,operc))then
               R(iR)=get_HybrF_full(this,opera%tau-Element%tau,opera%orbital,opera%spin,b1,s1)
               iR=iR+1
            endif
         endif
      endif
      Element=>Element%next
   enddo

   enddo
   enddo
   !call Ausgabe(q,shape(q),"q",5,.false.)
   !call Ausgabe(r,shape(r),"r",5,.false.)
   !write(*,*) "s", s

   if (allocated(this%wormContainer)) then
      S = get_DetRatAdd(this%MInv_full%Mat, (this%Noper-size(this%wormContainer))/2-1, Q, R, S)
   else
      S = get_DetRatAdd(this%MInv_full%Mat, this%Noper/2-1, Q, R, S)
   end if
 
   !sign from inverse by partitioning (formula for last row/column)
   sig=dble(1-2*(mod(abs(FPos(1)-FPos(2)),2)))

   get_LogDetRatPairAdd_full = TLogDet(log = log(abs(S)),&
                                       sign = sign(1.0_KINDR, S) * sig)
end function get_LogDetRatPairAdd_full

!===============================================================================
!> Calculates the determinant ratio and sign for the operation of removing two
!! operators.
type(TLogDet) function get_LogDetRatPairRem(this,Operators,FPos,S)
!===============================================================================
   type(TTrace)                        :: this
!input
   integer                             :: FPos(2)
   type(TOperPointer)                  :: Operators(2)
!ouput
   real(KINDR), intent(out)            :: S
   real(KINDR)                         :: sig
!local
   integer                             :: iB

   if(Operators(1)%p%CA.eq.2)then
      iB=FPos(2)
      FPos(2)=FPos(1)
      FPos(2)=iB
   endif
 
   S = get_DetRatRem(this%MInv(Operators(1)%p%Orbital,Operators(1)%p%Spin)%Mat,&
                     this%NOSOper(Operators(1)%p%Orbital,Operators(1)%p%Spin)/2+1,FPos)

   !sign from inverse by partitioning (formula for last row/column)
   sig=dble(1-2*(mod(abs(FPos(1)-FPos(2)),2)))

   get_LogDetRatPairRem = TLogDet(log = log(abs(S)),&
                                  sign = sign(1.0_KINDR, S) * sig)
end function get_LogDetRatPairRem

!===============================================================================
!> Calculates the determinant ratio and sign for the operation of removing two
!! operators.
type(TLogDet) function get_LogDetRatPairRem_full(this,Operators,FPos,S)
!===============================================================================
   type(TTrace)                        :: this
!input
   integer                             :: FPos(2)
   type(TOperPointer)                  :: Operators(2)
!ouput
   real(KINDR), intent(out)            :: S
   real(KINDR)                         :: sig
!local
   integer                             :: iB

   if(Operators(1)%p%CA.eq.2)then
      iB=FPos(2)
      FPos(2)=FPos(1)
      FPos(2)=iB
   endif
 
   if (allocated(this%wormContainer)) then
      S = get_DetRatRem(this%MInv_full%Mat, (this%Noper-size(this%wormContainer))/2+1, FPos)
   else
      S = get_DetRatRem(this%MInv_full%Mat, this%Noper/2+1, FPos)
   end if

   !sign from inverse by partitioning (formula for last row/column)
   sig=dble(1-2*(mod(abs(FPos(1)-FPos(2)),2)))

   get_LogDetRatPairRem_full = TLogDet(log = log(abs(S)),&
                                       sign = sign(1.0_KINDR, S) * sig)
end function get_LogDetRatPairRem_full

!===============================================================================
!> Add offset to per-flavour indices for indexing into the full
!  hybridization matrix due to flavours occuring prior in the ordering.
subroutine add_ftau_offset(this,dstates,oper,fpos)
!===============================================================================
!input
   type(TTrace)       :: this
   type(TStates)      :: DStates
   type(TOperPointer) :: Oper(2)
!inout
   integer            :: FPos(2)
!local    
   integer            :: b1,s1
   integer            :: offset1, offset2
   
   !!! calc offset
   offset1=0
   findloop: do b1=1,dstates%nbands
   do s1=1,2
      if((oper(1)%p%orbital.eq.b1).and.(oper(1)%p%spin.eq.s1)) exit findloop
      offset1=offset1+this%noscaoper(b1,s1,Oper(1)%p%CA)
   enddo
   enddo findloop

   offset2=0
   findloop2: do b1=1,dstates%nbands
   do s1=1,2
      if((oper(2)%p%orbital.eq.b1).and.(oper(2)%p%spin.eq.s1))exit findloop2
      offset2=offset2+this%noscaoper(b1,s1,Oper(2)%p%CA)
   enddo
   enddo findloop2

   if (Oper(1)%p%CA == 1 .and. Oper(2)%p%CA == 2) then
      fpos(1)=fpos(1)+offset1
      fpos(2)=fpos(2)+offset2
   else if (Oper(1)%p%CA == 2 .and. Oper(2)%p%CA == 1) then
      fpos(1)=fpos(1)+offset2
      fpos(2)=fpos(2)+offset1
   else
      stop "add_ftau_offset: incorrect usage"
   end if
end subroutine add_ftau_offset
!===============================================================================
!> Calculates the weight ratio (up to sign equivalent to determinant ratio)
!! of the bath problem, when switching the position of a worm operator with  a
!! hybridization operator. This is essentially proposes a rank-1
!! update of the hybridization matrix. The value of the local trace stays the same, as a worm
!! operator and a hybridization operator only differ by the hybridization lines
!! connected

!! OperW is operator with no hyb lines attached
!! OperHyb is operator with hyb lines attached
!! detRat,u,vdag are the results from the matrix update
!! sgn is the sign of the proposed configuration

subroutine propose_WormReplace(this,DStates,OperW,OperHyb,u,vdag,detRat,sgn)
!===============================================================================
!input
   type(TTrace)               :: this
   type(TStates)              :: DStates
   !OperW should be a pointer from the wormContainer
   !otherwise this may bug in the sign calculation towards the end
   type(TOper),pointer        :: OperW, OperHyb
   !these should be allocated
   real(KINDR),pointer        :: u(:),vdag(:)

!local
   integer                    :: iN,iQ,iR,operPos
   integer                    :: iO, iS, offset
   integer                    :: NCnt
   real(KINDR), intent(out)   :: detRat, sgn
   type(TOper),pointer        :: Element
  
!lapack stuff
   real(KINDR),allocatable    :: res(:)
   integer                    :: mdim

   nullify(Element)
   
   !make sure these are 0
   u=0d0
   vdag=0d0

   !find position of selected operator in trace
   operPos=0
   Element=>this%first
   do iN=1,this%NOper
      if(Element%has_hyb) then
         if(Element%CA.eq.OperHyb%CA.and.&
            Element%Orbital.eq.OperHyb%Orbital.and.&
            Element%Spin.eq.OperHyb%Spin) operPos=operPos+1
      endif
      if(associated(OperHyb,Element)) exit
      Element=>Element%next
   enddo

   !find position of operator in full hybridization matrix
   offset=0
   if(this%b_offdiag) then 
      findloop: do iO=1,dstates%nbands
         do iS=1,2
            if((OperHyb%orbital.eq.iO).and.(OperHyb%Spin.eq.iS)) exit findloop
            offset=offset+this%noscaoper(iO,iS,OperHyb%CA)
         enddo
      enddo findloop
      operPos=operPos+offset
   endif
   
   !unit vector
   if(OperW%CA.eq.1) then
      vdag(operPos)=1d0
   else
      u(operPos)=1d0
   endif

   !count how many operators of the same flavor and type lie in between
   !TODO: cycle flavor lists only
   NCnt=0
   Element=>this%first
   do iN=1,this%NOper
      if(Element%has_hyb) then
         if(Element%tau>min(OperHyb%tau,OperW%tau).and.&
            Element%tau<max(OperHyb%tau,OperW%tau)) then
               if((Element%Orbital.eq.OperW%Orbital).and.&
                  (Element%Spin.eq.OperW%Spin).and.&
                  (Element%CA.eq.OperW%CA)) NCnt=NCnt+1
         end if
      end if
      Element=>Element%next
   end do

   if(OperW%tau<OperHyb%tau) NCnt=-NCnt

   !subtracting the old hyb part from our vectors, OperHyb refers to the orginal hybridization line
   if(this%b_offdiag) then
      iQ=1; iR=1
      !remember flavor ordering in full matrix 
      do iO=1,dstates%nbands
         do iS=1,2
            Element=>this%first
            do iN=1,this%NOper
               if(Element%has_hyb.and.Element%CA.ne.OperHyb%CA.and..not.associated(Element,OperHyb))then
                  if(Element%Orbital.eq.iO.and.Element%Spin.eq.iS) then
                     if(OperHyb%CA.eq.1) then   
                        u(iQ)=u(iQ)-get_HybrF_full(this,Element%tau-OperHyb%tau,Element%Orbital,Element%Spin,OperHyb%Orbital,OperHyb%Spin)
                        iQ=iQ+1
                     else
                         vdag(iR)=vdag(iR)-get_HybrF_full(this,OperHyb%tau-Element%tau,OperHyb%Orbital,OperHyb%Spin,Element%Orbital,Element%Spin)
                         iR=iR+1
                     endif
                  endif
               endif
               Element=>Element%next
            enddo
         enddo
      enddo
   else
       iQ=1; iR=1
       Element=>this%first
       do iN=1,this%NOper
          if(Element%has_hyb.and.Element%CA.ne.OperHyb%CA.and..not.associated(Element,OperHyb))then
             if(Element%Orbital.eq.OperHyb%Orbital.and.Element%Spin.eq.OperHyb%Spin)then
                !fill column for CA=1
                if(OperHyb%CA.eq.1) then   
                   u(iQ)=u(iQ)-get_HybrF(this,Element%tau-OperHyb%tau,Element%Orbital,Element%Spin)
                   iQ=iQ+1
                !fill row for CA=2
                else
                   vdag(iR)=vdag(iR)-get_HybrF(this,OperHyb%tau-Element%tau,Element%Orbital,Element%Spin)
                   iR=iR+1
                endif
             endif
          endif
          Element=>Element%next
       enddo
    endif
   
   !switch position of operators by reconnecting hyb line
   OperW%has_hyb=.true.
   OperHyb%has_hyb=.false.
      
   !adding the new hyb part to our vectors, OperW now refers to the new hybridization line
   if(this%b_offdiag) then
      iQ=1; iR=1
      !remember flavor ordering in full matrix 
      do iO=1,dstates%nbands
         do iS=1,2
            Element=>this%first
            do iN=1,this%NOper
               if(Element%has_hyb.and.Element%CA.ne.OperW%CA.and..not.associated(Element,OperW))then
                  if(Element%Orbital.eq.iO.and.Element%Spin.eq.iS) then
                     if(OperW%CA.eq.1) then   
                        u(iQ)=u(iQ)+get_HybrF_full(this,Element%tau-OperW%tau,Element%Orbital,Element%Spin,OperW%Orbital,OperW%Spin)
                        iQ=iQ+1
                     else
                         vdag(iR)=vdag(iR)+get_HybrF_full(this,OperW%tau-Element%tau,OperW%Orbital,OperW%Spin,Element%Orbital,Element%Spin)
                         iR=iR+1
                     endif
                  endif
               endif
               Element=>Element%next
            enddo
         enddo
      enddo
   else
       iQ=1; iR=1
       Element=>this%first
       do iN=1,this%NOper
          if(Element%has_hyb.and.Element%CA.ne.OperW%CA.and..not.associated(Element,OperW))then
             if(Element%Orbital.eq.OperW%Orbital.and.Element%Spin.eq.OperW%Spin)then
                !fill column for CA=1
                if(OperW%CA.eq.1) then   
                   u(iQ)=u(iQ)+get_HybrF(this,Element%tau-OperW%tau,Element%Orbital,Element%Spin)
                   iQ=iQ+1
                !fill row for CA=2
                else
                   vdag(iR)=vdag(iR)+get_HybrF(this,OperW%tau-Element%tau,Element%Orbital,Element%Spin)
                   iR=iR+1
                endif
             endif
          endif
          Element=>Element%next
       enddo
    endif

   !calculating determinant ratio
   if(this%b_offdiag) then
      !detRat = 1d0 + dot_product(vdag,matmul(this%MinV_full%Mat,u))
      mdim = size(this%MinV_full%Mat(:,1))
      allocate(res(mdim))
      call DGEMV('N',mdim,mdim,1d0,this%MinV_full%Mat,mdim,u,1,0d0,res,1)
      detRat = 1d0 + dot_product(vdag,res)
      deallocate(res)

   else
      !detRat = 1d0 + dot_product(vdag,matmul(this%MinV(OperW%Orbital,OperW%Spin)%Mat,u))
      mdim = size(this%MinV(OperW%Orbital,OperW%Spin)%Mat(:,1))
      allocate(res(mdim))
      call DGEMV('N',mdim,mdim,1d0,this%MinV(OperW%Orbital,OperW%Spin)%Mat,mdim,u,1,0d0,res,1)
      detRat = 1d0 + dot_product(vdag,res)
      deallocate(res)

   endif

   !here comes the sign calculation of the proposed configuration
   !make sure OperW is the specific wormContainer(i), otherwise get_WormSign fails
   Element=>OperW
   OperW=>OperHyb
  
   !sign of the proposed configuration
   !note we do not update the trace determinant first, thats why we do it explicititly
   sgn=-get_Sign(this, DStates)*sign(1.0_KINDR, DetRat)*((-1)**NCnt)

   !reattach the wormContainer
   OperW=>Element
   
   !revert hyb lines because we only proposed
   OperHyb%has_hyb=.true.
   OperW%has_hyb=.false.

end subroutine propose_WormReplace

!===============================================================================
!> actually carries out the replacement previously proposed by propose_WormReplace
!  updates the hybridisation matrix, the determiant and the sign of trace object

!! OperW is operator with no hyb lines attached
!! OperHyb is operator with hyb lines attached
!! detRat,u,vdag are the results from the matrix update
!! sgn is the sign of the proposed configuration

subroutine do_WormReplace(this,Dstates,OperW,OperHyb,u,vdag,detRat)
!===============================================================================
!input
   type(TTrace)               :: this
   type(TStates)              :: DStates
   !OperW should be a pointer from the wormContainer
   !otherwise this may bug in the sign calculation towards the end
   type(TOper),pointer        :: OperW, OperHyb
   !these should include the correct value for the replace
   real(KINDR),pointer        :: u(:),vdag(:)

!local
   integer                    :: iN,operPos
   integer                    :: NCnt,N
   integer                    :: iO, iS, offset
   real(KINDR)                :: detRat
   type(TOper),pointer        :: Element
!   real(KINDR),pointer        :: temp(:,:)

!lapack stuff
   real(KINDR),allocatable        :: xres(:),yres(:)
   real(KINDR),pointer            :: atemp(:,:)


   nullify(Element) !;nullify(temp)

   !number of hybridization operators of same flavor as selected
   if(this%b_offdiag) then
       N=(this%Noper-size(this%wormContainer))/2
   else
       N=this%NOSOper(OperHyb%Orbital,OperHyb%Spin)/2
   endif

   !find position of selected operator in trace
   operPos=0
   Element=>this%first
   do iN=1,this%NOper
      if(Element%has_hyb) then
         if(Element%CA.eq.OperHyb%CA.and.&
            Element%Orbital.eq.OperHyb%Orbital.and.&
            Element%Spin.eq.OperHyb%Spin) operPos=operPos+1
      endif
      if(associated(OperHyb,Element)) exit
      Element=>Element%next
   enddo
   
   !find position of operator in full hybridization matrix
   offset=0
   if(this%b_offdiag) then
      findloop: do iO=1,dstates%nbands
         do iS=1,2
            if((OperHyb%orbital.eq.iO).and.(OperHyb%Spin.eq.iS)) exit findloop
            offset=offset+this%noscaoper(iO,iS,OperHyb%CA)
         enddo
      enddo findloop
      operPos=operPos+offset
   endif

   !count how many operators of the same flavor and type lie in between
   !TODO: cycle flavor lists only
   NCnt=0
   Element=>this%first
   do iN=1,this%NOper
      if(Element%has_hyb) then
         if(Element%tau>min(OperHyb%tau,OperW%tau).and.&
            Element%tau<max(OperHyb%tau,OperW%tau)) then
               if((Element%Orbital.eq.OperW%Orbital).and.&
                  (Element%Spin.eq.OperW%Spin).and.&
                  (Element%CA.eq.OperW%CA)) NCnt=NCnt+1
         end if
      end if
      Element=>Element%next
   end do

   if(OperW%tau<OperHyb%tau) NCnt=-NCnt

   !switch position of operators by reconnecting hyb line
   OperW%has_hyb=.true.
   OperHyb%has_hyb=.false.

   !updating the matrix
   !minv using sherman morrison formula
   if(this%b_offdiag) then
      !allocate(temp(N,N))
      !temp=this%MInv_full%Mat&
      !     -matmul(this%MInv_full%Mat,matmul(spread(u,2,N)&
      !      *spread(vdag,1,N),this%MInv_full%Mat))/detRat
      allocate(atemp(N,N),xres(N),yres(N))
      atemp(:,:)=this%MInv_full%Mat(:,:)
      call DGEMV('N',N,N,1d0,this%MInv_full%Mat,N,u,1,0d0,xres,1)
      call DGEMV('T',N,N,1d0,this%MInv_full%Mat,N,vdag,1,0d0,yres,1)
      call DGER(N,N,-1./detRat,xres,1,yres,1,atemp,N)
      deallocate(xres,yres)
      deallocate(this%MInv_full%Mat)
      this%MInv_full%Mat=>atemp

   else
      !allocate(temp(N,N))
      !temp=this%MInv(OperW%Orbital,OperW%Spin)%Mat&
      !     -matmul(this%MInv(OperW%Orbital,OperW%Spin)%Mat,matmul(spread(u,2,N)&
      !      *spread(vdag,1,N),this%MInv(OperW%Orbital,OperW%Spin)%Mat))/detRat
      allocate(atemp(N,N),xres(N),yres(N))
      atemp(:,:)=this%MInv(OperW%Orbital,OperW%Spin)%Mat(:,:)
      call DGEMV('N',N,N,1d0,this%MInv(OperW%Orbital,OperW%Spin)%Mat,N,u,1,0d0,xres,1)
      call DGEMV('T',N,N,1d0,this%MInv(OperW%Orbital,OperW%Spin)%Mat,N,vdag,1,0d0,yres,1)
      call DGER(N,N,-1./detRat,xres,1,yres,1,atemp,N)
      deallocate(xres,yres)
      deallocate(this%MInv(OperW%Orbital,OperW%Spin)%Mat)
      this%MInv(OperW%Orbital,OperW%Spin)%Mat=>atemp
   endif     

   !re-sort MInv
   !we use u and vdag as temporary storage
   !advance column
   if(OperW%CA.eq.2) then
      if(this%b_offdiag) then
         u=this%MInv_full%Mat(:,operPos)
         if(NCnt>0) then
            do iN=operPos+1,operPos+NCnt
               this%MInv_full%Mat(:,iN-1)=this%MInv_full%Mat(:,iN)
            enddo
         else
            do iN=operPos-1,operPos+NCnt,-1
               this%MInv_full%Mat(:,iN+1)=this%MInv_full%Mat(:,iN)
            enddo
         endif
         this%MInv_full%Mat(:,operPos+NCnt)=u
      else
         u=this%MInv(OperW%Orbital,OperW%Spin)%Mat(:,operPos)
         if(NCnt>0) then
            do iN=operPos+1,operPos+NCnt
               this%MInv(OperW%Orbital,OperW%Spin)%Mat(:,iN-1)=this%MInv(OperW%Orbital,OperW%Spin)%Mat(:,iN)
            enddo
         else
            do iN=operPos-1,operPos+NCnt,-1
               this%MInv(OperW%Orbital,OperW%Spin)%Mat(:,iN+1)=this%MInv(OperW%Orbital,OperW%Spin)%Mat(:,iN)
            enddo
         endif
         this%MInv(OperW%Orbital,OperW%Spin)%Mat(:,operPos+NCnt)=u
      endif
   !else advance row to restore row order in hyb matrix
   else
      if(this%b_offdiag) then
         vdag=this%MInv_full%Mat(operPos,:)
         if(NCnt>0) then
            do iN=operPos+1,operPos+NCnt
               this%MInv_full%Mat(iN-1,:)=this%MInv_full%Mat(iN,:)
            enddo
         else
            do iN=operPos-1,operPos+NCnt,-1
               this%MInv_full%Mat(iN+1,:)=this%MInv_full%Mat(iN,:)
            enddo
         endif
         this%MInv_full%Mat(operPos+NCnt,:)=vdag
      else
         vdag=this%MInv(OperW%Orbital,OperW%Spin)%Mat(operPos,:)
         if(NCnt>0) then
            do iN=operPos+1,operPos+NCnt
               this%MInv(OperW%Orbital,OperW%Spin)%Mat(iN-1,:)=this%MInv(OperW%Orbital,OperW%Spin)%Mat(iN,:)
            enddo
         else
            do iN=operPos-1,operPos+NCnt,-1
               this%MInv(OperW%Orbital,OperW%Spin)%Mat(iN+1,:)=this%MInv(OperW%Orbital,OperW%Spin)%Mat(iN,:)
            enddo
         endif
         this%MInv(OperW%Orbital,OperW%Spin)%Mat(operPos+NCnt,:)=vdag
      endif
   endif

   !we save the determinant of the current configuration
   !-1^NCnt to correct sign of determinant by row/column moved to correct location
   this%Det%log = this%Det%log + log(abs(DetRat))
   this%Det%sign = this%Det%sign * sign(1.0_KINDR, DetRat) * ((-1)**NCnt)
   
   ! TODO: implement update to do this faster
   this%prewinfcountvalid = .false.
   this%NPairs = -1
   
   !update nosoper and noscaoper for offidagonal replaces
   this%NOSOper(OperW%Orbital,OperW%Spin)=this%NOSOper(OperW%Orbital,OperW%Spin)+1
   this%NOSOper(OperHyb%Orbital,OperHyb%Spin)=this%NOSOper(OperHyb%Orbital,OperHyb%Spin)-1
   this%NOSCAOper(OperW%Orbital,OperW%Spin,OperW%CA)=this%NOSCAOper(OperW%Orbital,OperW%Spin,OperW%CA)+1
   this%NOSCAOper(OperHyb%Orbital,OperHyb%Spin,OperHyb%CA)=this%NOSCAOper(OperHyb%Orbital,OperHyb%Spin,OperHyb%CA)-1
   

   !reassign worm pointers
   !make sure OperW is the specific wormContainer(i), otherwise get_WormSign fails
   OperW=>OperHyb


end subroutine do_WormReplace

!===============================================================================
!> Calculates the sign of the current configuration in the following form
!> sign = sign(Det)*sign(TraceNew)*sign(Order), where:
!>
!> sign(Det) -- Sign of the determinant:
!> This is given by the determinant of the hybridization matrix.
!>
!> sign(TraceNew) -- Sign of the current trace:
!> This is given by the local trace calculation get_Trace directly
!>
!> sign(Order) -- Sign of Block Ordering and Alternating Pairs
!> Block Ordering: We store our determinant in flavor blocks. This way we need 
!> to sort our trace into flavor blocks as well.
!> Alternating Pairs: The hybridization matrix cannot distinguish between a pair
!> c,cdag and cdag,c. This is why we require a pair-wise ordering of the form c,cdag.
!> MAKE SURE TO UPDATE this%Trace AND this%Det before calling this function
!===============================================================================
real(KINDR) function get_Sign(this,DStates)
   type(TTrace)                        :: this
   type(TStates)                       :: DStates

   if (this%b_offdiag) then
      get_Sign = this%Trace%sign&
                 * this%Det%sign&
                 * real(sign(1, recalc_sign(this)), KINDR)
   else
      get_Sign = this%Trace%sign&
                 * this%Det%sign&
                 * sign(1.0_KINDR, get_BlockSign(this,DStates))&
                 * sign(1.0_KINDR, get_permaltca(this,DStates))

   endif
   !worm sign
   if(allocated(this%wormContainer)) then
      get_Sign=get_Sign*get_WormSign(this)
   endif
end function get_Sign

integer function recalc_sign(this)
    type(TTrace), intent(in) :: this
    real(KINDR), allocatable :: operplace(:)
    real(KINDR)              ::  tmp
    type(TOper),pointer :: Element
    integer :: i, j, perm


    if (allocated(this%wormContainer)) then
       allocate(operplace(this%NOper-size(this%wormContainer)))
    else
       allocate(operplace(this%NOper))
    end if

    operplace=0
    Element => this%first
    i = 1
    do while(associated(Element))
         if(Element%has_hyb) then
             operplace(i) = (3-Element%CA) * 100000 + Element%Orbital * 10000 + &
                            Element%Spin * 1000 + Element%Tau
             i = i + 1
         endif
         Element=>Element%next
    end do

    perm = 0
    do i = 1, size(operplace)
       do j = 1, size(operplace)
          if (operplace(i) .gt. operplace(j)) then
             tmp = operplace(i)
             operplace(i) = operplace(j)
             operplace(j) = tmp
             perm = perm + 1
          end if
       end do
    end do

    do i = 1, size(operplace)/2-1
       perm = perm + i
    end do
    !write(*,*) "perm", perm
    !recalc_sign = 2 * mod(perm, 2) - 1
    recalc_sign = (-1)**(perm)
    deallocate(operplace)
end function

!===============================================================================
!> Calculates the alternating sign due to the permutation of operators in one flavor 
!! sub-block of the trace, such that operators are in order c,c^dag to obtain 
!! alternating order this was originally done in get_DetRatPairRem and get_DetRatPairAdd 
!! but since it really does belong to the local trace and is not connected to the hyb matrix
!! we write it in an extra function --> especially with regard to worms
!===============================================================================
real(KINDR) function get_permaltca(this,DStates)
    type(TTrace)                        :: this
    type(TStates)                       :: DStates
    integer                             :: TNFlav(DStates%NBands,2),permaltca
    type(TOper),pointer                 :: Element
    integer                             :: Ntmp
    
    permaltca=0
    !we consider only non-worms for the ordering procedure
    TNFlav=this%NOSOper
    Element=>this%first
    do while(associated(Element))
         if(Element%has_hyb) then
             TNFlav(Element%Orbital,Element%Spin)=TNFlav(Element%Orbital,Element%Spin)-1
             permaltca=permaltca+mod(TNFlav(Element%Orbital,Element%Spin),3-Element%CA)
         endif
         Element=>Element%next
    end do
   
    !fixing Ntmp by subtracting worms
    if(allocated(this%wormContainer)) then
       !we subtract the worms from NOper if we want to remove operators with hyb
       Ntmp=this%NOper-size(this%wormContainer)
    else
       Ntmp=this%NOper
    endif
   
    !we keep the c^dag, c ordering out of historic reasons
    !for an odd number of pairs we get an extra sign to get the required 
    !c,c^dag ordering, remember in the above we ordered c^dag,c
    if(mod(Ntmp/2,2).eq.1) permaltca=permaltca+1
   
    get_permaltca=dble((-1)**(permaltca)) 

end function get_permaltca

!===============================================================================
!> Calculates the sign necessary to sort the hyb matrix in blocks
!this function is extracted from get_DetRatPairAdd
!===============================================================================
real(KINDR) function get_BlockSign(this,DStates)
   type(TTrace)                        :: this
   type(TStates)                       :: DStates
   type(TOper),pointer                 :: Element
   integer                             :: TNFlav(DStates%NBands,2),perm

   Element=>this%first
   perm=0
   TNFlav=this%NOSOper
   do while(associated(Element))
      if(Element%has_hyb) then
         TNFlav(Element%Orbital,Element%Spin)=TNFlav(Element%Orbital,Element%Spin)-1
         if(Element%Spin.eq.2)then
            perm=sum(TNFlav(1:Element%Orbital,1))+perm
         else
            perm=sum(TNFlav(1:Element%Orbital-1,2))+perm
         endif
         perm=sum(TNFlav(1:Element%Orbital-1,Element%Spin))+perm
      endif
      Element=>Element%next
   enddo
   get_BlockSign=dble((-1)**perm)
   
   
end function get_BlockSign

!===============================================================================
!> Calculates the Worm sign
!===============================================================================
real(KINDR) function get_WormSign(this)
   type(TTrace)                        :: this
   type(TOper),pointer                 :: Element
   integer                             :: perm,i,j
   
   perm = 0
   !front sort
   Element=>this%first
   do while(associated(Element))
      if(Element%has_hyb) then
         do i=1,size(this%wormContainer)
            if(Element%tau>this%wormContainer(i)%p%tau) perm=perm+1
         enddo
      endif
      Element=>Element%next
   enddo
   
   !reverse the time order of worm operators to bring them into order c_1 c^dag_2 c_3 c^dag_4 ...
   !for IE-SIGMA the order 2314 is equivalent to 1234
   !for IE-CHI the order 231456 is equivalent to 123456
   do i=1,size(this%wormContainer)-1
      do j=i+1,size(this%wormContainer)
            if(this%wormContainer(i)%p%tau<this%wormContainer(j)%p%tau) perm=perm+1
      enddo
   enddo
   
   get_WormSign=dble((-1)**(perm))
     
end function get_WormSign

!===============================================================================
!> N creation and annihilation operators are generated at random times
!  force_diagonal controls if operators are inserted flavor-pairwise
!  operators are stored in c, cdag order
subroutine pair_OperAdd(this,NBands,Oper,N,taudiff_factor,force_diagonal)
!===============================================================================
!input
    type(TTrace)                        :: this
    integer                             :: NBands
    integer, intent(in)                 :: N
    logical                             :: force_diagonal
!output
    type(TOperPointer)                  :: Oper(N)
    real(KINDR), intent(out)            :: taudiff_factor
!local    
    real(KINDR)                         :: tau_ref, tau_min, tau_width
    type(TOper),pointer                 :: Element
    integer                             :: i

  taudiff_factor = 1._KINDR
  do i=1,size(Oper)
      Element=>Oper(i)%p
      Element%CA=mod(i,2)+1

      if (modulo(i, 2) /= 0) then
         tau_width = this%tauwin_max - this%tauwin_min
         tau_ref = this%tauwin_min + tau_width * grnd()
         Element%tau = tau_ref
         taudiff_factor = taudiff_factor * tau_width
      else
         tau_min = max(tau_ref - this%taudiff_max, this%tauwin_min)
         tau_width = min(tau_ref + this%taudiff_max, this%tauwin_max) - tau_min
         Element%tau = tau_min + tau_width * grnd()
         taudiff_factor = taudiff_factor * tau_width
      endif
      
      if((mod(i,2).eq.0).and.force_diagonal)then
         Element%Orbital=Oper(i-1)%p%Orbital
         Element%Spin=Oper(i-1)%p%Spin
      else
         Element%Orbital=randint(1, NBands)
         Element%Spin=randint(1, 2)
      endif     
  enddo
      
end subroutine pair_OperAdd

!===============================================================================
!> N creation and annihilation operators are generated at random times
!  force_diagonal controls if operators are inserted flavor-pairwise
!  operators are stored in c, cdag order
subroutine unrestricted_pair_OperAdd(beta,NBands,Oper,N,force_diagonal)
!===============================================================================
!input
    real(KINDR)                         :: beta
    integer                             :: NBands
    integer, intent(in)                 :: N
    logical                             :: force_diagonal
!output
    type(TOperPointer)                  :: Oper(N)
!local
    type(TOper),pointer                 :: Element
    integer                             :: i

  do i=1,size(Oper)
      Element=>Oper(i)%p
      Element%CA=mod(i,2)+1
      Element%CA=randint(1, 2)

      Element%tau=grnd()*beta

      if((mod(i,2).eq.0).and.force_diagonal)then
         Element%Orbital=Oper(i-1)%p%Orbital
         Element%Spin=Oper(i-1)%p%Spin
      else
         Element%Orbital=randint(1, NBands)
         Element%Spin=randint(1, 2)
      endif
  enddo
end subroutine unrestricted_pair_OperAdd

!===============================================================================
!> N creation and annihilation operators are generated at random times
!> but with predefined orbital spin
subroutine flavor_OperAdd(beta,Oper,N,Orbital,Spin)
!===============================================================================
!input    
    real(KINDR)                         :: beta
    integer, intent(in)                 :: N
    integer                             :: Orbital(N), Spin(N)
!output
    type(TOperPointer)                  :: Oper(N)
!local    
    type(TOper),pointer                 :: Element
    integer                             :: i
   
  do i=1,size(Oper)
      Element=>Oper(i)%p
      Element%CA=mod(i,2)+1
      
      Element%tau=grnd()*beta
      
      Element%Orbital=Orbital(i)
      Element%Spin=Spin(i)
  enddo
      
end subroutine flavor_OperAdd

!===============================================================================
!> N creation and annihilation operators are generated, where each pair is at random times
!> but with predefined orbital spin
subroutine density_OperAdd(beta,equal_time_offset,NBands,Oper,N,force_diagonal)
!===============================================================================
!input
    real(KINDR)                         :: beta, equal_time_offset
    integer                             :: NBands
    integer, intent(in)                 :: N
    logical                             :: force_diagonal
!output
    type(TOperPointer)                  :: Oper(N)
!local    
    type(TOper),pointer                 :: Element
    integer                             :: i
   
  do i=1,size(Oper)
      Element=>Oper(i)%p
      Element%CA=mod(i,2)+1
      
      if(mod(i,2).eq.1) then
         Element%tau=grnd()*beta
         Element%Orbital=randint(1, NBands)
         Element%Spin=randint(1, 2)
      else
         !make the density
         Element%tau=Oper(i-1)%p%tau+equal_time_offset
         if(force_diagonal) then
            Element%Orbital=Oper(i-1)%p%Orbital
            Element%Spin=Oper(i-1)%p%Spin
         else
            Element%Orbital=randint(1, NBands)
            Element%Spin=randint(1, 2)
         endif
      endif
  enddo
      
end subroutine density_OperAdd

!===============================================================================
!> sets spins and orbitals of all c c^dag operators in Oper array according to
!> scalar quantity 'flavor'
!> if 'flavor' is not supplied, random band spin patterns will be generated
!> replaces flavor_OperAdd subroutine
subroutine set_OperFlavor(Oper,N,Nbands,flavor)
!===============================================================================
!input
    integer, intent(in)                 :: Nbands, N
    integer, intent(in)                 :: flavor
!output
    type(TOperPointer)                  :: Oper(N)
!local
    type(TOper),pointer                 :: Element
    integer                             :: i, bs(N), b(N), s(N)

    !convert flavor index to band spin pattern of length N
    call index2component_general(Nbands, N, flavor, bs, b, s)

    do i=1,N
       Element=>Oper(i)%p
       Element%CA=mod(i,2)+1
    
       Element%Orbital=b(i)
       Element%Spin=s(i)
    enddo

end subroutine set_OperFlavor

!> sets spins and orbitals for all c c^dag operators according to
!> scalar quantity 'flavor'
!> preceeding the flavors corresponding to the operators is an
!> additional dangling index of the umatrix
!> the 'inner' three operators are set randomly, as they are not
!> controlled by component sampling
subroutine set_IEOperFlavor(Oper,N,Nbands,flavor,u_o,u_s,Sector)
!===============================================================================
!input
    integer, intent(in)                 :: Nbands, N
    integer, intent(in)                 :: flavor
    integer,intent(in)                  :: Sector
!output
    integer, intent(out)                 :: u_o(4),u_s(4)
    type(TOperPointer)                  :: Oper(N)
!local
    type(TOper),pointer                 :: Element
    integer                             :: i, bs(N), b(N), s(N), Nred
    
    
    if(Sector==SectorGSigma &
      .or. Sector==SectorQQ) Nred=2
    if(Sector==SectorH4 &
      .or. Sector==SectorQ4 &
      .or. Sector==SectorNQQdag &
      .or. Sector==SectorQQdd &
      .or. Sector==SectorUcaca &
      .or. Sector==SectorUccaa &
      .or. Sector==SectorQUDdag) Nred=4
    
    !convert flavor index to band spin pattern of length Nred
    call index2component_general(Nbands, Nred, flavor, bs, b, s)

    !set to out-of-bounds to throw error for wrong use
    u_o=0
    u_s=0

    if(Sector==SectorGSigma &
      .or. Sector==SectorH4) then 
    !only first index belongs to umatrix
       u_o(1)=b(1)
       u_s(1)=s(1)
    elseif(Sector==SectorQQ &
      .or. Sector==SectorQ4 &
      .or. Sector==SectorNQQdag &
      .or. Sector==SectorQQdd &
      .or. Sector==SectorUcaca &
      .or. Sector==SectorUccaa &
      .or. Sector==SectorQUDdag) then
    !all free indices belong to umatrix
       do i=1,Nred
          u_o(i)=b(i)
          u_s(i)=s(i)
       enddo
    endif

    do i=1,N
       Element=>Oper(i)%p
       Element%CA=mod(i,2)+1
    
       !inner operators of symmetric ie
       if(Sector==SectorQQ &
         .or. Sector==SectorQ4 &
         .or. Sector==SectorNQQdag &
         .or. Sector==SectorQQdd &
         .or. Sector==SectorUcaca &
         .or. Sector==SectorUccaa &
         .or. Sector==SectorQUDdag) then 

         Element%Orbital=randint(1, Nbands)
         Element%Spin=randint(1, 2)

       !inner operators of asymmetric ie
       elseif(Sector==SectorGSigma &
         .or. Sector==SectorH4) then

         if(i.le.3) then 
            Element%Orbital=randint(1, Nbands)
            Element%Spin=randint(1, 2)
          !outer operators of asymetric ie
          else
             Element%Orbital=b(i-2)
             Element%Spin=s(i-2)
          endif

       else
          stop 'attempting to set IE times for non-ie sector'
       endif
    enddo

end subroutine set_IEOperFlavor

!===============================================================================
!> sets band spinpattern in Oper to a pairwise bandspin pattern for c c^dag ...
!> to be used after set_OperFlavor call
subroutine set_OperDiagonal(Oper,N)
!===============================================================================
!input
   integer, intent(in)                 :: N
!output
   type(TOperPointer)                  :: Oper(N)
!local
   integer                             :: i
  
   do i=2,size(Oper),2
      Oper(i)%p%Orbital=Oper(i-1)%p%Orbital
      Oper(i)%p%Spin=Oper(i-1)%p%Spin
   enddo

end subroutine set_OperDiagonal

!===============================================================================
!> sets imaginary times of all c c^dag operators in Oper array between 0 and beta
subroutine set_OperTime(Oper,N,beta,equal_time_offset,Sector)
!===============================================================================
!input
   integer, intent(in)                 :: N, Sector
   real(KINDR), intent(in)             :: beta, equal_time_offset
!output
   type(TOperPointer)                  :: Oper(N)
!local
   type(TOper), pointer                :: Element
   integer                             :: i

   !N operators with N random times 
   !--> Sector 1,2,4 (Z,1P-GF,2P-GF)
   if(Sector .eq. SectorZ &
     .or. Sector .eq. SectorG &
     .or. Sector .eq. SectorG4) then

      do i=1,N
         Element=>Oper(i)%p
         Element%tau=grnd()*beta
         if(i>1) then
            if(Oper(i-1)%p%tau.eq.Element%tau) then 
               Element%tau=Element%tau+equal_time_offset
            endif
         endif
      enddo
  
   !N operators with N-2 random times
   !--> Sector 3,5 (1P-IE,2P-IE)
   elseif(Sector .eq. SectorGSigma .or. Sector .eq. SectorH4) then

      do i=1,N
         Element=>Oper(i)%p
         Element%tau=grnd()*beta
      enddo

      !catching rare events close to beginning of trace      
      if(Oper(2)%p%tau .lt. 2d0*equal_time_offset) Oper(2)%p%tau=2d0*equal_time_offset

      !equal time object
      Oper(3)%p%tau=Oper(2)%p%tau-equal_time_offset
      Oper(1)%p%tau=Oper(2)%p%tau-2d0*equal_time_offset

   !4 operators with 2 random times
   !--> Sector 6,7 (P2,P2PP)
   elseif(Sector .eq. SectorP2 .or. Sector .eq. SectorP2pp) then

      Oper(1)%p%tau=grnd()*beta
      Oper(4)%p%tau=grnd()*beta
      
      !catching rare events close to beginning and end of trace      
      if(Oper(1)%p%tau .lt. equal_time_offset) Oper(1)%p%tau=equal_time_offset
      if(Oper(4)%p%tau .gt. beta-equal_time_offset) Oper(4)%p%tau=beta-equal_time_offset

      if(Sector .eq. SectorP2) then
         Oper(2)%p%tau=Oper(1)%p%tau-equal_time_offset
         Oper(3)%p%tau=Oper(4)%p%tau+equal_time_offset
      else
         Oper(2)%p%tau=Oper(4)%p%tau+equal_time_offset
         Oper(3)%p%tau=Oper(1)%p%tau-equal_time_offset
      endif
      
   !4 operators with 3 random times
   !--> Sector 8,9 (P3,P3PP)
   elseif(Sector .eq. SectorP3 .or. Sector .eq. SectorP3pp) then

      Oper(1)%p%tau=grnd()*beta
      Oper(4)%p%tau=grnd()*beta
      
      !catching rare events close to end of trace      
      if(Oper(4)%p%tau .gt. beta-equal_time_offset) Oper(4)%p%tau=beta-equal_time_offset

      if(Sector .eq. SectorP3) then 
         Oper(2)%p%tau=grnd()*beta
         Oper(3)%p%tau=Oper(4)%p%tau+equal_time_offset
      else
         Oper(2)%p%tau=Oper(4)%p%tau+equal_time_offset
         Oper(3)%p%tau=grnd()*beta
      endif

   elseif(Sector .eq. SectorQQ) then

      do i=1,N
         Element=>Oper(i)%p
         Element%tau=grnd()*beta
      enddo

      !catching rare events close to beginning of trace      
      if(Oper(2)%p%tau .lt. 2d0*equal_time_offset) Oper(2)%p%tau=2d0*equal_time_offset
      
      if(Oper(6)%p%tau .lt. 2d0*equal_time_offset) Oper(6)%p%tau=2d0*equal_time_offset

      !equal time object 1
      Oper(3)%p%tau=Oper(2)%p%tau-equal_time_offset
      Oper(1)%p%tau=Oper(2)%p%tau-2d0*equal_time_offset
      
      !equal time object 2
      Oper(4)%p%tau=Oper(6)%p%tau-equal_time_offset
      Oper(5)%p%tau=Oper(6)%p%tau-2d0*equal_time_offset
    
   elseif(Sector .eq. SectorQ4) then

      do i=1,N
         Element=>Oper(i)%p
         Element%tau=grnd()*beta
      enddo

      !catching rare events close to beginning of trace      
      if(Oper(2)%p%tau .lt. 2d0*equal_time_offset) Oper(2)%p%tau=2d0*equal_time_offset
      if(Oper(6)%p%tau .lt. 2d0*equal_time_offset) Oper(6)%p%tau=2d0*equal_time_offset
      if(Oper(8)%p%tau .lt. 2d0*equal_time_offset) Oper(8)%p%tau=2d0*equal_time_offset
      if(Oper(12)%p%tau .lt. 2d0*equal_time_offset) Oper(12)%p%tau=2d0*equal_time_offset

      !equal time object 1
      Oper(3)%p%tau=Oper(2)%p%tau-equal_time_offset
      Oper(1)%p%tau=Oper(2)%p%tau-2d0*equal_time_offset
      
      !equal time object 2
      Oper(4)%p%tau=Oper(6)%p%tau-equal_time_offset
      Oper(5)%p%tau=Oper(6)%p%tau-2d0*equal_time_offset
      
      !equal time object 3
      Oper(9)%p%tau=Oper(8)%p%tau-equal_time_offset
      Oper(7)%p%tau=Oper(8)%p%tau-2d0*equal_time_offset
      
      !equal time object 4
      Oper(10)%p%tau=Oper(12)%p%tau-equal_time_offset
      Oper(11)%p%tau=Oper(12)%p%tau-2d0*equal_time_offset
    
    elseif(Sector .eq. SectorNQQdag) then

      Oper(2)%p%tau = grnd() * beta
      Oper(4)%p%tau = grnd() * beta
      Oper(8)%p%tau = grnd() * beta

      !catching rare events close to beginning of trace      
      if(Oper(2)%p%tau .lt. 2d0*equal_time_offset) Oper(2)%p%tau=2d0*equal_time_offset
      if(Oper(4)%p%tau .lt. 2d0*equal_time_offset) Oper(4)%p%tau=2d0*equal_time_offset
      if(Oper(8)%p%tau .lt. equal_time_offset) Oper(8)%p%tau=equal_time_offset

      !equal time object 1
      Oper(3)%p%tau=Oper(2)%p%tau-equal_time_offset
      Oper(1)%p%tau=Oper(2)%p%tau-2d0*equal_time_offset
      
      !equal time object 2
      Oper(6)%p%tau=Oper(4)%p%tau-equal_time_offset
      Oper(5)%p%tau=Oper(4)%p%tau-2d0*equal_time_offset
      
      !equal time object 3
      Oper(7)%p%tau=Oper(8)%p%tau - equal_time_offset ! or maybe the other way round?
      !  7:  annih          8: creat

    elseif(Sector .eq. SectorQQdd) then

      Oper(2)%p%tau = grnd() * beta
      Oper(4)%p%tau = grnd() * beta
      Oper(6)%p%tau = grnd() * beta


      !catching rare events close to beginning of trace      
      if(Oper(2)%p%tau .lt. 2d0*equal_time_offset) Oper(2)%p%tau=2d0*equal_time_offset
      if(Oper(4)%p%tau .lt. 2d0*equal_time_offset) Oper(4)%p%tau=2d0*equal_time_offset
      if(Oper(6)%p%tau .lt. equal_time_offset) Oper(6)%p%tau=equal_time_offset

      !equal time object 1
      Oper(3)%p%tau=Oper(2)%p%tau-equal_time_offset
      Oper(1)%p%tau=Oper(2)%p%tau-2d0*equal_time_offset
      
      !equal time object 2
      Oper(7)%p%tau=Oper(4)%p%tau-equal_time_offset
      Oper(5)%p%tau=Oper(4)%p%tau-2d0*equal_time_offset
      
      !equal time object 3
      Oper(8)%p%tau=Oper(6)%p%tau - equal_time_offset 

    elseif(Sector .eq. SectorUcaca) then

      Oper(2)%p%tau = grnd() * beta
      Oper(4)%p%tau = grnd() * beta

      ! catching rare events close to beginning of trace
      if(Oper(2)%p%tau .lt. equal_time_offset) Oper(2)%p%tau = equal_time_offset
      if(Oper(4)%p%tau .lt. equal_time_offset) Oper(4)%p%tau = equal_time_offset

      ! equal time object 1
      Oper(1)%p%tau = Oper(2)%p%tau - equal_time_offset

      ! equal time object 2
      Oper(3)%p%tau = Oper(4)%p%tau - equal_time_offset

    elseif(Sector .eq. SectorUccaa) then

      Oper(2)%p%tau = grnd() * beta
      Oper(3)%p%tau = grnd() * beta

      ! catching rare events close to beginning of trace
      if(Oper(2)%p%tau .lt. equal_time_offset) Oper(2)%p%tau = equal_time_offset
      if(Oper(3)%p%tau .lt. equal_time_offset) Oper(3)%p%tau = equal_time_offset

      ! equal time object 1
      Oper(1)%p%tau = Oper(3)%p%tau - equal_time_offset

      ! equal time object 2
      Oper(4)%p%tau = Oper(2)%p%tau - equal_time_offset
      
   elseif(Sector .eq. SectorQUDdag) then
     
     Oper(2)%p%tau = grnd() * beta
     Oper(4)%p%tau = grnd() * beta

      ! catching rare events close to beginning of trace
     if(Oper(2)%p%tau .lt. 2d0*equal_time_offset) Oper(2)%p%tau = 2d0*equal_time_offset
     if(Oper(4)%p%tau .lt. equal_time_offset) Oper(4)%p%tau = equal_time_offset

      ! equal time object 1
     Oper(3)%p%tau = Oper(2)%p%tau - equal_time_offset
     Oper(1)%p%tau = Oper(2)%p%tau - 2d0*equal_time_offset

   endif



end subroutine set_OperTime

function array_nonzero_elements(array)
  implicit none
  real(KINDR) :: array(:)
  integer     :: in1,array_nonzero_elements

  array_nonzero_elements = 0

  do in1 = 1,size(array) 
     if (array(in1) /= 0.0_KINDR ) array_nonzero_elements = array_nonzero_elements + 1
  end do

  return
end function



!===============================================================================
!> check for rare cases where there is an operator between equal time operators
logical function check_EqualTime(this,Sector)
!===============================================================================
!input
 type(TTrace)                            :: this
 integer, intent(in)                     :: Sector
   
   check_EqualTime = .false. 

   !no equal time in Z and GF spaces
   if(Sector .eq. SectorZ &
     .or. Sector .eq. SectorG &
     .or. Sector .eq. SectorG4) then
      check_EqualTime = .true.

   !N operators with N-2 random times
   !--> Sector 3,5 (1P-IE,2P-IE)
   elseif(Sector .eq. SectorGSigma .or. Sector .eq. SectorH4) then
      
       if((associated(this%wormContainer(2)%p%prev,this%wormContainer(3)%p)).and.&
          (associated(this%wormContainer(3)%p%prev,this%wormContainer(1)%p))) check_EqualTime = .true.

   !4 operators with 2 random times
   !--> Sector 6,7 (P2,P2PP)
   elseif(Sector .eq. SectorP2) then
   
       if((associated(this%wormContainer(1)%p%prev,this%wormContainer(2)%p)).and.&
          (associated(this%wormContainer(3)%p%prev,this%wormContainer(4)%p))) check_EqualTime = .true.
   
   elseif(Sector .eq. SectorP2pp) then

       if((associated(this%wormContainer(1)%p%prev,this%wormContainer(3)%p)).and.&
          (associated(this%wormContainer(2)%p%prev,this%wormContainer(4)%p))) check_EqualTime = .true.
   
   !4 operators with 3 random times
   !--> Sector 8,9 (P3,P3PP)
   elseif(Sector .eq. SectorP3) then

       if(associated(this%wormContainer(3)%p%prev,this%wormContainer(4)%p)) check_EqualTime = .true.
   
   elseif(Sector .eq. SectorP3pp) then
       
       if(associated(this%wormContainer(2)%p%prev,this%wormContainer(4)%p)) check_EqualTime = .true.
   
   elseif(Sector .eq. SectorQQ) then
      
       if((associated(this%wormContainer(2)%p%prev,this%wormContainer(3)%p)).and.&
          (associated(this%wormContainer(3)%p%prev,this%wormContainer(1)%p)).and.&
          (associated(this%wormContainer(6)%p%prev,this%wormContainer(4)%p)).and.&
          (associated(this%wormContainer(4)%p%prev,this%wormContainer(5)%p))) check_EqualTime = .true.

   elseif(Sector .eq. SectorQ4) then
      !TODO: make this faster
       if((associated(this%wormContainer(2)%p%prev,this%wormContainer(3)%p)).and.&
          (associated(this%wormContainer(3)%p%prev,this%wormContainer(1)%p)).and.&
          (associated(this%wormContainer(6)%p%prev,this%wormContainer(4)%p)).and.&
          (associated(this%wormContainer(4)%p%prev,this%wormContainer(5)%p)).and.&
          (associated(this%wormContainer(8)%p%prev,this%wormContainer(9)%p)).and.&
          (associated(this%wormContainer(9)%p%prev,this%wormContainer(7)%p)).and.&
          (associated(this%wormContainer(12)%p%prev,this%wormContainer(10)%p)).and.&
          (associated(this%wormContainer(10)%p%prev,this%wormContainer(11)%p))) check_EqualTime = .true.
    elseif(Sector .eq. SectorNQQdag) then
      if((associated(this%wormContainer(2)%p%prev, this%wormContainer(3)%p)).and.&
         (associated(this%wormContainer(3)%p%prev, this%wormContainer(1)%p)).and.&
         (associated(this%wormContainer(4)%p%prev, this%wormContainer(6)%p)).and.&
         (associated(this%wormContainer(6)%p%prev, this%wormContainer(5)%p)).and.&
         (associated(this%wormContainer(8)%p%prev, this%wormContainer(7)%p))) check_EqualTime = .true.
    elseif(Sector .eq. SectorQQdd) then
      if((associated(this%wormContainer(2)%p%prev, this%wormContainer(3)%p)).and.&
         (associated(this%wormContainer(3)%p%prev, this%wormContainer(1)%p)).and.&
         (associated(this%wormContainer(4)%p%prev, this%wormContainer(7)%p)).and.&
         (associated(this%wormContainer(7)%p%prev, this%wormContainer(5)%p)).and.&
         (associated(this%wormContainer(6)%p%prev, this%wormContainer(8)%p))) check_EqualTime = .true.
     elseif(Sector .eq. SectorUcaca) then
       if((associated(this%wormContainer(2)%p%prev, this%wormContainer(1)%p)).and.&
          (associated(this%wormContainer(4)%p%prev, this%wormContainer(3)%p))) check_EqualTime = .true.
     elseif(Sector .eq. SectorUccaa) then
       if((associated(this%wormContainer(2)%p%prev, this%wormContainer(4)%p)).and.&
          (associated(this%wormContainer(3)%p%prev, this%wormContainer(1)%p))) check_EqualTime = .true.
     elseif(Sector .eq. SectorQUDdag) then
       if((associated(this%wormContainer(2)%p%prev, this%wormContainer(3)%p)).and. &
          (associated(this%wormContainer(3)%p%prev, this%wormContainer(1)%p))) check_EqualTime = .true.
   endif

   if (.not. check_EqualTime) write(*,*) 'equal time check yields FALSE'

end function check_EqualTime

!===============================================================================
!> Operators previously generated in the Oper array are time-ordered 
!  and their position within the trace is determined
!  May only be used for insertions inside of the window.
logical function process_OperAdd(this,Oper,FPos,N)
!===============================================================================
type(TTrace)                           :: this
!input
   integer                             :: N
!output
   type(TOperPointer)                  :: Oper(N)
   integer                             :: FPos(2)
!local
   integer                             :: i,j,TFPos(2)
   type(TOper),pointer                 :: Element,mEl

   process_OperAdd = .true.
   !time order newly generated operators
   do i=1,size(Oper)
      if(Oper(i)%p%tau.gt.this%beta) stop "error tau gt beta"
      do j=i+1,size(Oper)
         if(Oper(i)%p%tau.eq.Oper(j)%p%tau .and. (Oper(i)%p%has_hyb .or. Oper(j)%p%has_hyb)) then
            process_OperAdd = .false.
            return
         end if
         if(Oper(i)%p%tau>Oper(j)%p%tau) then
            Element=>Oper(i)%p
            Oper(i)%p=>Oper(j)%p
            Oper(j)%p=>Element
         endif
      enddo
   enddo

   ! generate initial links between operators:
   ! link ops to be added in time order so these links can be retained
   ! if they happen to land adjacently
   do i=1,size(Oper)
      if(i.gt.1)then
         Oper(i)%p%prev=>Oper(i-1)%p
         j = i - 1
         do
            if (Oper(j)%p%Orbital == Oper(i)%p%Orbital .and. Oper(j)%p%Spin == Oper(i)%p%Spin) then
               Oper(i)%p%fprev => Oper(j)%p
               exit
            end if
            j = j - 1
            if (j <= 0) then
               Oper(i)%p%fprev => null()
               exit
            end if
         end do
      else
         Oper(i)%p%prev=>null()
         Oper(i)%p%fprev => null()
      endif
      if(i.lt.size(Oper))then
         Oper(i)%p%next=>Oper(i+1)%p
         j = i + 1
         do
            if (Oper(j)%p%Orbital == Oper(i)%p%Orbital .and. Oper(j)%p%Spin == Oper(i)%p%Spin) then
               Oper(i)%p%fnext => Oper(j)%p
               exit
            end if
            j = j + 1
            if (j > size(Oper)) then
               Oper(i)%p%fnext => null()
               exit
            end if
         end do
      else
         Oper(i)%p%next=>null()
         Oper(i)%p%fnext => null()
      endif
   enddo

   i=1
   ! start at the first operator after the lower window edge
   call ensure_valid_PreWinFCount(this)
   TFPos = this%PreWinFCount(Oper(i)%p%Orbital, Oper(i)%p%Spin, :)
   if (associated(this%fprewin(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p)) then
      Element => this%fprewin(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p%fnext
   else
      Element => this%ffirst(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p
   end if
   
   ! find position in the list by going forward through flavour-specific list (wide stride) and
   ! then backward through the intervals in between
   do while (i <= size(Oper))
      ! go through the flavour-list and count until we have passed the insertion position
      if (associated(Element)) then
         if (Oper(i)%p%tau == Element%tau .and. (Element%has_hyb .or. Oper(i)%p%has_hyb)) then
            process_OperAdd = .false.
            return
         end if
         if (Oper(i)%p%tau >= Element%tau) then
            !make sure FPOS stays correct for hyb function
            !by only considering operators with hybs
            if (Element%has_hyb) TFPos(Element%CA) = TFPos(Element%CA) + 1
            Element => Element%fnext
            cycle
         end if
      end if

      TFPos(Oper(i)%p%CA) = TFPos(Oper(i)%p%CA) + 1
      FPos(Oper(i)%p%CA) = TFPos(Oper(i)%p%CA)

      ! set next and previous flavour-list element for new operator
      if (associated(Element)) then
         if (.not. associated(Oper(i)%p%fnext)) then
            Oper(i)%p%fnext => Element
         else if (Oper(i)%p%fnext%tau >= Element%tau) then
            Oper(i)%p%fnext => Element
         end if

         if (.not. associated(Oper(i)%p%fprev)) then
            Oper(i)%p%fprev => Element%fprev
         else if (associated(Element%fprev)) then
            if (Oper(i)%p%fprev%tau < Element%fprev%tau) Oper(i)%p%fprev => Element%fprev
         end if
         mEl => Element%prev
      else
         if (.not. associated(Oper(i)%p%fprev)) then
            Oper(i)%p%fprev => this%flast(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p
         else if (associated(this%flast(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p)) then
            if (Oper(i)%p%fprev%tau < this%flast(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p%tau)&
               Oper(i)%p%fprev => this%flast(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p
         end if
         mEl => this%last
      end if

      ! go back through the list past the insertion position
      do
         if (.not. associated(mEl)) exit
         if (mEl%tau <= Oper(i)%p%tau) then
            if (associated(mEl) .and. Oper(i)%p%tau == mEl%tau .and. (mEl%has_hyb .or. Oper(i)%p%has_hyb)) then
               process_OperAdd = .false.
               return
            end if
            exit
         else
            mEl => mEl%prev
         end if
      end do

      ! set previous and next list element for the new operator
      if (.not. associated(Oper(i)%p%prev)) then
         Oper(i)%p%prev => mEl
      else if (associated(mEl)) then
         if (Oper(i)%p%prev%tau < mEl%tau)&
            Oper(i)%p%prev => mEl
      end if
      if (.not. associated(Oper(i)%p%next)) then
         if (associated(mEl)) then
            Oper(i)%p%next => mEl%next
         else
            Oper(i)%p%next => this%first
         end if
      else
         if (associated(mEl)) then
            if (associated(mEl%next)) then
               if (Oper(i)%p%next%tau >= mEl%next%tau)&
                  Oper(i)%p%next => mEl%next
            end if
         else
            if (associated(this%first)) then
               if (Oper(i)%p%next%tau >= this%first%tau)&
                  Oper(i)%p%next => this%first
            end if
         end if
      end if

      if (i < size(Oper)) then
         ! keep count and position if the next operator has the same flavour

         ! FIXME: If more than two operators are inserted and there
         ! are two operators of the same flavour that are separated in
         ! tau by one of different flavour, the count is not the same
         ! as with the original method. Check if this needs to be
         ! corrected.
         if (Oper(i)%p%Orbital /= Oper(i+1)%p%Orbital&
             .or. Oper(i)%p%Spin /= Oper(i+1)%p%Spin) then
            TFPos = this%PreWinFCount(Oper(i+1)%p%Orbital, Oper(i+1)%p%Spin, :)
            if (associated(this%fprewin(Oper(i+1)%p%Orbital, Oper(i+1)%p%Spin)%p)) then
               Element => this%fprewin(Oper(i+1)%p%Orbital, Oper(i+1)%p%Spin)%p%fnext
            else
               Element => this%ffirst(Oper(i+1)%p%Orbital, Oper(i+1)%p%Spin)%p
            end if
         end if
      end if
      i = i + 1
   enddo
end function process_OperAdd

! !===============================================================================
! !> Operators previously generated in the Oper array are time-ordered 
! !  and their position within the trace is determined
! subroutine process_OperAdd_global(this,DStates,Oper,FPos,Pos,N)
! !===============================================================================
! type(TTrace)                           :: this
! !input
!    type(TStates)                       :: DStates
!    integer                             :: N
! !output
!    type(TOperPointer)                  :: Oper(N)
!    integer                             :: FPos(2),Pos(N)
! !local
!    integer                             :: i,j,k,TiSSt(1:size(this%States(:,1))),iOper
!    type(TOper),pointer                 :: Element
!    integer                             :: TFPos(DStates%Nbands,2,2)
   
!    TiSSt=this%States(:,1)
   
!    !time order newly generated operators
!    do i=1,size(Oper)
!       if(Oper(i)%p%tau.gt.this%beta) stop "error tau gt beta"
!       do j=i+1,size(Oper)
!          if(Oper(i)%p%tau.eq.Oper(j)%p%tau) Oper(j)%p%tau=Oper(j)%p%tau+this%equal_time_offset
!          if(Oper(i)%p%tau>Oper(j)%p%tau) then
!             Element=>Oper(i)%p
!             Oper(i)%p=>Oper(j)%p
!             Oper(j)%p=>Element
!          endif
!       enddo
!    enddo

!    ! generate initial links between operators
!    do i=1,size(Oper)
!       if(i.gt.1)then
!          Oper(i)%p%prev=>Oper(i-1)%p
!          j = i - 1
!          do
!             if (Oper(j)%p%Orbital == Oper(i)%p%Orbital .and. Oper(j)%p%Spin == Oper(i)%p%Spin) then
!                Oper(i)%p%fprev => Oper(j)%p
!                exit
!             end if
!             j = j - 1
!             if (j <= 0) then
!                Oper(i)%p%fprev => null()
!                exit
!             end if
!          end do
!       else
!          Oper(i)%p%prev=>null()
!          Oper(i)%p%fprev => null()
!       endif
!       if(i.lt.size(Oper))then
!          Oper(i)%p%next=>Oper(i+1)%p
!          j = i + 1
!          do
!             if (Oper(j)%p%Orbital == Oper(i)%p%Orbital .and. Oper(j)%p%Spin == Oper(i)%p%Spin) then
!                Oper(i)%p%fnext => Oper(j)%p
!                exit
!             end if
!             j = j + 1
!             if (j > size(Oper)) then
!                Oper(i)%p%fnext => null()
!                exit
!             end if
!          end do
!       else
!          Oper(i)%p%next=>null()
!          Oper(i)%p%fnext => null()
!       endif
!    enddo
!    i=1
!    TFPos=0
!    Element=>this%first
!    iOper=1
   
!    ! find position in the list
!    do while(i.lt.size(Oper)+1.and.associated(Element))
!       if(Element%has_hyb) then
!          if(Oper(i)%p%tau.eq.Element%tau)Oper(i)%p%tau=Oper(i)%p%tau+this%equal_time_offset
!       endif
!       if(Oper(i)%p%tau.lt.Element%tau)then
!          TFPos(Oper(i)%p%Orbital,Oper(i)%p%Spin,Oper(i)%p%CA)=TFPos(Oper(i)%p%Orbital,Oper(i)%p%Spin,Oper(i)%p%CA)+1
!          FPos(Oper(i)%p%CA)=TFPos(Oper(i)%p%Orbital,Oper(i)%p%Spin,Oper(i)%p%CA)
!          Oper(i)%p%next=>Element
!          if (Oper(i)%p%Orbital == Element%Orbital .and. Oper(i)%p%Spin == Element%Spin)&
!             Oper(i)%p%fnext => Element
!          if(i.gt.1)then
!             if(associated(Oper(i)%p%prev,Oper(i-1)%p))Oper(i-1)%p%next=>Oper(i)%p
!             if(associated(Oper(i)%p%fprev,Oper(i-1)%p))&
!                Oper(i-1)%p%fnext=>Oper(i)%p
!          endif
!          i=i+1
!          iOper=iOper+1
!          if(iOper.le.size(Pos))then
!             Pos(iOper)=Pos(iOper-1)
!             Pos(iOper)=Pos(iOper)+1
!          endif
!       else
!          !make sure FPOS stays correct for hyb function
!          !by only considering operators with hybs
!          if(Element%has_hyb) then
!             TFPos(Element%Orbital,Element%Spin,Element%CA)=TFPos(Element%Orbital,Element%Spin,Element%CA)+1
!          endif
!          Oper(i)%p%prev=>Element
!          if (Oper(i)%p%Orbital == Element%Orbital .and. Oper(i)%p%Spin == Element%Spin)&
!             Oper(i)%p%fprev=>Element
!          Element=>Element%next
!       endif
!    enddo
!    ! set the last flavour positions
!    do while(i.lt.size(Oper)+1)
!       TFPos(Oper(i)%p%Orbital,Oper(i)%p%Spin,Oper(i)%p%CA)=TFPos(Oper(i)%p%Orbital,Oper(i)%p%Spin,Oper(i)%p%CA)+1
!       FPos(Oper(i)%p%CA)=TFPos(Oper(i)%p%Orbital,Oper(i)%p%Spin,Oper(i)%p%CA)
!       i=i+1
!       iOper=iOper+1
!       if(iOper.le.size(Pos))then
!          Pos(iOper)=Pos(iOper-1)
!          Pos(iOper)=Pos(iOper)+1
!       endif
!    enddo
! end subroutine process_OperAdd_global

!===============================================================================
!> Operators previously generated in the Oper array are time-ordered 
!  and their position within the trace is determined
subroutine process_OperAdd_global(this,Oper,FPos,N)
!===============================================================================
type(TTrace)                           :: this
!input
   integer                             :: N
!output
   type(TOperPointer)                  :: Oper(N)
   integer                             :: FPos(2)
!local
   integer                             :: i,j,TFPos(2)
   type(TOper),pointer                 :: Element,mEl
   
   !time order newly generated operators
   do i=1,size(Oper)
      if(Oper(i)%p%tau.gt.this%beta) stop "error tau gt beta"
      do j=i+1,size(Oper)
         if(Oper(i)%p%tau.eq.Oper(j)%p%tau) Oper(j)%p%tau=Oper(j)%p%tau+this%equal_time_offset
         if(Oper(i)%p%tau>Oper(j)%p%tau) then
            Element=>Oper(i)%p
            Oper(i)%p=>Oper(j)%p
            Oper(j)%p=>Element
         endif
      enddo
   enddo

   ! generate initial links between operators:
   ! link ops to be added in time order so these links can be retained
   ! if they happen to land adjacently
   do i=1,size(Oper)
      if(i.gt.1)then
         Oper(i)%p%prev=>Oper(i-1)%p
         j = i - 1
         do
            if (Oper(j)%p%Orbital == Oper(i)%p%Orbital .and. Oper(j)%p%Spin == Oper(i)%p%Spin) then
               Oper(i)%p%fprev => Oper(j)%p
               exit
            end if
            j = j - 1
            if (j <= 0) then
               Oper(i)%p%fprev => null()
               exit
            end if
         end do
      else
         Oper(i)%p%prev=>null()
         Oper(i)%p%fprev => null()
      endif
      if(i.lt.size(Oper))then
         Oper(i)%p%next=>Oper(i+1)%p
         j = i + 1
         do
            if (Oper(j)%p%Orbital == Oper(i)%p%Orbital .and. Oper(j)%p%Spin == Oper(i)%p%Spin) then
               Oper(i)%p%fnext => Oper(j)%p
               exit
            end if
            j = j + 1
            if (j > size(Oper)) then
               Oper(i)%p%fnext => null()
               exit
            end if
         end do
      else
         Oper(i)%p%next=>null()
         Oper(i)%p%fnext => null()
      endif
   enddo

   i=1
   TFPos=0
   Element => this%ffirst(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p
   
   ! find position in the list by going forward through flavour-specific list (wide stride) and
   ! then backward through the intervals in between
   do while (i <= size(Oper))
      ! second nested if clause is non-short-circuit version of
      ! if (.not. associated(Element) .or. (associated(Element) .and. Oper(i)%p%tau < Element%tau))
      ! by prepending what used to be the else clause and cycling at its end
      if (associated(Element)) then
         if (Oper(i)%p%tau == Element%tau)&
            Oper(i)%p%tau = Oper(i)%p%tau + this%equal_time_offset
         if (Oper(i)%p%tau >= Element%tau) then
            !make sure FPOS stays correct for hyb function
            !by only considering operators with hybs
            if (Element%has_hyb) TFPos(Element%CA) = TFPos(Element%CA) + 1
            Element => Element%fnext
            cycle
         end if
      end if
      TFPos(Oper(i)%p%CA) = TFPos(Oper(i)%p%CA) + 1
      FPos(Oper(i)%p%CA) = TFPos(Oper(i)%p%CA)

      if (associated(Element)) then
         if (.not. associated(Oper(i)%p%fnext)) then
            Oper(i)%p%fnext => Element
         else if (Oper(i)%p%fnext%tau >= Element%tau) then
            Oper(i)%p%fnext => Element
         end if

         if (.not. associated(Oper(i)%p%fprev)) then
            Oper(i)%p%fprev => Element%fprev
         else if (associated(Element%fprev)) then
            if (Oper(i)%p%fprev%tau < Element%fprev%tau) Oper(i)%p%fprev => Element%fprev
         end if
         mEl => Element%prev
      else
         if (.not. associated(Oper(i)%p%fprev)) then
            Oper(i)%p%fprev => this%flast(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p
         else if (associated(this%flast(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p)) then
            if (Oper(i)%p%fprev%tau < this%flast(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p%tau)&
               Oper(i)%p%fprev => this%flast(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p
         end if
         mEl => this%last
      end if

      do
         if (.not. associated(mEl)) exit
         if (mEl%tau <= Oper(i)%p%tau) then
            if (associated(mEl) .and. Oper(i)%p%tau == mEl%tau) then
               Oper(i)%p%tau = Oper(i)%p%tau + this%equal_time_offset
            end if
            exit
         else
            mEl => mEl%prev
         end if
      end do

      if (.not. associated(Oper(i)%p%prev)) then
         Oper(i)%p%prev => mEl
      else if (associated(mEl)) then
         if (Oper(i)%p%prev%tau < mEl%tau)&
            Oper(i)%p%prev => mEl
      end if
      if (.not. associated(Oper(i)%p%next)) then
         if (associated(mEl)) then
            Oper(i)%p%next => mEl%next
         else
            Oper(i)%p%next => this%first
         end if
      else
         if (associated(mEl)) then
            if (associated(mEl%next)) then
               if (Oper(i)%p%next%tau >= mEl%next%tau)&
                  Oper(i)%p%next => mEl%next
            end if
         else
            if (associated(this%first)) then
               if (Oper(i)%p%next%tau >= this%first%tau)&
                  Oper(i)%p%next => this%first
            end if
         end if
      end if

      if (i < size(Oper)) then
         if (Oper(i)%p%Orbital /= Oper(i+1)%p%Orbital&
             .or. Oper(i)%p%Spin /= Oper(i+1)%p%Spin) then
            ! FIXME: If more than two operators are inserted and there
            ! are two operators of the same flavour that are separated in
            ! tau by one of different flavour, the count is not the same
            ! as with the original method. Check if this needs to be
            ! corrected.
            TFPos = 0
            Element => this%ffirst(Oper(i+1)%p%Orbital, Oper(i+1)%p%Spin)%p
         end if
      end if
      i = i + 1
   enddo

end subroutine process_OperAdd_global

!===============================================================================
!> Insert one or more properly prepared operators from the time-ordered Oper
!  array into the trace and adjust pointers to window edge operators if
!  necessary.
!
!  Proper preparation means that the (f)prev/(f)next pointers of all
!  operators point to the correct (flavour-)neighbours (by tau) in the
!  list after insertion. Operator arrays preprocessed by
!  process_OperAdd or removed (as a whole) from a correctly built
!  trace whose state is not different from the one right after the removal of
!  that operator array are properly prepared.
subroutine insert_Oper(this,Oper)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TOperPointer)                  :: Oper(:)
!local
   integer                             :: i

   ! change next and prev pointers to integrate Oper(i)s into list and set
   ! first, last, ffirst, and flast if necessary
   if(.not.associated(Oper(1)%p%prev))this%first=>Oper(1)%p
   if(.not.associated(Oper(size(Oper))%p%next))this%last=>Oper(size(Oper))%p
   do i=1,size(Oper)
      if (.not. associated(Oper(i)%p%fprev)) then
         this%ffirst(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p => Oper(i)%p
      else
         Oper(i)%p%fprev%fnext => Oper(i)%p
      end if
      if (.not. associated(Oper(i)%p%fnext)) then
         this%flast(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p => Oper(i)%p
      else
         Oper(i)%p%fnext%fprev => Oper(i)%p
      end if
      if (associated(Oper(i)%p%prev)) Oper(i)%p%prev%next => Oper(i)%p
      if (associated(Oper(i)%p%next)) Oper(i)%p%next%prev => Oper(i)%p
   enddo

   ! set (f)prewin and postwin pointers if Oper(i)s are inserted right
   ! before or after the window (this is only necessary for
   ! unrestricted moves)
   do i = 1, size(Oper)
      if (Oper(i)%p%tau < this%tauwin_min) then
         if (associated(this%fprewin(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p)) then
            if (Oper(i)%p%tau >= this%fprewin(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p%tau)&
               this%fprewin(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p => Oper(i)%p
         else
            this%fprewin(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p => Oper(i)%p
         end if
         if (associated(this%prewin)) then
            if (Oper(i)%p%tau >= this%prewin%tau)&
               this%prewin => Oper(i)%p
         else
            this%prewin => Oper(i)%p
         end if
      end if
   end do
   do i = size(Oper), 1, -1
      if (Oper(i)%p%tau > this%tauwin_max) then
         if (associated(this%postwin)) then
            if (Oper(i)%p%tau <= this%postwin%tau)&
               this%postwin => Oper(i)%p
         else
            this%postwin => Oper(i)%p
         end if
      end if
   end do
   
   !NOper goes over all operators including worms
   this%NOper = this%NOper + size(Oper)
   
   !we only count NOSOper over hybridization operators
   do i = 1, size(Oper)
      if (Oper(i)%p%has_hyb) then
         this%NOSOper(Oper(i)%p%Orbital, Oper(i)%p%Spin) =&
            this%NOSOper(Oper(i)%p%Orbital, Oper(i)%p%Spin) + 1
         this%NOSCAOper(Oper(i)%p%Orbital, Oper(i)%p%Spin, Oper(i)%p%CA) =&
            this%NOSCAOper(Oper(i)%p%Orbital, Oper(i)%p%Spin, Oper(i)%p%CA) + 1
      end if
   end do

end subroutine insert_Oper

!===============================================================================
!> Picks a randomly chosen pair of an annihilator and creator of the
!  same flavor inside the window or the worm.
!
!  Returns false if no proposals for a removal request can be generated.
!  Pointers inserted into the Oper array will be ordered by the operators' tau.
!
!  logical to switch between remove worm or random operator pairs
!  FPOS only for hyb operators
logical function gen_OperRemove(this,Oper,FPos,N,hybpairremfactor,with_hyb)
!===============================================================================
   type(TTrace)                        :: this
!input
   integer                             :: N
   logical                             :: with_hyb
!output
   integer                             :: FPos(2)
   type(TOperPointer)                  :: Oper(N)
   real(KINDR), intent(out)            :: hybpairremfactor
!local
   integer                             :: i,j
   integer                             :: NPair
   type(TOper),pointer                 :: Element, mEl

   gen_OperRemove=.true.

   if (with_hyb) then
      call ensure_valid_NPairs(this)
      if (this%NPairs <= 0) then
         gen_OperRemove = .false.
         return
      end if

      call ensure_valid_PreWinFCount(this)

      NPair = randint(1, this%NPairs)

      if (associated(this%prewin)) then
         Element => this%prewin%next
      else
         Element => this%first
      end if

      ! set pointers to the operators in the selected pair
      first_operator: do
         if (.not. associated(Element)) stop "npairs count logic error detected gen_operremove(2)"
         if (Element%has_hyb) then
            if (this%b_offdiag) then
               mEl => Element%next
            else
               mEl => Element%fnext
            end if
            do
               if (.not. associated(mEl)) exit
               if (mEl%tau > this%tauwin_max) exit
               if (mEl%tau > Element%tau + this%taudiff_max) exit

               if (mEl%has_hyb .and. mEl%CA /= Element%CA) then

                  NPair = NPair - 1
                  if (NPair <= 0) exit first_operator

               end if
               if (this%b_offdiag) then
                  mEl => mEl%next
               else
                  mEl => mEl%fnext
               end if
            end do
         end if

         Element => Element%next
      end do first_operator

      Oper(1)%p => Element
      Oper(2)%p => mEl

      ! calculate correct fpos for selected operators
      i = 1
      do
         Element => Element%fprev
         if (.not. associated(Element)) exit
         if (Element%tau < this%tauwin_min) exit
         if (Element%CA == Oper(1)%p%CA .and. Element%has_hyb) i = i + 1
      end do
      FPos(Oper(1)%p%CA) = this%PreWinFCount(Oper(1)%p%Orbital, Oper(1)%p%Spin, Oper(1)%p%CA) + i

      i = 1
      do
         mEl => mEl%fprev
         if (.not. associated(mEl)) exit
         if (mEl%tau < this%tauwin_min) exit
         if (mEl%CA == Oper(2)%p%CA .and. mEl%has_hyb) i = i + 1
      end do
      FPos(Oper(2)%p%CA) = this%PreWinFCount(Oper(2)%p%Orbital, Oper(2)%p%Spin, Oper(2)%p%CA) + i

      ! calculate insertion proposal probability factor due to the taus
      if (Oper(1)%p%CA == 2) then
         mEl => Oper(1)%p
      else
         mEl => Oper(2)%p
      end if

      hybpairremfactor = real(this%NPairs, KINDR)&
                         / ((this%tauwin_max - this%tauwin_min)&
                            * (min(mEl%tau + this%taudiff_max, this%tauwin_max)&
                               - max(mEl%tau - this%taudiff_max, this%tauwin_min)))

      return
   else  ! with_hyb == .false.
      if (.not. allocated(this%wormContainer)) then
         gen_OperRemove = .false.
         return
      else
         if (size(this%wormContainer) < N) then
            gen_OperRemove = .false.
            return
         end if

         do i = 1, N
            Oper(i)%p => this%wormContainer(i)%p
            do j = i - 1, 1, -1
               if (Oper(j)%p%tau > this%wormContainer(i)%p%tau) then
                  Oper(j+1)%p => Oper(j)%p
                  Oper(j)%p => this%wormContainer(i)%p
               else
                  exit
               end if
            end do
         end do

         return
      end if
   end if

end function gen_OperRemove

!===============================================================================
!> Picks two randomly chosen operators or the worm and checks if 
!!the configuration without those operators is compatible with the good quantum
!numbers. 
!!Returns true if it is compatible, false otherwise.
!!logical to switch between remove worm or random operator pairs
!! FPOS only for hyb operators
logical function gen_OperRemove_global(this,DStates,Oper,FPos,N,with_hyb)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
   integer                             :: N,Ntmp
   logical                             :: with_hyb
!output
   integer                             :: FPos(N)
   type(TOperPointer)                  :: Oper(N)
!local
   integer                             :: i,j,Pos(N),iPos,TFPos(DStates%NBands,2,2)
   type(TOper),pointer                 :: Element
!worm positions
   integer,allocatable                 :: wpos(:)

   gen_OperRemove_global=.true.
   Ntmp=this%NOper
   
   if(allocated(this%wormContainer)) then
      allocate(wpos(size(this%wormContainer)))
   endif

   !find worms
   if(allocated(this%wormContainer)) then
      if(with_hyb) then
         !we subtract the worms from NOper if we want to remove operators with
         !hyb
         Ntmp=this%NOper-size(this%wormContainer)
      endif
      !find position of all worms
      j=1
      Element=>this%first
      do i=1,this%NOper
         if(.not.Element%has_hyb) then
            wpos(j)=i
            j=j+1
         endif
         Element=>Element%next
      end do
   endif
   
      
   !in case we want to remove more operators than present
   if(Ntmp.lt.N)then
      gen_OperRemove_global=.false.
      if(allocated(wpos)) deallocate(wpos)
      return
   endif
   
   Pos=0
   do i=1,N
      !TODO: select creation/annihilation operator pair to avoid trivial rejects
      !making sure randomly chosen operators are not worms
      iPos=randint(1, Ntmp)
      !correcting for worm position
      if(with_hyb) then
         if(allocated(this%wormContainer)) then
               do j=1,size(this%wormContainer)
                  if(iPos>=wpos(j)) iPos=iPos+1
               enddo
         endif
      else
         iPos=wpos(i)
      endif

! time order the elements to be removed
      do j=i-1,1,-1
         if(iPos.eq.Pos(j))then
            gen_OperRemove_global=.false.
            return
         endif
         if(iPos.gt.Pos(j))exit
         Pos(j+1)=Pos(j)
      enddo
      Pos(j+1)=iPos
   enddo
   
   if(allocated(wpos)) deallocate(wpos)

   j=1
   TFPos=0
   Element=>this%first
   do i=1,this%NOper
      if(Element%has_hyb) then
         TFPos(Element%Orbital,Element%Spin,Element%CA)=TFPos(Element%Orbital,Element%Spin,Element%CA)+1
      endif
      if(Pos(j).eq.i)then
         FPos(Element%CA)=TFPos(Element%Orbital,Element%Spin,Element%CA)
         Oper(j)%p=>Element
         if(j.lt.N)j=j+1
      endif
      Element=>Element%next
   enddo 

end function gen_OperRemove_global


!===============================================================================
!> Remove the operators specified in the time-ordered Oper array parameter from
!  the trace and adjust pointers to window edge operators if necessary .
subroutine remove_Oper(this,Oper)
!===============================================================================
   type(TTrace)                        :: this
   type(TOperPointer)                  :: Oper(:)
!local
   type(TOper),pointer                 :: Element
   integer                             :: i,j,iB,iS

   ! set (f)prewin and postwin pointers if Oper(i)s are removed from
   ! right before or after the window (this is only necessary for
   ! unrestricted moves)
   do i = size(Oper), 1, -1
      if (associated(this%fprewin(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p, Oper(i)%p))&
         this%fprewin(Oper(i)%p%Orbital, Oper(i)%p%Spin)%p => Oper(i)%p%fprev
      if (associated(this%prewin, Oper(i)%p))&
         this%prewin => Oper(i)%p%prev
   end do
   do i = 1, size(Oper)
      if (associated(this%postwin, Oper(i)%p))&
         this%postwin => Oper(i)%p%next
   end do

!fixing the start of the trace
   do i=1,size(Oper)
      if(associated(this%first,Oper(i)%p))then
         this%first=>Oper(i)%p%next
      else
         exit
      endif
   enddo
!fixing the "next" elements in the trace
   if(associated(this%first).and.i.le.size(Oper))then
      Element=>Oper(i)%p%prev
      do j=i,size(Oper)
         if(.not.associated(Element%next,Oper(j)%p))&
            Element=>Oper(j)%p%prev
         Element%next=>Oper(j)%p%next
      enddo   
   endif

   ! correct references to operators being removed at flavour-list beginning
   ! and flavour-neighbours
   do iB = 1, size(this%ffirst, 1)
      do iS = 1, size(this%ffirst, 2)
         do i = 1, size(Oper)
            if (Oper(i)%p%Orbital /= iB .or. Oper(i)%p%Spin /= iS) cycle
            if (associated(this%ffirst(iB, iS)%p, Oper(i)%p)) then
               this%ffirst(iB, iS)%p => Oper(i)%p%fnext
            else
               exit
            end if
         end do
         if (i > size(Oper)) cycle
         Element => Oper(i)%p%fprev
         do j = i, size(Oper)
            if (Oper(j)%p%Orbital /= iB .or. Oper(j)%p%Spin /= iS) cycle
            if (.not. associated(Element)) cycle
            if (.not. associated(Element%fnext, Oper(j)%p))&
               Element => Oper(j)%p%fprev
            Element%fnext => Oper(j)%p%fnext
         end do
      end do
   end do

!fixing the end of the trace
   do i=size(Oper),1,-1
      if(associated(this%last,Oper(i)%p))then
         this%last=>Oper(i)%p%prev
      else
         exit
      endif
   enddo
!fixing the "prev" elements in the trace
   if(associated(this%last).and.i.ge.1)then
      Element=>Oper(i)%p%next
      do j=i,1,-1
         if(.not.associated(Element%prev,Oper(j)%p))&
            Element=>Oper(j)%p%next
         Element%prev=>Oper(j)%p%prev
      enddo   
   endif

   ! correct references to operators being removed at flavour-list end
   ! and flavour-neighbours
   do iB = 1, size(this%ffirst, 1)
      do iS = 1, size(this%ffirst, 2)
         do i = size(Oper), 1, -1
            if (Oper(i)%p%Orbital /= iB .or. Oper(i)%p%Spin /= iS) cycle
            if (associated(this%flast(iB, iS)%p, Oper(i)%p)) then
               this%flast(iB, iS)%p => Oper(i)%p%fprev
            else
               exit
            end if
         end do
         if (i < 1) cycle
         Element => Oper(i)%p%fnext
         do j = i, 1, -1
            if (Oper(j)%p%Orbital /= iB .or. Oper(j)%p%Spin /= iS) cycle
            if (.not. associated(Element)) cycle
            if (.not. associated(Element%fprev, Oper(j)%p))&
               Element => Oper(j)%p%fnext
            Element%fprev => Oper(j)%p%fprev
         end do
      end do
   end do

   !setting the right values for the number of operators in the trace
   !NOper goes over all operators including worms
   this%NOper = this%NOper - size(Oper)
   
   !TODO: only counting hybs to NOSOper, is this smart?
   do i = 1, size(Oper)
      if (Oper(i)%p%has_hyb) then
         this%NOSOper(Oper(i)%p%Orbital, Oper(i)%p%Spin) =&
            this%NOSOper(Oper(i)%p%Orbital, Oper(i)%p%Spin) - 1
         this%NOSCAOper(Oper(i)%p%Orbital, Oper(i)%p%Spin, Oper(i)%p%CA) =&
            this%NOSCAOper(Oper(i)%p%Orbital, Oper(i)%p%Spin, Oper(i)%p%CA) - 1
      end if
   end do

end subroutine remove_Oper

!===============================================================================
!> clears the entire trace by selecting arbitrary operator pairs consecutively
!===============================================================================
subroutine clear_Trace(this,DStates)
   type(TTrace)                        :: this
   type(TStates)                       :: DStates
!local
   type(TOperPointer),allocatable      :: Oper(:)
   type(TOper),pointer                 :: Element
   integer                             :: i,N,CA,iB,iS
   
   N=this%NOper
   
   allocate(Oper(N))
   Element=>this%first
   
   do i=1,this%NOper
      Oper(i)%p=>Element
      Element=>Element%next
   enddo
   
   do i=1,N,2
      call remove_Oper(this,Oper(i:i+1))
   enddo

   this%prewin => null()
   this%postwin => null()
   do iB = 1, DStates%NBands
      do iS = 1, 2
         this%fprewin(iB, iS)%p => null()
      enddo
   enddo
   this%NPairs = 0
   this%window_position = 0
   this%last_prewin => null()
   call update_PreWinFCount(this)
   call set_window(this)

   !updating the empty trace
   this%Trace = get_Trace_EB(this, DStates, global=.true.)
   call update_trace_EB(this, global=.true.)
         
   do CA=1,size(Oper)
      if (associated(Oper(CA)%p)) then
         if(allocated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(allocated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(allocated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(allocated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         deallocate(Oper(CA)%p)
      end if
   enddo
   
   if(associated(this%MInv))then
      do iB=1,size(this%MInv(:,1))
      do iS=1,2
         deallocate(this%MInv(iB,iS)%Mat)
      enddo
      enddo
      deallocate(this%Minv)
   endif
   
   allocate(this%MInv(DStates%NBands,2))
   if(associated(this%MInv))then
      do iB=1,size(this%MInv(:,1))
      do iS=1,2
         allocate(this%MInv(iB,iS)%Mat(0,0))
      enddo
      enddo
   endif
   
   if(associated(this%MInv_full))then
      if (associated(this%MInv_full%Mat)) deallocate(this%MInv_full%Mat)
      deallocate(this%Minv_full)
   endif

   allocate(this%MInv_full)
   allocate(this%MInv_full%Mat(0,0))


   !dealloctaing worm container
   if(allocated(this%wormContainer)) then
      deallocate(this%wormContainer)
   endif
      
   deallocate(Oper)
   
   !we have to reset the det manually
   this%Det = TLogDet(log = 0.0_KINDR, sign = 1.0_KINDR)

   this%cfgsign = -get_Sign(this, DStates)
               
end subroutine clear_Trace

!===============================================================================
!> Calculates the value of the bosonic trace.
real(KINDR) function get_BosonicTrace(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   integer                             :: iB
   type(TOper),pointer                 :: Element1,Element2
   real(KINDR)                         :: BWeight(DStates%NBands),tau

   if(this%Phonon.eq.1)then
      BWeight=0d0
      Element1=>this%first
      do while(associated(Element1))
         Element2=>Element1%next
         do while(associated(Element2))
!            if(Element2%Orbital.eq.Element1%Orbital)then
               tau=Element2%tau-Element1%tau
               BWeight(Element1%Orbital)=BWeight(Element1%Orbital)+&
                   (-1d0)**(Element1%CA/2)*(-1d0)**(Element2%CA/2)*&
                   get_lin_screening_function(this,tau,Element1%Orbital,Element1%Spin,Element2%Orbital,Element2%Spin)
!            endif
            Element2=>Element2%next
         enddo
         Element1=>Element1%next
      enddo
      get_BosonicTrace=1d0
      do iB=1,DStates%NBands
         BWeight(iB)=exp(BWeight(iB))
         get_BosonicTrace=get_BosonicTrace*BWeight(iB)
      enddo
   else
      get_BosonicTrace=1d0
   endif

end function get_BosonicTrace

pure subroutine split_log(vec, logfactor)
   real(KINDR), intent(inout) :: vec(:), logfactor
   real(KINDR)                :: mv

   mv = maxval(abs(vec))
   if (mv > 0.0_KINDR) vec = vec / mv
   logfactor = logfactor + log(mv)
end subroutine split_log

!===============================================================================
!> Calculates the value of the local trace for the current
!  configuration and stores and caches states and superstates in the
!  operators.
type(TLogTr) function get_Trace_EB(this,DStates,FixEndR,global)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
   type(TOper),pointer,optional        :: FixEndR
   logical, optional                   :: global
!local
   type(TOper),pointer                 :: ElementR, ElementL, StartR, StartL, EndR
   integer                             :: sst, nsst, NStates_sst, NStates_nsst
   integer                             :: iNTruncSt, louter_sst_size, outer_index
   integer                             :: vec_outer_ind
   real(KINDR)                         :: bra_t(DStates%NStatesMax), bra_logmax
! bra after time evolution
   real(KINDR)                         :: bra_tp(DStates%NStatesMax)
! ket before time evolution
   real(KINDR)                         :: ket_t(DStates%NStatesMax), ket_logmax
! ket after time evolution
   real(KINDR)                         :: ket_tp(DStates%NStatesMax)
   real(KINDR)                         :: tauR, tauL, braket, braket_logmax
   logical                             :: CalcNew, FullR, FullL, tainted

   get_Trace_EB = TLogTr(log = -huge(0.0_KINDR), sign = 0.0_KINDR)
   this%tparttrace(:) = TLogTr(log = -huge(0.0_KINDR), sign = 0.0_KINDR)

   if (this%sst_to_statesindex(this%outer_sst) == -1) then
      stop "unreachable: invalid outer sst"
   end if

   if (present(global)) then
      if (global) then
         CalcNew = .true.
      else
         CalcNew = .false.
      end if
   else
      CalcNew = .false.
   end if
   CalcNew = CalcNew .or. this%outer_sst /= this%outer_sst_old&
             .or. this%outer_state /= this%outer_state_old
   outer_index = this%sst_to_statesindex(this%outer_sst)
   if (this%b_statesampling) then
      louter_sst_size = 1
   else
      louter_sst_size = this%States(outer_index, 2)
   end if

   ! The trace will be recalculated between StartR and StartL or from
   ! tau = 0 instead of StartR if FullR is set and from tau = beta
   ! instead of StartL if FullL is set.
   !
   ! If the move is not a global move, it is never necessary to
   ! recalculate anything outside of the window, so operators close to
   ! the window edge are more suitable as starting points for the
   ! StartR/StartL search than the ends. Note that currently the
   ! operators right next to the window edge might be marked calc_new
   ! due to proposed removal of an operator inside of the window, but
   ! a recalculation of the state at that operator would still not be
   ! necessary. If no ket is stored before the current window edge,
   ! also check the operator before the last window edge.
   FullR = .true.
   if (.not. CalcNew) then
      if (associated(this%prewin)) then
         if (this%prewin%ket_state) then
            StartR => this%prewin
            FullR = .false.
         else if (associated(this%last_prewin)) then
            if (this%last_prewin%ket_state) then
               StartR => this%last_prewin
               FullR = .false.
            else
               StartR => this%first
            end if
         else
            StartR => this%first
         end if
      else
         StartR => this%first
      end if
   else
      StartR => this%first
   end if
   EndR => null()

   FullL = .true.
   if (.not. CalcNew) then
      if (associated(this%postwin)) then
         if (.not. this%postwin%ket_state) then
            StartL => this%postwin
            FullL = .false.
         else
            StartL => this%last
         end if
      else
         StartL => this%last
      end if
   else
      StartL => this%last
   end if

   ! try to move StartR/StartL further inward to avoid as many
   ! calculations as possible; note that StartL will always be
   ! associated if StartR is
   if (associated(StartR)) then
      ! FullR/FullL must be unset if it is possible to continue the
      ! calculation from the state stored in the StartR/StartL
      ! operator
      if (.not. CalcNew .and. .not. StartR%calc_new .and. StartR%ket_state) then
         FullR = .false.
         do
            if (.not. associated(StartR%next)) exit
            if (StartR%next%calc_new .or. .not. StartR%next%ket_state) exit
            StartR => StartR%next
         end do
      end if

      if (.not. CalcNew .and. .not. StartL%calc_new .and. .not. StartL%ket_state) then
         FullL = .false.
         do
            if (.not. associated(StartL%prev)) exit
            if (StartL%prev%calc_new .or. StartL%prev%ket_state) exit
            StartL => StartL%prev
         end do
      end if

      ! no operators from window edge to opposite end of trace; have
      ! StartR/L point to different operators
      if (associated(StartR, StartL)) then
         if (associated(StartL%next)) then
            StartL => StartL%next
         else
            StartR => StartR%prev
         end if
      end if

      ! The calculation is started after StartR first and continues to
      ! calculate kets until EndR is encountered. EndR is supposed to
      ! balance the amounts of calculated ket and bra states inside
      ! the window by calculating bras after the last changed operator
      ! as this reduces the need for recalculation after the following
      ! move. This is a significant performance improvement and it
      ! might be possible to increase it even further by calculating
      ! from both sides ~"to the middle of the window", but this would
      ! require more extra code and since the window is always shifted
      ! to higher taus I chose a bias in favor of kets (which are
      ! calculated not only to the first changed operator but also
      ! between the changed operators).
      !
      ! For global moves this does not work, so there is an extra
      ! condition to stop calculation of kets when the first operator
      ! after the window is encountered.
      if (.not. present(FixEndR)) then
         EndR => StartL
         if (.not. CalcNew .and. .not. EndR%calc_new) then
            do
               if (.not. associated(EndR%prev)) exit
               if (EndR%prev%calc_new) exit
               EndR => EndR%prev
            end do
         end if

         if (associated(this%postwin) .and. associated(EndR)) then
            if (EndR%tau > this%postwin%tau) EndR => this%postwin
         end if
      else
         EndR => FixEndR
         if (.not. associated(EndR)) then
            StartL => this%last
            FullL = .true.
         else
            if (StartR%tau >= EndR%tau) then
               StartR => EndR%prev
               if (.not. associated(StartR)) then
                  StartR => this%first
                  FullR = .true.
               end if
            end if
            if (StartL%tau < EndR%tau) then
               StartL => EndR
            end if
         end if
      end if

      StartR%cache_written = .false.
      StartL%cache_written = .false.
   end if

   ! store for update_trace_eb
   if (.not. FullR) then
      this%UpdateR => StartR
   else
      this%UpdateR => null()
   end if
   if (.not. FullL) then
      this%UpdateL => StartL
   else
      this%UpdateL => null()
   end if
   this%LastEndR => EndR

   ! call check_preconditions()

   StateLoop: do iNTruncSt = this%outer_state, this%outer_state + louter_sst_size - 1
      if (this%b_statesampling) then
         vec_outer_ind = 1
      else
         vec_outer_ind = iNTruncSt
      end if
!=====================================================
! start from the suitable starting point at lesser tau
!=====================================================
      ElementR => StartR

      if (FullR) then
         sst = this%outer_sst
         NStates_sst = this%States(outer_index, 2)
         ket_logmax = 0.0_KINDR
         ket_t = 0.0_KINDR
         ket_t(iNTruncSt) = 1.0_KINDR
         tauR = 0.0_KINDR
         tainted = .true.
      else
         sst = ElementR%following_sst
         NStates_sst = DStates%SubStates(sst)%NStates
         ket_logmax = ElementR%slogmax(vec_outer_ind)
         ket_t(1:NStates_sst) = ElementR%state(1:NStates_sst, vec_outer_ind)
         tauR = ElementR%tau
         tainted = ElementR%calc_new
         if (vec_outer_ind == louter_sst_size) ElementR%calc_new = .false.
         ElementR => ElementR%next
      end if

      ElementRLoop: do while (associated(ElementR))
         if (associated(ElementR, EndR)) exit
         if (ElementR%calc_new) tainted = .true.

         nsst = DStates%SubStates(sst)%Connect(ElementR%Orbital,ElementR%Spin,ElementR%CA)
         NStates_nsst = DStates%SubStates(nsst)%NStates !rows

         if (NStates_sst > 1) then
! exp(-tau H)ket
            ket_tp(1:NStates_sst) =&
               exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                     - DStates%substates(sst)%Eval(0)) * (ElementR%tau - tauR))&
               * ket_t(1:NStates_sst)
            ket_logmax = ket_logmax&
               -(DStates%substates(sst)%Eval(0) - this%Egs) * (ElementR%tau - tauR)
            call split_log(ket_tp(1:NStates_sst), ket_logmax)

! O(^+)ket
!ket_t(0:DStates%SubStates(nsst)%NStates-1)=&
!&matmul(DStates%SubStates(sst)%Psis_EB(ElementR%Orbital,ElementR%Spin,ElementR%CA)%Op,&
!        ket_tp(0:DStates%SubStates(sst)%NStates-1))

            call DGEMV('N',NStates_nsst,NStates_sst,1d0,&
                       DStates%SubStates(sst)%Psis_EB(ElementR%Orbital,ElementR%Spin,ElementR%CA)%Op,&
                       NStates_nsst,&
                       ket_tp(1:NStates_sst),1,0d0,&
                       ket_t(1:NStates_nsst),1)
            call split_log(ket_t(1:NStates_nsst), ket_logmax)
         else
            ket_logmax = ket_logmax&
               -(DStates%substates(sst)%Eval(0) - this%Egs) * (ElementR%tau - tauR)
            ket_t(1:NStates_nsst)&
               = ket_t(1)&
                 * DStates%SubStates(sst)%Psis_EB(ElementR%Orbital,ElementR%Spin,ElementR%CA)%Op(:, 1)
         end if

         ! once a changed operator is encountered, the calculated
         ! states are stored in cache and not in state
         if (.not. tainted) then
            if (vec_outer_ind == 1) then
               if (this%b_fix_prealloc_stvec) then
                  ElementR%state(:, :) = 0.0_KINDR
                  ElementR%slogmax(:) = -huge(ElementR%slogmax(1))
               else
               if (allocated(ElementR%state)) then
                  deallocate(ElementR%state)
                  deallocate(ElementR%slogmax)
               end if
               allocate(ElementR%state(NStates_nsst, louter_sst_size))
               allocate(ElementR%slogmax(louter_sst_size))
               end if

               ElementR%preceding_sst = sst
               ElementR%following_sst = nsst
               ElementR%ket_state = .true.
               ElementR%cache_written = .false.
               this%UpdateR => ElementR
            end if
            if (vec_outer_ind == louter_sst_size) then
               ElementR%calc_new = .false.
            end if

            ElementR%state(1:NStates_nsst, vec_outer_ind) = ket_t(1:NStates_nsst)
            ElementR%slogmax(vec_outer_ind) = ket_logmax
         else
            if (vec_outer_ind == 1) then
               if (this%b_fix_prealloc_stvec) then
                  ElementR%cache(:, :) = 0.0_KINDR
                  ElementR%clogmax(:) = -huge(ElementR%clogmax(1))
               else
               if (allocated(ElementR%cache)) then
                  deallocate(ElementR%cache)
                  deallocate(ElementR%clogmax)
               end if
               allocate(ElementR%cache(NStates_nsst, louter_sst_size))
               allocate(ElementR%clogmax(louter_sst_size))
               end if

               ElementR%cache_pre_sst = sst
               ElementR%cache_post_sst = nsst
               ElementR%ket_cache = .true.
               ElementR%cache_written = .true.
            end if
            if (vec_outer_ind == louter_sst_size) then
               ElementR%calc_new = .false.
            end if

            ElementR%cache(1:NStates_nsst, vec_outer_ind) = ket_t(1:NStates_nsst)
            ElementR%clogmax(vec_outer_ind) = ket_logmax
         end if

         tauR = ElementR%tau
         sst = nsst
         NStates_sst = NStates_nsst
         ElementR => ElementR%next
      enddo ElementRLoop

      ! to use as stopping condition in L loop
      if (associated(ElementR)) then
         ElementR => ElementR%prev
      else
         ElementR => this%last
      end if

!======================================================
! start from the suitable starting point at greater tau
!======================================================
      ElementL => StartL

      if (FullL) then
         sst = this%outer_sst
         NStates_sst = this%States(outer_index, 2)
         bra_logmax = 0.0_KINDR
         bra_t = 0.0_KINDR
         bra_t(iNTruncSt) = 1.0_KINDR
         tauL = this%beta
         tainted = .true.
      else
         sst = ElementL%preceding_sst
         NStates_sst = DStates%SubStates(sst)%NStates
         bra_logmax = ElementL%slogmax(vec_outer_ind)
         bra_t(1:NStates_sst) = ElementL%state(1:NStates_sst, vec_outer_ind)
         tauL = ElementL%tau
         tainted = ElementL%calc_new
         if (vec_outer_ind == louter_sst_size) ElementL%calc_new = .false.
         ElementL => ElementL%prev
      end if

      ElementLLoop: do while (associated(ElementL))
         if (associated(ElementL, ElementR)) exit
         if (ElementL%calc_new) tainted = .true.

         nsst = DStates%SubStates(sst)%Connect(ElementL%Orbital,ElementL%Spin,3-ElementL%CA)
         NStates_nsst = DStates%SubStates(nsst)%NStates !rows

         if (NStates_sst > 1) then
! exp(-tau H)bra
            bra_tp(1:NStates_sst) =&
               exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                     - DStates%substates(sst)%Eval(0)) * (-(ElementL%tau - tauL)))&
               * bra_t(1:NStates_sst)
            bra_logmax = bra_logmax&
               -(DStates%substates(sst)%Eval(0) - this%Egs) * (-(ElementL%tau - tauL))
            call split_log(bra_tp(1:NStates_sst), bra_logmax)

! O(^+)bra
!bra_t(0:DStates%SubStates(nsst)%NStates-1)=&
!matmul(DStates%SubStates(sst)%Psis_EB(ElementL%Orbital,ElementL%Spin,3-ElementL%CA)%Op,&
!       bra_tp(0:DStates%SubStates(sst)%NStates-1))

            call DGEMV('N',NStates_nsst,NStates_sst,1d0,&
                       DStates%SubStates(sst)%Psis_EB(ElementL%Orbital,ElementL%Spin,3-ElementL%CA)%Op,&
                       NStates_nsst,&
                       bra_tp(1:NStates_sst),1,0d0,&
                       bra_t(1:NStates_nsst),1)
            call split_log(bra_t(1:NStates_nsst), bra_logmax)
         else
            bra_logmax = bra_logmax&
               -(DStates%substates(sst)%Eval(0) - this%Egs) * (-(ElementL%tau - tauL))
            bra_t(1:NStates_nsst)&
               = bra_t(1)&
                 * DStates%SubStates(sst)%Psis_EB(ElementL%Orbital,ElementL%Spin,3-ElementL%CA)%Op(:, 1)
         end if

         if (.not. tainted) then
            if (vec_outer_ind == 1) then
               if (this%b_fix_prealloc_stvec) then
                  ElementL%state(:, :) = 0.0_KINDR
                  ElementL%slogmax(:) = -huge(ElementL%slogmax(1))
               else
               if (allocated(ElementL%state)) then
                  deallocate(ElementL%state)
                  deallocate(ElementL%slogmax)
               end if
               allocate(ElementL%state(NStates_nsst, louter_sst_size))
               allocate(ElementL%slogmax(louter_sst_size))
               end if

               ElementL%preceding_sst = nsst
               ElementL%following_sst = sst
               ElementL%ket_state = .false.
               ElementL%cache_written = .false.
               this%UpdateL => ElementL
            end if
            if (vec_outer_ind == louter_sst_size) then
               ElementL%calc_new = .false.
            end if

            ElementL%state(1:NStates_nsst, vec_outer_ind) = bra_t(1:NStates_nsst)
            ElementL%slogmax(vec_outer_ind) = bra_logmax
         else
            if (vec_outer_ind == 1) then
               if (this%b_fix_prealloc_stvec) then
                  ElementL%cache(:, :) = 0.0_KINDR
                  ElementL%clogmax(:) = -huge(ElementL%clogmax(1))
               else
               if (allocated(ElementL%cache)) then
                  deallocate(ElementL%cache)
                  deallocate(ElementL%clogmax)
               end if
               allocate(ElementL%cache(NStates_nsst, louter_sst_size))
               allocate(ElementL%clogmax(louter_sst_size))
               end if

               ElementL%cache_post_sst = sst
               ElementL%cache_pre_sst = nsst
               ElementL%ket_cache = .false.
               ElementL%cache_written = .true.
            end if
            if (vec_outer_ind == louter_sst_size) then
               ElementL%calc_new = .false.
            end if

            ElementL%cache(1:NStates_nsst, vec_outer_ind) = bra_t(1:NStates_nsst)
            ElementL%clogmax(vec_outer_ind) = bra_logmax
         end if

         tauL = ElementL%tau
         sst = nsst
         NStates_sst = NStates_nsst
         ElementL => ElementL%prev
      enddo ElementLLoop

      bra_tp(1:NStates_sst) =&
         exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
               - DStates%substates(sst)%Eval(0)) * (tauL - tauR))&
         * bra_t(1:NStates_sst)
      bra_logmax = bra_logmax&
            -(DStates%substates(sst)%Eval(0) - this%Egs) * (tauL - tauR)
      call split_log(bra_tp(1:NStates_sst), bra_logmax)

      braket = dot_product(ket_t(1:NStates_sst), bra_tp(1:NStates_sst))
      braket_logmax = ket_logmax + bra_logmax + log(abs(braket))

      this%tparttrace(iNTruncSt) = TLogTr(log = braket_logmax,&
                                          sign = sign(1.0_KINDR, braket))
      get_Trace_EB = get_Trace_EB + this%tparttrace(iNTruncSt)

   enddo StateLoop

   ! if the move is not a global move, states up to the current window
   ! edge were calculated and stored in state, so we can start from
   ! there
   if (.not. CalcNew) this%last_prewin => this%prewin

   !call print_trace(this,9992,get_trace_eb)

! contains
!    subroutine check_preconditions()
!       type(TOper), pointer :: ptr
! !      logical              :: passedStartR, passedStartL

!       if (this%NOper == 0 .and. (associated(StartR) .or. associated(StartL)))&
!          stop "Start(R/L) associated to something with 0 operators in trace"
!       if (this%NOper > 0 .and. (.not. associated(StartL) .or. .not. associated(StartR)))&
!          stop "Start(R/L) incorrectly null"
!       if (this%NOper > 0 .and. CalcNew&
!           .and. (.not. associated(StartR, this%first)&
!                  .or. .not. associated(StartL, this%last)))&
!          stop "Start(R/L) not at edge for global move"

!       ptr => this%first
!       ! before StartR
!       do
!          if (.not. associated(ptr))&
!             stop "unreachable: StartR or StartL not encountered"
!          if (associated(ptr, StartL))&
!             stop "unreachable: StartL encountered before StartR"
!          if (FullR .and. associated(ptr, StartR)) exit
!          if (.not. ptr%ket_state .or. ptr%calc_new .or. CalcNew)&
!             stop "unreachable: unsuitable StartR"
!          if (associated(StartR)) exit
!          ptr => ptr%next
!       end do

!       ! between StartR and EndR
!       do
!          if (.not. associated(ptr))&
!             stop "unreachable: StartL not encountered"
!          if (associated(ptr, StartL) .and. .not. FullL)&
!             stop "unreachable: StartL encountered before EndR"
!          if (FullR .and. associated(ptr, StartR)) exit
!          if (.not. ptr%ket_state .or. ptr%calc_new .or. CalcNew)&
!             stop "unreachable: unsuitable StartR"
!          if (associated(StartR)) exit
!          ptr => ptr%next
!       end do


!       if (associated(StartR)) then
!          if (.not. StartR%ket_state&
!          .or. StartR%calc_new&
!          .or. .not. allocated(StartR%state)) then
!             if (allocated(StartR%prev)) then
!                if (.not. StartR%prev%ket_state&
!                .or. StartR%prev%calc_new&
!                .or. .not. allocated(StartR%prev%state))&
!                stop "StartR restart impossible"
!             end if
!          end if
!          if (StartR%tau >= StartL%tau) stop "incorrect ordering of StartR and StartL"
!          if (StartL%ket_state&
!          .or. StartL%calc_new&
!          .or. .not. allocated(StartL%state)) then
!             if (associated(StartL%next)) then
!                if (StartL%next%ket_state&
!                .or. StartL%next%calc_new&
!                .or. .not. allocated(StartL%next%state))&
!                stop "StartL restart impossible"
!             end if
!          end if
!       end if
!    end subroutine check_preconditions

end function get_Trace_EB



!===============================================================================
!> To be called after accepting a move. Updates the states and
!  superstates stored in the operators and marks beta/2 states as
!  invalidated if necessary.
subroutine update_trace_EB(this, global)
!===============================================================================
   type(TTrace)                        :: this
   logical, optional                   :: global
!local
   type(TOper),pointer                 :: ElementR, StartR, StartL, EndR
   real(KINDR), allocatable            :: temp(:, :), templm(:)
   logical                             :: CalcNew, FullR, FullL

   this%parttrace(:)=this%tparttrace(:)
   if (this%b_statesampling) then
      this%outer_sst_size = 1
   else
      this%outer_sst_size = this%States(this%sst_to_statesindex(this%outer_sst), 2)
   end if

   if (present(global)) then
      if (global) then
         CalcNew = .true.
      else
         CalcNew = .false.
      end if
   else
      CalcNew = .false.
   end if
   CalcNew = CalcNew .or. this%outer_sst /= this%outer_sst_old&
             .or. this%outer_state /= this%outer_state_old

   if (CalcNew .or. (this%tauwin_min <= this%beta/2.0 .and. this%tauwin_max >= this%beta/2.0)) then
      this%ket_b2_invalid = .true.
      this%bra_b2_invalid = .true.
   else if (this%tauwin_min >= this%beta/2.0) then
      this%bra_b2_invalid = .true.
   else
      this%ket_b2_invalid = .true.
   end if

   this%last_prewin => this%prewin

   ! start from pointer stored by get_trace_eb (it is critical that
   ! only elements that were visited during the last execution of
   ! get_trace_eb are revisited here because others may contain wrong
   ! cache information)
   if (.not. associated(this%UpdateR)) then
      FullR = .true.
      StartR => this%first
   else
      FullR = .false.
      StartR => this%UpdateR
   end if
   if (.not. associated(this%UpdateL)) then
      FullL = .true.
      StartL => this%last
   else
      FullL = .false.
      StartL => this%UpdateL
   end if
   EndR => this%LastEndR

   if (.not. FullR) then
      ElementR => StartR%next
   else
      ElementR => StartR
   end if

   if (FullL) then
      EndR => null()
   else
      EndR => StartL
   end if

   ElementLoop: do
      if (.not. associated(ElementR)) exit
      if (associated(ElementR, EndR)) exit

      if (.not. ElementR%cache_written) then
         stop "update_trace_eb assertion error: ElementR%cache_written == false"
      else
         ElementR%following_sst = ElementR%cache_post_sst
         ElementR%preceding_sst = ElementR%cache_pre_sst
         ElementR%ket_state = ElementR%ket_cache

         call move_alloc(ElementR%state, temp)
         call move_alloc(ElementR%cache, ElementR%state)
         call move_alloc(temp, ElementR%cache)

         call move_alloc(ElementR%slogmax, templm)
         call move_alloc(ElementR%clogmax, ElementR%slogmax)
         call move_alloc(templm, ElementR%clogmax)

         ElementR%cache_written = .false.
      end if

      ElementR => ElementR%next
   enddo ElementLoop

   if (associated(ElementR)) then
      if (ElementR%cache_written)&
         stop "update_trace_eb assertion error: ElementR%cache_written == true"
   end if
end subroutine update_trace_EB

!===============================================================================
!> Calculates and stores states at beta/2 if necessary
subroutine update_b2_states_eb(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   type(TOper),pointer                 :: Element, StartR, StartL
   integer                             :: sst, nsst, NStates_sst, NStates_nsst
   integer                             :: iNTruncSt, outer_index, vec_outer_ind
! bra before time evolution
   real(KINDR)                         :: bra_t(DStates%NStatesMax), bra_logmax
! bra after time evolution
   real(KINDR)                         :: bra_tp(DStates%NStatesMax)
! ket before time evolution
   real(KINDR)                         :: ket_t(DStates%NStatesMax), ket_logmax
! ket after time evolution
   real(KINDR)                         :: ket_tp(DStates%NStatesMax)
   real(KINDR)                         :: tauR, tauL, betahalf
   logical                             :: FullR, FullL

   if (.not. this%ket_b2_invalid .and. .not. this%bra_b2_invalid) return
   betahalf = this%beta/2.0
   outer_index = this%sst_to_statesindex(this%outer_sst)

   ! Determine the operator closest to beta/2 that has stored ket states.
   if (this%ket_b2_invalid) then
      FullR = .true.
      if (associated(this%prewin)) then
         if (this%prewin%tau <= betahalf .and. this%prewin%ket_state) then
            StartR => this%prewin
            FullR = .false.
         else if (associated(this%last_prewin)) then
            if (this%last_prewin%tau <= betahalf&
                .and. this%last_prewin%ket_state) then
               StartR => this%last_prewin
               FullR = .false.
            else
               StartR => this%first
            end if
         else
            StartR => this%first
         end if
      else
         StartR => this%first
      end if
   end if

   if (this%bra_b2_invalid) then
      FullL = .true.
      if (associated(this%postwin)) then
         if (this%postwin%tau > betahalf .and. .not. this%postwin%ket_state) then
            StartL => this%postwin
            FullL = .false.
         else
            StartL => this%last
         end if
      else
         StartL => this%last
      end if
   end if

   if (this%ket_b2_invalid .and. associated(StartR)) then
      if (StartR%ket_state .and. StartR%tau <= betahalf) then
         FullR = .false.
         do
            if (.not. associated(StartR%next)) exit
            if (StartR%next%tau > betahalf .or. .not. StartR%next%ket_state) exit
            StartR => StartR%next
         end do
      end if
   end if

   if (this%bra_b2_invalid .and. associated(StartL)) then
      if (.not. StartL%ket_state .and. StartL%tau > betahalf) then
         FullL = .false.
         do
            if (.not. associated(StartL%prev)) exit
            if (StartL%prev%tau <= betahalf .or. StartL%prev%ket_state) exit
            StartL => StartL%prev
         end do
      end if
   end if

   StateLoop: do iNTruncSt = this%outer_state, this%outer_state + this%outer_sst_size - 1
      if (this%b_statesampling) then
         vec_outer_ind = 1
      else
         vec_outer_ind = iNTruncSt
      end if
!=====================================================
! start from the suitable starting point at lesser tau
!=====================================================
      if (this%ket_b2_invalid) then
         if (vec_outer_ind == 1 .and. associated(this%ket_b2)) then
            deallocate(this%ket_b2)
            deallocate(this%ket_b2_logmax)
         end if

         Element => StartR

         if (FullR) then
            sst = this%outer_sst
            NStates_sst = this%States(outer_index, 2)
            ket_logmax = 0.0_KINDR
            ket_t = 0.0_KINDR
            ket_t(iNTruncSt) = 1.0_KINDR
            tauR = 0.0_KINDR
         else
            sst = Element%following_sst
            NStates_sst = DStates%SubStates(sst)%NStates
            ket_logmax = Element%slogmax(vec_outer_ind)
            ket_t(1:NStates_sst) = Element%state(1:NStates_sst, vec_outer_ind)
            tauR = Element%tau
            Element => Element%next
         end if

         ElementRLoop: do while(associated(Element))
            nsst = Element%following_sst
            NStates_nsst = DStates%SubStates(nsst)%NStates !rows

            !=============== store ket at beta/2 if it was passed ==============
            if (element%tau > betahalf .and. tauR <= betahalf) exit

! exp(-tau H)ket
            ket_tp(1:NStates_sst) =&
               exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                     - DStates%substates(sst)%Eval(0)) * (Element%tau - tauR))&
               *ket_t(1:NStates_sst)
            ket_logmax = ket_logmax&
               -(DStates%substates(sst)%Eval(0) - this%Egs) * (Element%tau - tauR)
            call split_log(ket_tp(1:NStates_sst), ket_logmax)

! O(^+)ket
!ket_t(0:DStates%SubStates(nsst)%NStates-1)=&
!matmul(DStates%SubStates(sst)%Psis_EB(ElementR%Orbital,ElementR%Spin,ElementR%CA)%Op,&
!       ket_tp(0:DStates%SubStates(sst)%NStates-1))

            call DGEMV('N',NStates_nsst,NStates_sst,1d0,&
                       DStates%SubStates(sst)%Psis_EB(Element%Orbital,Element%Spin,Element%CA)%Op,&
                       NStates_nsst,&
                       ket_tp(1:NStates_sst),1,0d0,&
                       ket_t(1:NStates_nsst),1)
            call split_log(ket_t(1:NStates_nsst), ket_logmax)

            tauR = Element%tau
            sst = nsst
            NStates_sst = NStates_nsst

            !========= store ket at beta/2 if it is not going to be passed ========
            if (.not. associated(Element%next)) exit

            Element => Element%next
         enddo ElementRLoop

         this%sst_ket_b2 = sst
         if (vec_outer_ind == 1) then
            allocate(this%ket_b2(NStates_sst, this%outer_sst_size))
            allocate(this%ket_b2_logmax(this%outer_sst_size))
         end if
         ! evolve state to beta/2
         ket_t(1:NStates_sst) =&
            exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                  - DStates%substates(sst)%Eval(0)) * (betahalf - tauR))&
            * ket_t(1:NStates_sst)
         call split_log(ket_t(1:NStates_sst), ket_logmax)

         this%ket_b2(:, vec_outer_ind) = ket_t(1:NStates_sst)
         this%ket_b2_logmax(vec_outer_ind) = ket_logmax&
            -(DStates%substates(sst)%Eval(0) - this%Egs) * (betahalf - tauR)
      end if

!======================================================
! start from the suitable starting point at greater tau
!======================================================
      if (this%bra_b2_invalid) then
         if (vec_outer_ind == 1 .and. associated(this%bra_b2)) then
            deallocate(this%bra_b2)
            deallocate(this%bra_b2_logmax)
         end if

         Element => StartL

         if (FullL) then
            sst = this%outer_sst
            NStates_sst = this%States(outer_index, 2)
            bra_logmax = 0.0_KINDR
            bra_t = 0.0_KINDR
            bra_t(iNTruncSt) = 1.0_KINDR
            tauL = this%beta
         else
            sst = Element%preceding_sst
            NStates_sst = DStates%SubStates(sst)%NStates
            bra_logmax = Element%slogmax(vec_outer_ind)
            bra_t(1:NStates_sst) = Element%state(1:NStates_sst, vec_outer_ind)
            tauL = Element%tau
            Element => Element%prev
         end if

         ElementLLoop: do while(associated(Element))
            nsst = Element%preceding_sst
            NStates_nsst = DStates%SubStates(nsst)%NStates !rows

            !=============== store bra at beta/2 if it was passed ==============
            if (element%tau <= betahalf .and. tauL >= betahalf) exit

! exp(-tau H)bra
            bra_tp(1:NStates_sst) =&
               exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                     - DStates%substates(sst)%Eval(0)) * (-(Element%tau - tauL)))&
               * bra_t(1:NStates_sst)
            bra_logmax = bra_logmax&
               -(DStates%substates(sst)%Eval(0) - this%Egs) * (-(Element%tau - tauL))
            call split_log(bra_tp(1:NStates_sst), bra_logmax)

! O(^+)bra
!bra_t(0:DStates%SubStates(nsst)%NStates-1)=&
!matmul(DStates%SubStates(sst)%Psis_EB(ElementL%Orbital,ElementL%Spin,3-ElementL%CA)%Op,&
!       bra_tp(0:DStates%SubStates(sst)%NStates-1))

            call DGEMV('N',NStates_nsst,NStates_sst,1d0,&
                       DStates%SubStates(sst)%Psis_EB(Element%Orbital,Element%Spin,3-Element%CA)%Op,&
                       NStates_nsst,&
                       bra_tp(1:NStates_sst),1,0d0,&
                       bra_t(1:NStates_nsst),1)
            call split_log(bra_t(1:NStates_nsst), bra_logmax)

            tauL = Element%tau
            sst = nsst
            NStates_sst = NStates_nsst

            !========= store bra at beta/2 if it is not going to be passed ========
            if (.not. associated(Element%prev)) exit

            Element => Element%prev
         enddo ElementLLoop

         this%sst_bra_b2 = sst
         if (vec_outer_ind == 1) then
            allocate(this%bra_b2(NStates_sst, this%outer_sst_size))
            allocate(this%bra_b2_logmax(this%outer_sst_size))
         end if
         ! evolve state to beta/2
         bra_t(1:NStates_sst) =&
            exp((DStates%substates(sst)%Eval(0:NStates_sst-1)&
                 - DStates%substates(sst)%Eval(0)) * (betahalf - tauL))&
            * bra_t(1:NStates_sst)
         call split_log(bra_t(1:NStates_sst), bra_logmax)

         this%bra_b2(:, vec_outer_ind) = bra_t(1:NStates_sst)
         this%bra_b2_logmax(vec_outer_ind) = bra_logmax&
            + (DStates%substates(sst)%Eval(0) - this%Egs) * (betahalf - tauL)
      end if
   enddo StateLoop

   this%ket_b2_invalid = .false.
   this%bra_b2_invalid = .false.
end subroutine update_b2_states_eb

!===============================================================================
subroutine meas_ntau_n0(this, DStates, ntau_n0)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!inout
   real(KINDR)                         :: ntau_n0(:, :, :, :, 0:)
!local
   type(TOper),pointer                 :: Element, KetEnd
   integer                             :: sst, nsst, NStates_sst, NStates_nsst
   integer                             :: iNTruncSt, outer_index, taupt, tauptmax, taustop, taucount
   integer                             :: iS1, iB1, iS2, iB2, vec_outer_ind
   real(KINDR)                         :: bra_t(DStates%NStatesMax), bra_t_logmax
   real(KINDR)                         :: ket_t(DStates%NStatesMax), ket_t_logmax
   real(KINDR)                         :: bra_int(DStates%NStatesMax), ket_int_logmax
   real(KINDR)                         :: ket_int(DStates%NStatesMax), bra_int_logmax
   real(KINDR)                         :: state_tp(DStates%NStatesMax), temp(DStates%NStatesMax)
   real(KINDR)                         :: tauR, tauL, tauD, tauint, tauref


   outer_index = this%sst_to_statesindex(this%outer_sst)
   tauptmax = size(ntau_n0(1, 1, 1, 1, :)) - 2
   tauD = this%beta/real(tauptmax + 1, KINDR)

   this%Trace = get_Trace_EB(this, DStates, this%postwin)
   call update_trace_EB(this)
   tauref = this%tauwin_max
   KetEnd => this%postwin

   RefOrbLoop: do iB1 = 1, DStates%NBands
      RefSpinLoop: do iS1 = 1, 2

         StateLoop: do iNTruncSt = this%outer_state, this%outer_state + this%outer_sst_size - 1
            if (this%b_statesampling) then
               vec_outer_ind = 1
            else
               vec_outer_ind = iNTruncSt
            end if

            taucount = 0
!=====================================================
! start with section at greater tau with missing kets
!=====================================================
            if (associated(KetEnd)) then
               Element => KetEnd%prev
            else
               Element => this%last
            end if

            ! load ket and general data for first section (between last ket and first bra)
            if (.not. associated(Element)) then
               sst = this%outer_sst
               NStates_sst = this%States(outer_index, 2)
               ket_t_logmax = 0.0_KINDR
               ket_t = 0.0_KINDR
               ket_t(iNTruncSt) = 1.0_KINDR
               tauR = 0.0_KINDR
               Element => this%first
            else
               sst = Element%following_sst
               NStates_sst = DStates%SubStates(sst)%NStates
               if (.not. Element%ket_state) stop "unreachable: wrong state type (1) in meas_ntau_n0"
               ket_t_logmax = Element%slogmax(vec_outer_ind)
               ket_t(1:NStates_sst) = Element%state(1:NStates_sst, vec_outer_ind)
               tauR = Element%tau
               Element => Element%next
            end if

            taupt = 0

            ! "insert n(0)" at tauref
            nsst = DStates%SubStates(sst)%Connect(iB1, iS1, 2)
            if (nsst == -1) cycle RefSpinLoop
            NStates_nsst = DStates%SubStates(nsst)%NStates

            ket_t(1:NStates_sst) =&
               exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                     - DStates%substates(sst)%Eval(0)) * (tauref - tauR))&
               * ket_t(1:NStates_sst)
            ket_t_logmax = ket_t_logmax&
               - (DStates%substates(sst)%Eval(0) - this%Egs) * (tauref - tauR)
            call split_log(ket_t(1:NStates_sst), ket_t_logmax)
            tauR = tauref

            call DGEMV('N', NStates_nsst, NStates_sst, 1d0,&
                       DStates%SubStates(sst)%Psis_EB(iB1, iS1, 2)%Op, NStates_nsst,&
                       ket_t, 1, 0d0,&
                       temp, 1)
            call DGEMV('N', NStates_sst, NStates_nsst, 1d0,&
                       DStates%SubStates(nsst)%Psis_EB(iB1, iS1, 1)%Op, NStates_sst,&
                       temp, 1, 0d0,&
                       ket_t, 1)
            call split_log(ket_t(1:NStates_sst), ket_t_logmax)


            ElementRLoop: do

               ! load stored bra
               if (.not. associated(Element)) then
                  bra_t_logmax = 0.0_KINDR
                  bra_t = 0.0_KINDR
                  bra_t(iNTruncSt) = 1.0_KINDR
                  tauL = this%beta
               else
                  if (Element%ket_state) stop "unreachable: wrong state type (2) in meas_ntau_n0"
                  bra_t_logmax = Element%slogmax(vec_outer_ind)
                  bra_t(1:NStates_sst) = Element%state(1:NStates_sst, vec_outer_ind)
                  tauL = Element%tau
               end if

               do
                  if (taupt > tauptmax) exit ElementRLoop
                  if (tauL < tauref + real(taupt, KINDR) * tauD) exit

                  taucount = taucount + 1

                  tauint = tauref + real(taupt, KINDR) * tauD
                  if (tauR > tauint .or. tauL < tauint) stop "unreachable: wrong section (1) in meas_ntau"

                  bra_int(1:NStates_sst)&
                     = bra_t(1:NStates_sst)&
                       * exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                               - DStates%substates(sst)%Eval(0)) * (-(tauint - tauL)))
                  bra_int_logmax = bra_t_logmax&
                     - (DStates%substates(sst)%Eval(0) - this%Egs) * (-(tauint - tauL))
                  call split_log(bra_int(1:NStates_sst), bra_int_logmax)

                  ket_int(1:NStates_sst)&
                     = ket_t(1:NStates_sst)&
                       * exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                               - DStates%substates(sst)%Eval(0)) * (tauint - tauR))
                  ket_int_logmax = ket_t_logmax&
                     - (DStates%substates(sst)%Eval(0) - this%Egs) * (tauint - tauR)
                  call split_log(ket_int(1:NStates_sst), ket_int_logmax)

                  ! write n(tau)n(0) for one tau-point
                  FlavourLoop: do iB2 = 1, DStates%NBands
                     do iS2 = 1, 2
                        nsst = DStates%SubStates(sst)%Connect(iB2, iS2, 2)
                        if (nsst == -1) cycle
                        NStates_nsst = DStates%SubStates(nsst)%NStates

                        call DGEMV('N', NStates_nsst, NStates_sst, 1d0,&
                                   DStates%SubStates(sst)%Psis_EB(iB2, iS2, 2)%Op, NStates_nsst,&
                                   ket_int, 1, 0d0,&
                                   state_tp, 1)

                        call DGEMV('N', NStates_sst, NStates_nsst, 1d0,&
                                   DStates%SubStates(nsst)%Psis_EB(iB2, iS2, 1)%Op, NStates_sst,&
                                   state_tp, 1, 0d0,&
                                   temp, 1)

                        ntau_n0(iB1, iS1, iB2, iS2, taupt)&
                           = ntau_n0(iB1, iS1, iB2, iS2, taupt)&
                             - dot_product(bra_int(1:NStates_sst), temp(1:NStates_sst))&
                               * trval(TLogTr(log = ket_int_logmax + bra_int_logmax,&
                                              sign = this%cfgsign) / this%Trace)
                     end do
                  end do FlavourLoop

                  taupt = taupt + 1
               end do

               if (.not. associated(Element)&
                   .or. tauref + real(taupt, KINDR) * tauD > this%beta) exit

               nsst = Element%following_sst
               NStates_nsst = DStates%SubStates(nsst)%NStates !rows

! exp(-tau H)ket
               state_tp(1:NStates_sst) =&
                  exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                        - DStates%substates(sst)%Eval(0)) * (Element%tau - tauR))&
                  * ket_t(1:NStates_sst)
               ket_t_logmax = ket_t_logmax&
                  -(DStates%substates(sst)%Eval(0) - this%Egs) * (Element%tau - tauR)
               call split_log(state_tp(1:NStates_sst), ket_t_logmax)

! O(^+)ket
!ket_t(0:DStates%SubStates(nsst)%NStates-1)=&
!&matmul(DStates%SubStates(sst)%Psis_EB(Element%Orbital,Element%Spin,Element%CA)%Op,&
!        ket_tp(0:DStates%SubStates(sst)%NStates-1))

               call DGEMV('N',NStates_nsst,NStates_sst,1d0,&
                          DStates%SubStates(sst)%Psis_EB(Element%Orbital,Element%Spin,Element%CA)%Op,&
                          NStates_nsst,&
                          state_tp,1,0d0,&
                          ket_t,1)
               call split_log(ket_t(1:NStates_nsst), ket_t_logmax)

               tauR = Element%tau
               sst = nsst
               NStates_sst = NStates_nsst
               Element => Element%next
            enddo ElementRLoop

            taustop = taupt - 1

!======================================================
! do sections at lesser tau with missing bras
!======================================================
            if (associated(KetEnd)) then
               Element => KetEnd
            else
               Element => null()
            end if

            ! load bra and general data for first section (between last ket and first bra again)
            if (.not. associated(Element)) then
               sst = this%outer_sst
               NStates_sst = this%States(outer_index, 2)
               bra_t_logmax = 0.0_KINDR
               bra_t = 0.0_KINDR
               bra_t(iNTruncSt) = 1.0_KINDR
               tauL = this%beta
               Element => this%last
            else
               sst = Element%preceding_sst
               NStates_sst = DStates%SubStates(sst)%NStates
               if (Element%ket_state) stop "unreachable: wrong state type (3) in meas_ntau_n0"
               bra_t_logmax = Element%slogmax(vec_outer_ind)
               bra_t(1:NStates_sst) = Element%state(1:NStates_sst, vec_outer_ind)
               tauL = Element%tau
               Element => Element%prev
            end if

            taupt = tauptmax

            ! first section + density insertion
! exp(-tau H)bra
            bra_t(1:NStates_sst) =&
               exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                     - DStates%substates(sst)%Eval(0)) * (-(tauref - tauL)))&
               * bra_t(1:NStates_sst)
            bra_t_logmax = bra_t_logmax&
               - (DStates%substates(sst)%Eval(0) - this%Egs) * (-(tauref - tauL))
            call split_log(bra_t(1:NStates_sst), bra_t_logmax)
            tauL = tauref

            ! "insert n(0)"
            nsst = DStates%SubStates(sst)%Connect(iB1, iS1, 2)
            NStates_nsst = DStates%SubStates(nsst)%NStates
            call DGEMV('N', NStates_nsst, NStates_sst, 1d0,&
                       DStates%SubStates(sst)%Psis_EB(iB1, iS1, 2)%Op, NStates_nsst,&
                       bra_t, 1, 0d0,&
                       temp, 1)
            call DGEMV('N', NStates_sst, NStates_nsst, 1d0,&
                       DStates%SubStates(nsst)%Psis_EB(iB1, iS1, 1)%Op, NStates_sst,&
                       temp, 1, 0d0,&
                       bra_t, 1)
            call split_log(bra_t(1:NStates_sst), bra_t_logmax)


            ElementLLoop: do

               ! load stored ket
               if (.not. associated(Element)) then
                  ket_t_logmax = 0.0_KINDR
                  ket_t = 0.0_KINDR
                  ket_t(iNTruncSt) = 1.0_KINDR
                  tauR = 0.0_KINDR
               else
                  if (.not. Element%ket_state)&
                     stop "unreachable: wrong state type (4) in meas_ntau_n0"
                  ket_t_logmax = Element%slogmax(vec_outer_ind)
                  ket_t(1:NStates_sst) = Element%state(1:NStates_sst, vec_outer_ind)
                  tauR = Element%tau
               end if

               do
                  if (taupt <= taustop) exit ElementLLoop
                  if (tauR > tauref - real(tauptmax + 1 - taupt, KINDR) * tauD) exit

                  taucount = taucount + 1

                  tauint = tauref - real(tauptmax + 1 - taupt, KINDR) * tauD
                  if (tauR > tauint .or. tauL < tauint) stop "unreachable: wrong section (2) in meas_ntau"

                  bra_int(1:NStates_sst)&
                     = bra_t(1:NStates_sst)&
                       * exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                               - DStates%substates(sst)%Eval(0)) * (-(tauint - tauL)))
                  bra_int_logmax = bra_t_logmax&
                     - (DStates%substates(sst)%Eval(0) - this%Egs) * (-(tauint - tauL))
                  call split_log(bra_int(1:NStates_sst), bra_int_logmax)

                  ket_int(1:NStates_sst)&
                     = ket_t(1:NStates_sst)&
                       * exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                               - DStates%substates(sst)%Eval(0)) * (tauint - tauR))
                  ket_int_logmax = ket_t_logmax&
                     - (DStates%substates(sst)%Eval(0) - this%Egs) * (tauint - tauR)
                  call split_log(ket_int(1:NStates_sst), ket_int_logmax)

                  ! write n(tau)n(0) at one tau-point
                  FlavourLoop2: do iB2 = 1, DStates%NBands
                     do iS2 = 1, 2
                        nsst = DStates%SubStates(sst)%Connect(iB2, iS2, 2)
                        if (nsst == -1) cycle
                        NStates_nsst = DStates%SubStates(nsst)%NStates

                        call DGEMV('N', NStates_nsst, NStates_sst, 1d0,&
                                   DStates%SubStates(sst)%Psis_EB(iB2, iS2, 2)%Op, NStates_nsst,&
                                   ket_int, 1, 0d0,&
                                   state_tp, 1)

                        call DGEMV('N', NStates_sst, NStates_nsst, 1d0,&
                                   DStates%SubStates(nsst)%Psis_EB(iB2, iS2, 1)%Op, NStates_sst,&
                                   state_tp, 1, 0d0,&
                                   temp, 1)

                        ntau_n0(iB1, iS1, iB2, iS2, taupt)&
                           = ntau_n0(iB1, iS1, iB2, iS2, taupt)&
                             - dot_product(bra_int(1:NStates_sst), temp(1:NStates_sst))&
                               * trval(TLogTr(log = ket_int_logmax + bra_int_logmax,&
                                              sign = this%cfgsign) / this%Trace)
                     end do
                  end do FlavourLoop2

                  taupt = taupt - 1
               end do

               if (.not. associated(Element)&
                   .or. tauref - real(tauptmax + 1 - taupt, KINDR) * tauD < 0) exit

               nsst = Element%preceding_sst
               NStates_nsst = DStates%SubStates(nsst)%NStates !rows

! exp(-tau H)bra
               state_tp(1:NStates_sst) =&
                  exp(-(DStates%substates(sst)%Eval(0:NStates_sst-1)&
                        - DStates%substates(sst)%Eval(0)) * (-(Element%tau - tauL)))&
                  * bra_t(1:NStates_sst)
               bra_t_logmax = bra_t_logmax&
                  -(DStates%substates(sst)%Eval(0) - this%Egs) * (-(Element%tau - tauL))
               call split_log(state_tp(1:NStates_sst), bra_t_logmax)

! O(^+)bra
!bra_t(0:DStates%SubStates(nsst)%NStates-1)=&
!matmul(DStates%SubStates(sst)%Psis_EB(ElementL%Orbital,ElementL%Spin,3-ElementL%CA)%Op,&
!       bra_tp(0:DStates%SubStates(sst)%NStates-1))

               call DGEMV('N',NStates_nsst,NStates_sst,1d0,&
                          DStates%SubStates(sst)%Psis_EB(Element%Orbital,Element%Spin,3-Element%CA)%Op,&
                          NStates_nsst,&
                          state_tp,1,0d0,&
                          bra_t,1)
               call split_log(bra_t(1:NStates_nsst), bra_t_logmax)

               tauL = Element%tau
               sst = nsst
               NStates_sst = NStates_nsst
               Element => Element%prev
            enddo ElementLLoop

            if (taucount /= tauptmax + 1)&
               stop "unreachable: wrote incorrect number of entries in meas_ntau_n0"
         end do StateLoop
      end do RefSpinLoop
   end do RefOrbLoop
end subroutine meas_ntau_n0

!===============================================================================
!> Dump sufficient information about the operator sequence into arrays to
!  recreate Z-space configurations.
subroutine dump_mc_config_into_arrays(this, taus, orbs, spins, cas, hashybs, outer_sst, outer_state)
!===============================================================================
   type(TTrace)                                        :: this
!output
   real(KINDR), dimension(:), intent(out)              :: taus
   integer, dimension(:), intent(out)                  :: orbs, spins, cas, hashybs
   integer, intent(out)                                :: outer_sst, outer_state
!local
   type(TOper), pointer                                :: ptr
   integer                                             :: i

   if (size(taus) /= this%NOper&
       .or. size(orbs) /= this%NOper&
       .or. size(spins) /= this%NOper&
       .or. size(cas) /= this%NOper&
       .or. size(hashybs) /= this%NOper) return
   outer_sst = this%outer_sst
   outer_state = this%outer_state
   ptr => this%first
   do i = 1, this%NOper
      taus(i) = ptr%tau
      orbs(i) = ptr%orbital
      spins(i) = ptr%spin
      cas(i) = ptr%ca
      if (ptr%has_hyb) then
         hashybs(i) = -1
      else
         hashybs(i) = 0
      end if
      ptr => ptr%next
   end do
end subroutine dump_mc_config_into_arrays

!===============================================================================
!> Create Z-space configuration from data arrays.
subroutine set_mc_config_from_arrays(this, DStates, taus, orbs, spins, cas,&
                                     hashybs, outer_sst, outer_state)
!===============================================================================
   type(TTrace)                                        :: this
   type(TStates)                                       :: DStates
!input
   real(KINDR), dimension(:), intent(in)               :: taus
   integer, dimension(:), intent(in)                   :: orbs, spins, cas, hashybs
   integer, intent(in)                                 :: outer_sst, outer_state
!local
   type(TOper), pointer                                :: ptr, lastptr
   type(TOperPointer), allocatable                     :: flast(:, :)
   integer                                             :: i, j
   type(TSubMatrix), pointer                           :: FullHybr(:, :)
   real(KINDR), pointer                                :: FullHybr_offdiag(:, :)

   ! clear current configuration
   ptr => this%first
   do i = 1, this%NOper
      lastptr => ptr
      ptr => ptr%next
      call recycle_oper(this, lastptr)
   end do
   this%first => null()
   this%last => null()
   this%prewin => null()
   this%last_prewin => null()
   this%postwin => null()
   this%UpdateL => null()
   this%UpdateR => null()
   this%LastEndR => null()
   do j = 1, size(this%fprewin, 2)
      do i = 1, size(this%fprewin, 1)
         this%fprewin(i, j)%p => null()
         this%ffirst(i, j)%p => null()
         this%flast(i, j)%p => null()
         this%NOSOper(i, j) = 0
         this%NOSCAOper(i, j, 1) = 0
         this%NOSCAOper(i, j, 2) = 0
      end do
   end do


   ! create new configuration
   allocate(flast(size(this%ffirst, 1), size(this%ffirst, 2)))
   this%outer_sst = outer_sst
   this%outer_state = outer_state
   this%window_position = 0
   this%NOper = size(taus)
   lastptr => null()

   if (this%NOper > 0) then
      call fresh_oper(this, ptr)
      this%first => ptr

      do i = 1, this%NOper
         ptr%tau = taus(i)
         ptr%orbital = orbs(i)
         ptr%spin = spins(i)
         ptr%ca = cas(i)
         if (hashybs(i) == 0) then
            ptr%has_hyb = .false.
         else
            ptr%has_hyb = .true.
         end if

         if (.not. associated(flast(ptr%orbital, ptr%spin)%p)) then
            this%ffirst(ptr%orbital, ptr%spin)%p => ptr
         else
            flast(ptr%orbital, ptr%spin)%p%fnext => ptr
         end if
         ptr%prev => lastptr
         ptr%fprev => flast(ptr%orbital, ptr%spin)%p
         ptr%fnext => null()
         flast(ptr%orbital, ptr%spin)%p => ptr
         if (ptr%has_hyb) then
            this%NOSOper(ptr%orbital, ptr%spin) =&
               this%NOSOper(ptr%orbital, ptr%spin) + 1
            this%NOSCAOper(ptr%orbital, ptr%spin, ptr%ca) =&
               this%NOSCAOper(ptr%orbital, ptr%spin, ptr%ca) + 1
         end if

         if (i < this%NOper) then
            call fresh_oper(this, ptr%next)
         else
            ptr%next => null()
         end if
         lastptr => ptr
         ptr => ptr%next
      end do

      this%last => lastptr
      do j = 1, size(this%flast, 2)
         do i = 1, size(this%flast, 1)
            this%flast(i, j)%p => flast(i, j)%p
         end do
      end do
   end if


   ! set some global state
   this%window_position = 0
   call set_window(this)
   this%NPairs = -1
   call update_PreWinFCount(this)
   this%ket_b2_invalid = .true.
   this%bra_b2_invalid = .true.

   ! calculate and store weight
   this%Trace = get_Trace_EB(this, DStates, global=.true.)
   call update_trace_EB(this, global=.true.)
   this%BosonicTrace = get_BosonicTrace(this, DStates)
   if (this%b_offdiag) then
      FullHybr_offdiag => get_FullHybr_offdiag(this, DStates)
      call get_MatLogDetFull(size(FullHybr_offdiag, 1), FullHybr_offdiag, this%Det)
      if (associated(this%MInv_full%Mat)) deallocate(this%MInv_full%Mat)
      this%MInv_full%Mat => FullHybr_offdiag
   else
      FullHybr => get_InvFullHybr(this, DStates, this%Det)
      do j = 1, size(this%MInv, 2)
         do i = 1, size(this%MInv, 1)
            if (associated(this%MInv(i, j)%Mat)) deallocate(this%MInv(i, j)%Mat)
         end do
      end do
      if (associated(this%MInv)) deallocate(this%MInv)
      this%MInv => FullHybr
   end if
   this%cfgsign = -get_Sign(this, DStates)
end subroutine set_mc_config_from_arrays

!===============================================================================
subroutine fresh_oper(this, oper)
!===============================================================================
   type(TTrace)                      :: this
!output
   type(TOper), pointer, intent(out) :: oper

   if (this%iOperPool == 0) then
      allocate(oper)
      if (this%b_fix_prealloc_stvec) then
         allocate(oper%state(this%NStatesMax, this%NTruncStatesMax))
         allocate(oper%cache(this%NStatesMax, this%NTruncStatesMax))
         allocate(oper%slogmax(this%NTruncStatesMax))
         allocate(oper%clogmax(this%NTruncStatesMax))
      end if
      oper%calc_new = .true.
   else
      oper => this%OperPool(this%iOperPool)%p
      this%OperPool(this%iOperPool)%p => null()
      this%iOperPool = this%iOperPool - 1
      oper%calc_new = .true.
      oper%cache_written = .false.
   end if
end subroutine fresh_oper

!===============================================================================
subroutine duplicate_oper(this, opera, operb)
!===============================================================================
   type(TTrace)                      :: this
!output
   type(TOper), pointer, intent(in)  :: opera
   type(TOper), pointer, intent(out) :: operb

   call fresh_oper(this, operb)
   operb%next => opera%next
   operb%fnext => opera%fnext
   operb%prev => opera%prev
   operb%fprev => opera%fprev
   operb%tau = opera%tau
   operb%Orbital = opera%Orbital
   operb%Spin = opera%Spin
   operb%CA = opera%CA
   operb%has_hyb = opera%has_hyb
   operb%calc_new = opera%calc_new
   operb%following_sst = opera%following_sst
   operb%preceding_sst = opera%preceding_sst
   operb%cache_post_sst = opera%cache_post_sst
   operb%cache_pre_sst = opera%cache_pre_sst
   operb%ket_state = opera%ket_state
   operb%ket_cache = opera%ket_cache
   operb%cache_written = opera%cache_written
   if (allocated(opera%state)) then
      operb%state = opera%state(:, :)
   end if
   if (allocated(opera%cache)) then
      operb%cache = opera%cache(:, :)
   end if
   if (allocated(opera%slogmax)) then
      operb%slogmax = opera%slogmax(:)
   end if
   if (allocated(opera%clogmax)) then
      operb%clogmax = opera%clogmax(:)
   end if
end subroutine duplicate_oper

!===============================================================================
subroutine recycle_oper(this, oper)
!===============================================================================
   type(TTrace)                     :: this
!output
   type(TOper), pointer, intent(inout) :: oper

   if (.not. this%b_fix_prealloc_stvec) then
   if (allocated(oper%state)) deallocate(oper%state)
   if (allocated(oper%cache)) deallocate(oper%cache)
   if (allocated(oper%slogmax)) deallocate(oper%slogmax)
   if (allocated(oper%clogmax)) deallocate(oper%clogmax)
   end if
   oper%calc_new=.true.
   oper%cache_written = .false.
   this%iOperPool = this%iOperPool + 1
   this%OperPool(this%iOperPool)%p => oper
end subroutine recycle_oper

!===============================================================================
!> Small helper function for debugging. Writes the current configuration to the
!! file fort.666
subroutine print_Trace(this, file, spur)
!===============================================================================
   type(TTrace), intent(in)            :: this
!local
   type(TOper),pointer                 :: Element
   integer                             :: iOper, file
   real(kindr) :: spur

   iOper=1
   Element=>this%first
   write(file,*)this%NOper," elements in the trace"
   write(unit=file,fmt=*) "trace", spur
   do while(associated(Element))
      write(file,*)iOper,Element%tau,Element%Orbital,Element%Spin,Element%CA
      Element=>Element%next
      iOper=iOper+1
   enddo
end subroutine print_Trace

!===============================================================================
!> Small helper function for debugging. Writes the current configuration to the
!! file fort.666
subroutine print_Trace_screen(this)
!===============================================================================
   type(TTrace), intent(in)            :: this
!local
   type(TOper),pointer                 :: Element
   integer                             :: iOper

   iOper=1
   Element=>this%first
   write(*, "(I5, ' ops, win ', F6.1, ' - ', F6.1, ', NPairs ', I5)")&
      this%NOper, this%tauwin_min, this%tauwin_max, this%NPairs
   do while(associated(Element))
      write(*, "(I5, ' ', L5, ' ', F12.7, I2, I2, I2, ' ', L5,&
                &' ', I4, I4, ' ', L5, ' c: ', L5, ' ', I4, I4, ' ', L5)")&
                iOper, Element%calc_new,&
                Element%tau, Element%Orbital, Element%Spin, Element%CA, Element%has_hyb,&
                Element%preceding_sst, Element%following_sst, Element%ket_state,&
                Element%cache_written,&
                Element%cache_pre_sst, Element%cache_post_sst, Element%ket_cache
      Element=>Element%next
      iOper=iOper+1
   enddo
   
end subroutine print_Trace_screen

!===============================================================================
!> @brief
!> Performs an evolution.
!>
!> @param substate[in]
!> @param state_vector[inout]
!> @param state_vector_evolved[in]
!> @param slen[in] the length of the state arrays
!> @param tau[in]
!> @param egs[in]
!===============================================================================
subroutine time_evolve_in_eigenbasis(substate, state_vector, state_vector_evolved, slen, tau, egs)
!===============================================================================
!input
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

end subroutine time_evolve_in_eigenbasis


!===============================================================================
!> Small helper function for debugging. Writes the current configuration to the
!! file fort.666
subroutine print_Trace_flavour(this,dstates,orb,spin)
!===============================================================================
   type(TTrace)                        :: this
   type(TStates)                       :: DStates
!local
   type(TOper),pointer                 :: Element
   integer                             :: iOper,ib,is,orb,spin

   if(orb.eq.0)then
   do ib=1,DStates%nbands
   do is=1,2
      iOper=1
      Element=>this%ffirst(ib,is)%p
      write(*,*)this%NOSOper(ib,is)," elements of flavour     ", ib, is, " in the trace"
      do while(associated(Element))
         write(*,*)iOper,Element%tau,Element%Orbital,Element%Spin,Element%CA
         Element=>Element%fnext
         iOper=iOper+1
      enddo
   enddo
   enddo
   endif

   do ib=1,DStates%nbands
   do is=1,2
   if(ib.ne.orb) cycle
   if(is.ne.spin) cycle
      iOper=1
      Element=>this%ffirst(ib,is)%p
      write(*,*)this%NOSOper(ib,is)," elements of flavour     ", ib, is, " in the trace"
      do while(associated(Element))
         write(*,*)iOper,Element%tau,Element%Orbital,Element%Spin,Element%CA
         Element=>Element%fnext
         iOper=iOper+1
      enddo
   enddo
   enddo
   
end subroutine print_Trace_flavour






!=====================  here come the segment functions ! ======================
!===============================================================================
!===============================================================================
!===============================================================================
!===============================================================================

!!!qn checks

!===============================================================================
!===============================================================================
!===============================================================================
!===============================================================================


!===============================================================================
!> Checks if current configuration is allowed, considering also the saved outer
! state in segment picture
logical function check_qn_add_eachflavour(this,band,spin)
!===============================================================================
   type(TTrace)                  :: this
!input
   integer                       :: band,spin
!local
   integer                       :: ca
   type(TOper),pointer           :: Element

   check_qn_add_eachflavour=.true.

   !!! initial set of "ca" according to saved outer state
   if(this%initial_outerstate_old(band,spin).eqv..true.)then
      ca=1
   else
      ca=2
   endif

   Element=>this%ffirst(band,spin)%p
   !!! go through linked list for the flavour of the inserted operator, and
   !!! check if creators and annihilator appear in alternating order
   do while(associated(Element))
      if(element%ca.eq.ca)then
         check_qn_add_eachflavour=.false.
      else
         ca=element%ca
      endif
      !!! go to next element of same flavour
      Element=>Element%fnext
   enddo

end function check_qn_add_eachflavour


!===============================================================================
!> Calculates additional sign for segment calculation, that comes from operators
!ordering in the impurity states. In Matrix algorithm this sign is in the
!matrices of the operators
real(kindr) function get_sign_seg(this,oper)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TOperPointer)                  :: Oper(2)
!local
   type(TOper),pointer                 :: Element
   integer :: nperms
   logical :: meas_signum

   nperms=0
   meas_signum=.false.

   Element=>this%first
   !!! First we move the pair considered the the position where the operator
   !with smaller tau belongs to; this gives no sign since it is a boson.
   !Then the second operator is moved to its position, and the number of
   !permutations with other operators counted. -> sign
   do while(associated(element))
      if(associated(oper(2)%p,element))meas_signum=.false.
      if(meas_signum)then
         nperms=nperms+1
      endif
      if(associated(oper(1)%p,element))meas_signum=.true.
      Element=>Element%next
   enddo

   get_sign_seg=(-1)**(nperms)

end function get_sign_seg


!===============================================================================
!> compares the first significant figures of two large numbers
logical function check_equal_times(this,tau1,tau2)
!===============================================================================
   type(TTrace)                  :: this
   
   type(toper),pointer           :: element
   real(kindr) :: tau1, tau2

   check_equal_times=.false.

   element=>this%first
   do while(associated(element))

      if(element%tau.eq.tau1)then
         check_equal_times=.true.
         return
      endif
      if(element%tau.eq.tau2)then
         check_equal_times=.true.
         return
      endif

      element=>element%next
   enddo

end function check_equal_times 


!===============================================================================
!> compares the first significant figures of two large numbers
logical function update_noscaoper(this,dstates)
!===============================================================================
   type(TTrace)                  :: this
   type(TStates)                 :: DStates
   
   type(toper),pointer           :: element
   integer :: b,s

   update_noscaoper=.false.

   !!! update the noscaoper array
   this%NOSCAOper=0
   element=>this%first
   do while(associated(element))

      if(element%has_hyb) then
         this%NOSCAOper(element%orbital,element%spin,element%ca)=&
         this%NOSCAOper(element%orbital,element%spin,element%ca) + 1
      endif

      element=>element%next
   enddo

   !!! print out special cases
   do b=1,dstates%nbands
   do s=1,2
      
      !!! if there are more creators than annihilators
      if(abs(this%noscaoper(b,s,1)-this%noscaoper(b,s,2)).gt.2)then
         write(*,*) " "
         write(*,*) "hier!"
         write(*,*) "this%noscaoper ca=1",b,s,  this%noscaoper(b,s,1)
         write(*,*) "this%noscaoper ca=2",b,s,  this%noscaoper(b,s,2)
      endif
      
   enddo
   enddo

end function update_noscaoper 


!===============================================================================
!===============================================================================
!===============================================================================
!===============================================================================

!!!propose and insert / remove functions

!===============================================================================
!===============================================================================
!===============================================================================
!===============================================================================


!===============================================================================
!> not used yet
subroutine propose_insert(this,DStates,oper)
!===============================================================================
   type(TTrace)                  :: this
   type(TStates)                 :: DStates
   type(TOperPointer)            :: Oper(2),tmp
   real(kindr) :: rand

   rand=grnd()
   Oper(1)%p%tau=rand*this%beta
   Oper(1)%p%orbital=randint(1, DStates%NBands)
   Oper(1)%p%spin=randint(1, 2)
   Oper(1)%p%ca=1

   rand=grnd()
   Oper(2)%p%tau=rand*this%beta
   Oper(2)%p%orbital=Oper(1)%p%orbital
   Oper(2)%p%spin=Oper(1)%p%spin
   Oper(2)%p%ca=2

   if(Oper(1)%p%tau.gt.Oper(2)%p%tau)then
      tmp%p=>oper(1)%p
      oper(1)%p=>oper(2)%p
      oper(2)%p=>tmp%p
   endif
   if(Oper(1)%p%tau.eq.Oper(2)%p%tau)then
      Oper(2)%p%tau=Oper(2)%p%tau+this%equal_time_offset
   endif

end subroutine propose_insert

!===============================================================================
!> Inserts two operators into the full trace and the trace of that specific
!flavour.
subroutine my_insert_Oper(this,Oper)
!===============================================================================
   type(TTrace)                  :: this
!input
   type(TOperPointer)            :: Oper(2)
   type(toper),pointer           :: element
!local
   integer                       :: i

   !!! When is is the first operator of that flavour, begin a new list and
   !!! set the pointers to the beginning and the end of it.
   !!! Set the connections and null pointers, where no connections are.
   if(this%nosoper(oper(1)%p%orbital,oper(1)%p%spin).eq.0)then

      this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p=>Oper(1)%p
      this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p=>Oper(2)%p
      this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p%fnext=>this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p
      this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p%fprev=>this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p
      this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p%fprev=>null()
      this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p%fnext=>null()

   !!! When a list of that flavour exists
   else

      !!! loop over the 2 operators
      do i=1,2

         !!! is operator new last element of the list?
         if(Oper(i)%p%tau.ge.this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p%tau)then
            !!! be always careful if the operator has the same tau than an operator existing in the list! then shift.
            !if(Oper(i)%p%tau.eq.this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p%tau)Oper(i)%p%tau=Oper(i)%p%tau+this%equal_time_offset
            this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p%fnext=>Oper(i)%p
            Oper(i)%p%fprev=>this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p
            Oper(i)%p%fnext=>null()
            this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p=>Oper(i)%p
         !!! is operator new first element of the list?
         elseif(Oper(i)%p%tau.le.this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p%tau)then
            !if(Oper(i)%p%tau.eq.this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p%tau)Oper(i)%p%tau=Oper(i)%p%tau-this%equal_time_offset
            this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p%fprev=>Oper(i)%p
            Oper(i)%p%fprev=>null()
            Oper(i)%p%fnext=>this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p
            this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p=>Oper(i)%p
         !!! operator is in the middle
         else
            Element=>this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p
            FindPositionLoop1: do while(associated(Element))
               !if(Oper(i)%p%tau.eq.Element%tau)Oper(i)%p%tau=Oper(i)%p%tau+this%equal_time_offset
               if(Oper(i)%p%tau.lt.Element%tau)then
                  Element%fprev%fnext=>Oper(i)%p
                  Oper(i)%p%fprev=>Element%fprev
                  Element%fprev=>Oper(i)%p
                  Oper(i)%p%fnext=>Element
                  exit FindPositionLoop1
               endif

               Element=>Element%fnext
            enddo FindPositionLoop1
         endif

      enddo

   endif


   !!! the same done for the list of all flavours time ordered.
   if(this%noper.eq.0)then

      this%first=>Oper(1)%p
      this%last=>Oper(2)%p
      this%first%next=>this%last
      this%last%prev=>this%first
      this%first%prev=>null()
      this%last%next=>null()

   else

      do i=1,2

         if(Oper(i)%p%tau.ge.this%last%tau)then
            !if(Oper(i)%p%tau.eq.this%last%tau)Oper(i)%p%tau=Oper(i)%p%tau+this%equal_time_offset
            this%last%next=>Oper(i)%p
            Oper(i)%p%prev=>this%last
            Oper(i)%p%next=>null()
            this%last=>Oper(i)%p
         elseif(Oper(i)%p%tau.le.this%first%tau)then
            !if(Oper(i)%p%tau.eq.this%first%tau)Oper(i)%p%tau=Oper(i)%p%tau-this%equal_time_offset
            this%first%prev=>Oper(i)%p
            Oper(i)%p%prev=>null()
            Oper(i)%p%next=>this%first
            this%first=>Oper(i)%p
         else
            Element=>this%first
            FindPositionLoop2: do while(associated(Element))
               !if(Oper(i)%p%tau.eq.Element%tau)Oper(i)%p%tau=Oper(i)%p%tau+this%equal_time_offset
               if(Oper(i)%p%tau.lt.Element%tau)then
                  Element%prev%next=>Oper(i)%p
                  Oper(i)%p%prev=>Element%prev
                  Element%prev=>Oper(i)%p
                  Oper(i)%p%next=>Element
                  exit FindPositionLoop2
               endif

               Element=>Element%next
            enddo FindPositionLoop2
         endif

      enddo

   endif
   
   this%NOper=this%NOper+2
   !we only count NOSOper over hybridization operators
   if(Oper(1)%p%has_hyb.and.Oper(2)%p%has_hyb) then
      this%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)=this%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)+1
      this%NOSCAOper(Oper(1)%p%Orbital,Oper(1)%p%Spin,Oper(1)%p%CA) =&
         this%NOSCAOper(Oper(1)%p%Orbital,Oper(1)%p%Spin,Oper(1)%p%CA) + 1
      this%NOSOper(Oper(2)%p%Orbital,Oper(2)%p%Spin)=this%NOSOper(Oper(2)%p%Orbital,Oper(2)%p%Spin)+1
      this%NOSCAOper(Oper(2)%p%Orbital,Oper(2)%p%Spin,Oper(2)%p%CA) =&
         this%NOSCAOper(Oper(2)%p%Orbital,Oper(2)%p%Spin,Oper(2)%p%CA) + 1
   endif

end subroutine my_insert_Oper

!===============================================================================
!> Removes two operators from the full trace and the trace of that specific
!flavour.
subroutine my_remove_Oper(this,Oper)
!===============================================================================
   type(TTrace)                  :: this
!input
   type(TOperPointer)            :: Oper(2)
   type(toper),pointer           :: element
!local
   integer                       :: i

   !!! When operator pair is the only one of that flavour, reset list to empty list.
   if(this%nosoper(oper(1)%p%orbital,oper(1)%p%spin).eq.2)then
      this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p%fnext=>null()
      this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p%fprev=>null()
      this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p=>null()
      this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p=>null()
   else
      !!! loop over the 2 operators
      do i=1,2
         !!! if operator is the first of the list, make second of list 
         !!! to the first and set pointers "fnext" and "fprev" right.
         if(associated(Oper(i)%p,this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p))then
            Oper(i)%p%fnext%fprev=>null()
            this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p=>Oper(i)%p%fnext
            Oper(i)%p%fnext=>null()
         !!! if operator is the last of the list
         elseif(associated(Oper(i)%p,this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p))then
            Oper(i)%p%fprev%fnext=>null()
            this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p=>Oper(i)%p%fprev
            Oper(i)%p%fprev=>null()
         !!! if operator is in the middle of the list
         else
            Element=>this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p
            RemoveOperLoop: do while(associated(Element))
               if(Oper(i)%p%tau.eq.Element%tau)then
                  Element%fnext%fprev=>Element%fprev
                  Element%fprev%fnext=>Element%fnext
                  exit RemoveOperLoop
               endif
               Element=>Element%fnext
            enddo RemoveOperLoop
         endif
      enddo
   endif

   !!! the same for list of all flavours
   if(this%noper.eq.2)then
      this%first%next=>null()
      this%last%prev=>null()
      this%first=>null()
      this%last=>null()
   else
      do i=1,2
         if(associated(Oper(i)%p,this%first))then
            Oper(i)%p%next%prev=>null()
            this%first=>Oper(i)%p%next
            Oper(i)%p%next=>null()
         elseif(associated(Oper(i)%p,this%last))then
            Oper(i)%p%prev%next=>null()
            this%last=>Oper(i)%p%prev
            Oper(i)%p%prev=>null()
         else
            Element=>this%first
            RemoveOperLoop2: do while(associated(Element))
               if(Oper(i)%p%tau.eq.Element%tau)then
                  Element%next%prev=>Element%prev
                  Element%prev%next=>Element%next
                  exit RemoveOperLoop2
               endif
               Element=>Element%next
            enddo RemoveOperLoop2
         endif
      enddo
   endif

   this%NOper=this%NOper-2
   !TODO: only counting hybs to NOSOper, is this smart?
   if(Oper(1)%p%has_hyb.and.Oper(2)%p%has_hyb) then
      this%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)=this%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)-1
      this%NOSCAOper(Oper(1)%p%Orbital,Oper(1)%p%Spin,Oper(1)%p%CA) =&
         this%NOSCAOper(Oper(1)%p%Orbital,Oper(1)%p%Spin,Oper(1)%p%CA) - 1
      this%NOSOper(Oper(2)%p%Orbital,Oper(2)%p%Spin)=this%NOSOper(Oper(2)%p%Orbital,Oper(2)%p%Spin)-1
      this%NOSCAOper(Oper(2)%p%Orbital,Oper(2)%p%Spin,Oper(2)%p%CA) =&
         this%NOSCAOper(Oper(2)%p%Orbital,Oper(2)%p%Spin,Oper(2)%p%CA) - 1
   endif

end subroutine my_remove_Oper

!===============================================================================
!> Proposed the insert of a segment or antisegment in the trace
logical function propose_insert_seg(this,DStates,oper,lmax,seg,overbeta,check_qn_additionally)
!===============================================================================
   type(TTrace)                  :: this
   type(TStates)                 :: DStates
   type(TOperPointer)            :: Oper(2)
   real(kindr)                   :: lmax,tauprev
   logical :: seg,segoverbeta,overbeta,check_qn_additionally
   integer :: nosoper,orb,spin
   type(TOper),pointer  :: elementf,elementl,element

   !!! random choose if segment or antisegment to be inserted
   if(grnd()<0.5D0)then
      seg=.true.
   else
      seg=.false.
   endif

   !!! flag to check quantum numbers additionally, in case a dangerous 
   !!! insert is going to be proposed
   check_qn_additionally=.false.

   propose_insert_seg=.true.

   !!! does the segment / antisegment go over beta?
   overbeta=.false.

   !!! TODO: find out if this is a good idea ?!?
   this%outerstate_full=.false.
   this%outerstate_empty=.false.
   
   !!! randomly coose properties of operators
   Oper(1)%p%tau=grnd()*this%beta
   Oper(1)%p%orbital=randint(1, DStates%NBands)
   Oper(1)%p%spin=randint(1, 2)
   orb=Oper(1)%p%orbital
   spin=Oper(1)%p%spin

   nosoper=this%nosoper(oper(1)%p%orbital,oper(1)%p%spin)

   Oper(2)%p%orbital=Oper(1)%p%orbital
   Oper(2)%p%spin=Oper(1)%p%spin
   Oper(2)%p%tau=0d0

   if(seg.eqv..true.)then
      !!! insert segment

      oper(1)%p%ca=1
      oper(2)%p%ca=2

      if(nosoper.eq.0)then
         !!! no (anti)segments in the trace

         if(this%initial_outerstate_old(orb,spin).eqv..false.)then
            !!! outer state unoccupied

            lmax=this%beta

            Oper(2)%p%tau=oper(1)%p%tau+lmax*grnd()

            !if(oper(1)%p%tau.eq.oper(2)%p%tau)then
               !oper(2)%p%tau=oper(2)%p%tau+this%equal_time_offset
            !endif

            if(oper(2)%p%tau.gt.this%beta)then
               oper(2)%p%tau=oper(2)%p%tau-this%beta
               overbeta=.true.
            endif

            return

         else
            !!! outer state occupied
            !!! proposed insert of segment over another segment
            !!! reject

            propose_insert_seg=.false.
            return

         endif


      endif

      if(nosoper.gt.0)then
         !!! segments in the trace

         elementf=>this%ffirst(orb,spin)%p
         elementl=>this%flast(orb,spin)%p

         !!! check if there is a segment or antisegment at tau=0
         if(elementf%ca.eq.1)then
            segoverbeta=.false.
         else
            segoverbeta=.true.
         endif

         !!! in front of all segments
         if(oper(1)%p%tau.le.elementf%tau)then
            
            !!!! dangerous cases!
            !if(oper(1)%p%tau.eq.elementf%tau)then
               !oper(1)%p%tau=oper(1)%p%tau-this%equal_time_offset
               !check_qn_additionally=.true.
            !endif
            !if(oper(1)%p%tau.eq.0d0)then
               !oper(1)%p%tau=this%equal_time_offset
               !check_qn_additionally=.true.
            !endif

            !!! if there is a segment over beta, reject
            if(segoverbeta)then
               lmax=0d0
               propose_insert_seg=.false.
               return
            else
               !!! insert segment in front of all, no segoverbeta
               !!! lmax is the intervall that the second operator (= end of
               !!! segment) can be inserted. E.i. one actually had to insert the
               !!! second operator in the intervall [0,beta], but because of
               !!! quantum numbers it is possible to calculate the conditional
               !!! probability of inserting the second operator in the allowed
               !!! intervall lmax as lmax/beta.
               lmax=elementf%tau-oper(1)%p%tau
               Oper(2)%p%tau=oper(1)%p%tau+grnd()*lmax

               !!!! dangerous cases!
               !if(oper(2)%p%tau.eq.oper(1)%p%tau)then
                  !oper(2)%p%tau=oper(2)%p%tau+this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif
               !if(oper(2)%p%tau.eq.elementf%tau)then
                  !oper(2)%p%tau=oper(2)%p%tau-this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif

               return
            endif

         !!! after all segments
         elseif(oper(1)%p%tau.ge.elementl%tau)then

            !!!! dangerous!
            !if(oper(1)%p%tau.eq.elementl%tau)then
               !oper(1)%p%tau=oper(1)%p%tau+this%equal_time_offset
               !check_qn_additionally=.true.
            !endif

            !!! if there is already a segment, reject
            if(segoverbeta)then
               lmax=0d0
               propose_insert_seg=.false.
               return
            else
               !!! insert segment after last segment
               lmax=this%beta-oper(1)%p%tau+elementf%tau
               Oper(2)%p%tau=grnd()*lmax+oper(1)%p%tau

               if(oper(2)%p%tau.gt.this%beta)then
                  oper(2)%p%tau=oper(2)%p%tau-this%beta
                  overbeta=.true.
               endif

               !if(oper(2)%p%tau.eq.oper(1)%p%tau)then
                  !oper(2)%p%tau=oper(2)%p%tau+this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif
               !if(oper(2)%p%tau.eq.elementf%tau)then
                  !oper(2)%p%tau=oper(2)%p%tau-this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif

               return
            endif

         !!! insert segment in the middle
         else

            tauprev=0d0
            element=>elementf
            findloop:do while(associated(element))

               !if(oper(1)%p%tau.eq.tauprev)then
                  !oper(1)%p%tau=oper(1)%p%tau+this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif

               !!! find position of new segment, and its neighbouring segments
               if((oper(1)%p%tau.gt.tauprev).and.(oper(1)%p%tau.lt.element%tau))then
                  if(element%ca.eq.1)then
                     lmax=element%tau-oper(1)%p%tau
                     Oper(2)%p%tau=oper(1)%p%tau+grnd()*lmax

                     !if(oper(2)%p%tau.eq.element%tau)then
                        !oper(2)%p%tau=oper(2)%p%tau-this%equal_time_offset
                        !check_qn_additionally=.true.
                     !endif
                     return
                  else
                     lmax=0d0
                     propose_insert_seg=.false.
                     return
                  endif
               endif

               tauprev=element%tau
               element=>element%fnext

            enddo findloop

         endif

      endif


   !!! the same for antisegments
   else
      !if(www)write(*,*) "--> seg eq false"

      oper(1)%p%ca=2
      oper(2)%p%ca=1

      if(nosoper.eq.0)then
         !if(www)write(*,*) "--> nosoper eq 0"

         if(this%initial_outerstate_old(orb,spin).eqv..true.)then
            !if(www)write(*,*)"--> outer state occupied"

            lmax=this%beta

            Oper(2)%p%tau=oper(1)%p%tau+grnd()*lmax

            if(oper(2)%p%tau.gt.this%beta)then
               oper(2)%p%tau=oper(2)%p%tau-this%beta
               overbeta=.true.
            endif


            !if(www)write(*,*) "oper(1)%p%tau", oper(1)%p%tau
            !if(www)write(*,*) "oper(2)%p%tau", oper(2)%p%tau
            !if(www)write(*,*) "oper(1)%p%ca", oper(1)%p%ca
            !if(www)write(*,*) "oper(2)%p%ca", oper(2)%p%ca
            !if(www)write(*,*) "lmax", lmax
            !if(www)write(*,*)"1008"
            return

         else
            !if(www)write(*,*)"--> outer state unoccupied"

            propose_insert_seg=.false.
            !if(www)write(*,*)"qn rejection"
            !if(www)write(*,*)"10088"
            return

         endif

      endif

      if(nosoper.gt.0)then
         !if(www)write(*,*) "--> nosoper gt 0"

         !!! check if segment is allowed at that position
         elementf=>this%ffirst(oper(1)%p%orbital,oper(1)%p%spin)%p
         elementl=>this%flast(oper(1)%p%orbital,oper(1)%p%spin)%p

         if(elementf%ca.eq.1)then
            segoverbeta=.false.
         else
            segoverbeta=.true.
         endif

         if(oper(1)%p%tau.lt.elementf%tau)then

            !if(oper(1)%p%tau.eq.elementf%tau)then
               !oper(1)%p%tau=oper(1)%p%tau-this%equal_time_offset
               !check_qn_additionally=.true.
            !endif
            !if(oper(1)%p%tau.eq.0d0)then
               !oper(1)%p%tau=this%equal_time_offset
               !check_qn_additionally=.true.
            !endif

            if(.not.segoverbeta)then
               lmax=0d0
               !if(www)write(*,*) "in front of all"
               !if(www)write(*,*) "qn rejection"
               !if(www)write(*,*) "1009"
               propose_insert_seg=.false.
               return
            else
               !!! insert antisegment in front of all, segoverbeta
               lmax=elementf%tau-oper(1)%p%tau
               Oper(2)%p%tau=oper(1)%p%tau+grnd()*lmax

               !if(oper(2)%p%tau.eq.oper(1)%p%tau)then
                  !oper(2)%p%tau=oper(2)%p%tau+this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif
               !if(oper(2)%p%tau.eq.elementf%tau)then
                  !oper(2)%p%tau=oper(2)%p%tau-this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif

               !if(www)write(*,*) "in front of all"
               !if(www)write(*,*) "oper(1)%p%tau", oper(1)%p%tau
               !if(www)write(*,*) "oper(2)%p%tau", oper(2)%p%tau
               !if(www)write(*,*) "oper(1)%p%ca", oper(1)%p%ca
               !if(www)write(*,*) "oper(2)%p%ca", oper(2)%p%ca
               !if(www)write(*,*) "lmax", lmax
               !if(www)write(*,*) "1010"
               return
            endif

            
         elseif(oper(1)%p%tau.gt.elementl%tau)then
            
            !if(oper(1)%p%tau.eq.elementl%tau)then
               !oper(1)%p%tau=oper(1)%p%tau+this%equal_time_offset
               !check_qn_additionally=.true.
            !endif

            if(.not.segoverbeta)then
               lmax=0d0
               !if(www)write(*,*) "in back of all"
               !if(www)write(*,*) "qn rejection"
               !if(www)write(*,*) "1011"
               propose_insert_seg=.false.
               return
            else
               !!! insert antisegment after last operator, segoverbeta
               lmax=this%beta-oper(1)%p%tau+elementf%tau
               Oper(2)%p%tau=oper(1)%p%tau+grnd()*lmax

               !if(oper(2)%p%tau.lt.oper(1)%p%tau)then
                  !oper(2)%p%ca=2
                  !oper(1)%p%ca=1
               !endif

               if(oper(2)%p%tau.gt.this%beta)then
                  oper(2)%p%tau=oper(2)%p%tau-this%beta
                  overbeta=.true.
               endif

               !if(oper(2)%p%tau.eq.oper(1)%p%tau)then
                  !oper(2)%p%tau=oper(2)%p%tau+this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif
               !if(oper(2)%p%tau.eq.elementf%tau)then
                  !oper(2)%p%tau=oper(2)%p%tau-this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif

               !if(www)write(*,*) "in back of all"
               !if(www)write(*,*) "oper(1)%p%tau", oper(1)%p%tau
               !if(www)write(*,*) "oper(2)%p%tau", oper(2)%p%tau
               !if(www)write(*,*) "oper(1)%p%ca", oper(1)%p%ca
               !if(www)write(*,*) "oper(2)%p%ca", oper(2)%p%ca
               !if(www)write(*,*) "lmax", lmax
               !if(www)write(*,*) "1012"
               return
            endif

         else

            tauprev=0d0
            element=>elementf
            findloop2:do while(associated(element))

               !if(oper(1)%p%tau.eq.tauprev)then
                  !oper(1)%p%tau=oper(1)%p%tau+this%equal_time_offset
                  !check_qn_additionally=.true.
               !endif

               if((oper(1)%p%tau.gt.tauprev).and.(oper(1)%p%tau.lt.element%tau))then
                  if(element%ca.eq.2)then
                     lmax=element%tau-oper(1)%p%tau
                     Oper(2)%p%tau=oper(1)%p%tau+grnd()*lmax

                     !if(oper(2)%p%tau.eq.element%tau)then
                        !oper(2)%p%tau=oper(2)%p%tau-this%equal_time_offset
                        !check_qn_additionally=.true.
                     !endif

                     !if(www)write(*,*) "in middle"
                     !if(www)write(*,*) "oper(1)%p%tau", oper(1)%p%tau
                     !if(www)write(*,*) "oper(2)%p%tau", oper(2)%p%tau
                     !if(www)write(*,*) "oper(1)%p%ca", oper(1)%p%ca
                     !if(www)write(*,*) "oper(2)%p%ca", oper(2)%p%ca
                     !if(www)write(*,*) "lmax", lmax
                     !if(www)write(*,*) "1013"
                     return
                  else
                     lmax=0d0
                     !if(www)write(*,*) "in middle"
                     !if(www)write(*,*) "qn rejection"
                     !if(www)write(*,*) "1014"
                     propose_insert_seg=.false.
                     return
                  endif
               endif
                  
               tauprev=element%tau
               element=>element%fnext

            enddo findloop2

         endif

      endif

   endif

end function propose_insert_seg


!===============================================================================
subroutine get_fpos_insert(this,Oper,FPos)
!===============================================================================
   type(TTrace)                        :: this
!output
   type(TOperPointer)                  :: Oper(2)
   integer                             :: FPos(2)
!local
   type(TOper),pointer                 :: Element

   Fpos(:)=0
   Element=>this%first      ! iterates over c^+(tau_2)

   FINDLOOP: do while(associated(Element))
   if(associated(Element,Oper(1)%p))then
      FPos(Oper(1)%p%ca)=FPos(Oper(1)%p%ca)+1
      exit FINDLOOP
   endif
   if((Element%CA.eq.Oper(1)%p%ca).and.(Element%Orbital.eq.Oper(1)%p%Orbital).and.(Element%Spin.eq.Oper(1)%p%Spin))then
      if(Element%has_hyb) then
         FPos(Oper(1)%p%ca)=FPos(Oper(1)%p%ca)+1
      endif
   endif
   Element=>Element%next
   enddo FINDLOOP 

   Element=>this%first      ! iterates over c^+(tau_2)
   FINDLOOP2: do while(associated(Element))
   if(associated(Element,Oper(2)%p))then
      FPos(Oper(2)%p%ca)=FPos(Oper(2)%p%ca)+1
      exit FINDLOOP2
   endif
   if((Element%CA.eq.Oper(2)%p%ca).and.(Element%Orbital.eq.Oper(2)%p%Orbital).and.(Element%Spin.eq.Oper(2)%p%Spin))then
      if(Element%has_hyb) then
         FPos(Oper(2)%p%ca)=FPos(Oper(2)%p%ca)+1
      endif
   endif
   Element=>Element%next
   enddo FINDLOOP2

end subroutine get_fpos_insert

!===============================================================================
logical function propose_remove_seg(this,DStates,Oper,N,lmax,seg,overbeta)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
   integer                             :: N
   !logical                             :: with_hyb
!output
   type(TOperPointer)                  :: Oper(N)
!local
   integer                             :: i,ipos
   type(TOper),pointer  :: elementf,elementl,element
!worm positions or vertex positions
   real(kindr) :: lmax
   logical :: seg,segoverbeta,overbeta
   integer :: nosoper,orb,spin

   propose_remove_seg=.true.
   overbeta=.false.

   !!! TODO: find out if this is a good idea ?!?
   this%outerstate_full=.false.
   this%outerstate_empty=.false.

   !!! randomly choose flavour
   orb=randint(1, DStates%NBands)
   spin=randint(1, 2)

   nosoper=this%nosoper(orb,spin)

   !!! when chosen flavour has no operators, reject
   if(nosoper.eq.0)then
      propose_remove_seg=.false.
      return
   endif

   !!! remove segment or antisegment
   if(grnd()<0.5D0)then
      seg=.true.
   else
      seg=.false.
   endif

   elementf=>this%ffirst(orb,spin)%p
   elementl=>this%flast(orb,spin)%p

   !!! check if there is a segment or antisegment at tau=0
   if(elementf%ca.eq.1)then
      segoverbeta=.false.
   else
      segoverbeta=.true.
   endif

   !!! randomly choose position of segment
   ipos=randint(1, nosoper/2)

   if(seg.eqv..true.)then
      !!! remove segment

      if(nosoper.eq.2)then
         !!! only one segment in trace

         if(segoverbeta.eqv..true.)then
            !!! this segment goes over beta

            lmax=this%beta
            oper(1)%p=>elementf
            oper(2)%p=>elementl
            !!! removing this segment over beta, sets the state at tau=0 to
            !!! empty
            this%outerstate_empty=.true.
            overbeta=.true.

         else
            !!! antisegment over beta

            lmax=this%beta
            oper(1)%p=>elementf
            oper(2)%p=>elementl
            this%outerstate_empty=.true.

         endif

         return

      endif

      if(nosoper.gt.2)then
         !!! more than 2 segments in the trace

         if(ipos.eq.nosoper/2)then
            !!! last segment

            if(segoverbeta)then
               !!! remove last segment lasting over beta

               lmax=this%beta-elementl%tau+elementf%fnext%tau
               oper(1)%p=>elementf
               oper(2)%p=>elementl
               overbeta=.true.
               return

            else
               !!! remove last normal segment

               lmax=this%beta-elementl%fprev%tau+elementf%tau
               oper(1)%p=>elementl%fprev
               oper(2)%p=>elementl
               return

            endif

         elseif(ipos.eq.1)then
            !!! first segment

            if(segoverbeta)then
               !!! remove first pair, where segoverbeta

               lmax=elementf%fnext%fnext%fnext%tau-elementf%fnext%tau
               oper(1)%p=>elementf%fnext
               oper(2)%p=>elementf%fnext%fnext
               return

            else
               !!! remove first segment, no segoverbeta

               lmax=elementf%fnext%fnext%tau-elementf%tau
               oper(1)%p=>elementf
               oper(2)%p=>elementf%fnext
               return

            endif

         else
            !!! segment in the middle

            if(segoverbeta)then

               !!! find segment
               element=>elementf%fnext
               i=1
               findloop2: do while(associated(element))
                  if(i.eq.ipos)then
                     oper(1)%p=>element
                     oper(2)%p=>element%fnext
                     lmax=element%fnext%fnext%tau-element%tau
                     return
                  endif
                  element=>element%fnext%fnext
                  i=i+1
               enddo findloop2

            else

               !!! find segment
               element=>elementf
               i=1
               findloop: do while(associated(element))
                  if(i.eq.ipos)then
                     oper(1)%p=>element
                     oper(2)%p=>element%fnext
                     lmax=element%fnext%fnext%tau-element%tau
                     return
                  endif
                  element=>element%fnext%fnext
                  i=i+1
               enddo findloop

            endif

         endif

      endif

   !!! the same for antisegments
   else

      if(nosoper.eq.2)then

         if(segoverbeta)then

            lmax=this%beta
            oper(1)%p=>elementf
            oper(2)%p=>elementl
            this%outerstate_full=.true.
            return

         else

            lmax=this%beta
            oper(1)%p=>elementf
            oper(2)%p=>elementl
            this%outerstate_full=.true.
            overbeta=.true.
            return

         endif

      endif

      if(nosoper.gt.2)then

         if(ipos.eq.nosoper/2)then

            if(segoverbeta)then

               lmax=this%beta-elementl%fprev%tau+elementf%tau
               oper(1)%p=>elementl%fprev
               oper(2)%p=>elementl
               return

            else

               lmax=this%beta-elementl%tau+elementf%fnext%tau
               oper(1)%p=>elementf
               oper(2)%p=>elementl
               overbeta=.true.
               return

            endif

         elseif(ipos.eq.1)then

            if(segoverbeta)then
            
               lmax=elementf%fnext%fnext%tau-elementf%tau
               oper(1)%p=>elementf
               oper(2)%p=>elementf%fnext
               return

            else

               lmax=elementf%fnext%fnext%fnext%tau-elementf%fnext%tau
               oper(1)%p=>elementf%fnext
               oper(2)%p=>elementf%fnext%fnext
               return

            endif

         else

            if(segoverbeta)then

               element=>elementf
               i=1
               findloop4: do while(associated(element))
               if(i.eq.ipos)then
                  oper(1)%p=>element
                  oper(2)%p=>element%fnext
                  lmax=element%fnext%fnext%tau-element%tau
                  return
               endif
               
               element=>element%fnext%fnext
               i=i+1
               enddo findloop4

            else

               element=>elementf%fnext
               i=1
               findloop3: do while(associated(element))
               if(i.eq.ipos)then
                  oper(1)%p=>element
                  oper(2)%p=>element%fnext
                  lmax=element%fnext%fnext%tau-element%tau
                  return
               endif
               
               element=>element%fnext%fnext
               i=i+1
               enddo findloop3

            endif

         endif

      endif

   endif

end function propose_remove_seg










!===============================================================================
!===============================================================================
!===============================================================================
!===============================================================================

!!!calculation of trace

!===============================================================================
!===============================================================================
!===============================================================================
!===============================================================================














!===============================================================================
!> Sets an lookup array, to know if there is a segment or antisegment at tau=0
!TODO: is this function still necessary?!
subroutine get_state_at_zero(this,oper)
!===============================================================================
   type(TTrace)                        :: this
   type(TOper),pointer                 :: Element
   type(TOperPointer)                  :: Oper(2)

   this%initial_outerstate(:,:)=this%initial_outerstate_old(:,:)

   Element=>this%first
   do while(associated(element))
      if(element%ca.eq.1)then
         this%initial_outerstate(element%orbital,element%spin)=.true.
      else
         this%initial_outerstate(element%orbital,element%spin)=.false.
      endif
      Element=>Element%next
   enddo

   if(associated(oper(1)%p))then
      if(this%outerstate_empty.eqv..true.)then
         this%initial_outerstate(oper(1)%p%orbital,oper(1)%p%spin)=.false.
         !!! DANGEROUS: when function get_state_at_zero is called twice, it
         !!! gives a right and then a wrong result, since here this variable 
         !!! is set back to .false.!
         this%outerstate_empty=.false.
      endif
      if(this%outerstate_full.eqv..true.)then
         this%initial_outerstate(oper(1)%p%orbital,oper(1)%p%spin)=.true.
         this%outerstate_full=.false.
      endif

   endif

end subroutine get_state_at_zero

!===============================================================================
!> calculate the overlap betweeen two segments
real(kindr) function two_segments_overlap(taus,taue,t1,t2)
!===============================================================================
   real(kindr) :: taus, taue, t1, t2
   real(kindr) :: tli, tri

   if(taus.lt.t1)then
      tli=t1
   else
      tli=taus
   endif

   if(taue.lt.t2)then
      tri=taue
   else
      tri=t2
   endif

   if(tli.lt.tri)then
      two_segments_overlap=tri-tli
   else
      two_segments_overlap=0d0
   endif

end function two_segments_overlap


!===============================================================================
!> calculate the overlap created by inserting a segment or antisegment
subroutine get_overlap(this,dstates,band,spin,taus,taue,overlap)
!===============================================================================
   type(TTrace)                        :: this
   type(TStates)                       :: DStates
   integer :: band,spin
   real(kindr) :: taus,taue
   real(kindr) :: overlap(dstates%nbands,2)
   integer :: ib,is
   type(TOper),pointer                 :: Element

   overlap=0d0

   do ib=1,dstates%nbands
   do is=1,2
   
      if((ib.eq.band).and.(is.eq.spin))cycle
      
      if(.not.associated(this%ffirst(ib,is)%p))then
         if(this%initial_outerstate(ib,is).eqv..false.)then
            overlap(ib,is)=0d0
         else
            overlap(ib,is)=taue-taus
         endif
      else
         if(this%ffirst(ib,is)%p%ca.eq.1)then !! segmente

            element=>this%ffirst(ib,is)%p
            do while(associated(Element))

               overlap(ib,is)=overlap(ib,is)+two_segments_overlap(taus,taue,element%tau,element%fnext%tau)
               element=>element%fnext

               if(associated(element))element=>element%fnext

            enddo

         else
            !if(www)write(*,*) "antisegmente"

            if(this%nosoper(ib,is).gt.2)then

               element=>this%ffirst(ib,is)%p
               overlap(ib,is)=overlap(ib,is)+two_segments_overlap(taus,taue,0d0,element%tau)

               element=>element%fnext
               do while(associated(Element%fnext))

                  overlap(ib,is)=overlap(ib,is)+two_segments_overlap(taus,taue,element%tau,element%fnext%tau)

                  element=>element%fnext
                  if(associated(element))element=>element%fnext

               enddo

               overlap(ib,is)=overlap(ib,is)+two_segments_overlap(taus,taue,this%flast(ib,is)%p%tau,this%beta)

            endif

            if(this%nosoper(ib,is).eq.2)then

               element=>this%ffirst(ib,is)%p
               overlap(ib,is)=overlap(ib,is)+two_segments_overlap(taus,taue,0d0,element%tau)

               overlap(ib,is)=overlap(ib,is)+two_segments_overlap(taus,taue,this%flast(ib,is)%p%tau,this%beta)

            endif

         endif

      endif

   enddo
   enddo

end subroutine get_overlap


!===============================================================================
type(TLogTr) function get_trace_seg_onestate_add(this,DStates,oper,seg,overbeta)
!===============================================================================
   type(TTrace)                              :: this
!input
   type(TStates)                             :: DStates
   type(TOperPointer)                        :: Oper(2)
!local
   integer                                   :: ib1,is1,ib2,is2
   real(kindr)                               :: tmp,mup,intp
   logical                                   :: seg,overbeta
   real(kindr), dimension(DStates%NBands, 2) :: overlap, overlap1, overlap2

   overlap1=0d0
   overlap2=0d0

   call get_state_at_zero(this,oper)
   !write(*,*) "this%initial_outerstate ", this%initial_outerstate 

   if(.not.overbeta)then

      call get_overlap(this,dstates,oper(1)%p%orbital,oper(1)%p%spin,oper(1)%p%tau,oper(2)%p%tau,overlap)

      mup=(oper(2)%p%tau-oper(1)%p%tau)*this%muimp(oper(1)%p%orbital,oper(1)%p%spin)

   else

      call get_overlap(this,dstates,oper(1)%p%orbital,oper(1)%p%spin,0d0,oper(1)%p%tau,overlap1)
      call get_overlap(this,dstates,oper(1)%p%orbital,oper(1)%p%spin,oper(2)%p%tau,this%beta,overlap2)

      overlap=overlap1+overlap2

      mup=(this%beta-oper(2)%p%tau+oper(1)%p%tau)*this%muimp(oper(1)%p%orbital,oper(1)%p%spin)

   endif

   intp=0d0
   ib2=oper(1)%p%orbital
   is2=oper(1)%p%spin
   do ib1=1,dstates%nbands
   do is1=1,2
      intp=intp+overlap(ib1,is1)*(this%u_matrix(2*(ib1-1)+is1,2*(ib2-1)+is2,2*(ib1-1)+is1,2*(ib2-1)+is2)-this%u_matrix(2*(ib1-1)+is1,2*(ib2-1)+is2,2*(ib2-1)+is2,2*(ib1-1)+is1))
   enddo
   enddo

   tmp = mup + intp
   if(seg)then
      get_trace_seg_onestate_add = TLogTr(log = -tmp, sign = 1.0_KINDR)
   else
      get_trace_seg_onestate_add = TLogTr(log = tmp, sign = 1.0_KINDR)
   endif


end function get_trace_seg_onestate_add

!===============================================================================
type(TLogTr) function get_trace_seg_onestate_rem(this,DStates,oper,seg,overbeta)
!===============================================================================
   type(TTrace)                              :: this
!input
   type(TStates)                             :: DStates
   type(TOperPointer)                        :: Oper(2)
!local
   integer                                   :: ib1,is1,ib2,is2
   real(kindr)                               :: tmp,mup,intp
   logical                                   :: seg,overbeta
   real(kindr), dimension(DStates%NBands, 2) :: overlap, overlap1, overlap2

   overlap1=0d0
   overlap2=0d0

   call get_state_at_zero(this,oper)

   if(.not.overbeta)then

      call get_overlap(this,dstates,oper(1)%p%orbital,oper(1)%p%spin,oper(1)%p%tau,oper(2)%p%tau,overlap)

      mup=(oper(2)%p%tau-oper(1)%p%tau)*this%muimp(oper(1)%p%orbital,oper(1)%p%spin)

   else

      call get_overlap(this,dstates,oper(1)%p%orbital,oper(1)%p%spin,0d0,oper(1)%p%tau,overlap1)
      call get_overlap(this,dstates,oper(1)%p%orbital,oper(1)%p%spin,oper(2)%p%tau,this%beta,overlap2)

      overlap=overlap1+overlap2

      mup=(this%beta-oper(2)%p%tau+oper(1)%p%tau)*this%muimp(oper(1)%p%orbital,oper(1)%p%spin)

   endif



   intp=0d0
   ib2=oper(1)%p%orbital
   is2=oper(1)%p%spin
   do ib1=1,dstates%nbands
   do is1=1,2
      intp=intp+overlap(ib1,is1)*(this%u_matrix(2*(ib1-1)+is1,2*(ib2-1)+is2,2*(ib1-1)+is1,2*(ib2-1)+is2)-this%u_matrix(2*(ib1-1)+is1,2*(ib2-1)+is2,2*(ib2-1)+is2,2*(ib1-1)+is1))
   enddo
   enddo

   tmp = mup + intp
   if(seg)then
      get_trace_seg_onestate_rem = TLogTr(log = tmp, sign = 1.0_KINDR)
   else
      get_trace_seg_onestate_rem = TLogTr(log = -tmp, sign = 1.0_KINDR)
   endif

end function get_trace_seg_onestate_rem

!===============================================================================
!> Generate a spin flip update.
subroutine SpinFlipUpdate_mine(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   integer                             :: i, ib
   type(TOper),pointer                 :: Element,tmp
   logical ::  iao_tmp,a1,a2
   element=>null();tmp=>null()
 
   !!! flip all spins
   Element=>this%first
   do while(associated(element))
      if(element%spin.eq.1)then
         element%spin=2
      else
         element%spin=1
      endif
      element=>element%next
   enddo

   !!! update nosoper array
   do ib=1,dstates%nbands
      i=this%nosoper(ib,1)
      this%nosoper(ib,1)=this%nosoper(ib,2)
      this%nosoper(ib,2)=i
      i = this%NOSCAOper(ib, 1, 1)
      this%NOSCAOper(ib, 1, 1) = this%NOSCAOper(ib, 2, 1)
      this%NOSCAOper(ib, 2, 1) = i
      i = this%NOSCAOper(ib, 1, 2)
      this%NOSCAOper(ib, 1, 2) = this%NOSCAOper(ib, 2, 2)
      this%NOSCAOper(ib, 2, 2) = i
   enddo

   !!! update pointers ffirst and fflast
   do ib=1,dstates%nbands
      a1=associated(this%ffirst(ib,1)%p)
      a2=associated(this%ffirst(ib,2)%p)
      if(a1.and.a2)then
         tmp=>this%ffirst(ib,1)%p
         this%ffirst(ib,1)%p=>this%ffirst(ib,2)%p
         this%ffirst(ib,2)%p=>tmp
         tmp=>this%flast(ib,1)%p
         this%flast(ib,1)%p=>this%flast(ib,2)%p
         this%flast(ib,2)%p=>tmp
      elseif(a1.and..not.a2)then
         this%ffirst(ib,2)%p=>this%ffirst(ib,1)%p
         this%ffirst(ib,1)%p=>null()
         this%flast(ib,2)%p=>this%flast(ib,1)%p
         this%flast(ib,1)%p=>null()
      elseif(a2.and..not.a1)then
         this%ffirst(ib,1)%p=>this%ffirst(ib,2)%p
         this%ffirst(ib,2)%p=>null()
         this%flast(ib,1)%p=>this%flast(ib,2)%p
         this%flast(ib,2)%p=>null()
      endif

      !!! also exchange initial_outerstate_old
      iao_tmp=this%initial_outerstate_old(ib,2)
      this%initial_outerstate_old(ib,2)=this%initial_outerstate_old(ib,1)
      this%initial_outerstate_old(ib,1)=iao_tmp

   enddo

end subroutine SpinFlipUpdate_mine

!===============================================================================
!> Generate a spin flip update.
subroutine FlavourExchangeUpdate(this,b1,s1,b2,s2)
!===============================================================================
   type(TTrace)                        :: this
!input
   integer :: b1,s1,b2,s2
!local
   integer                             :: i
   type(TOper),pointer                 :: Element,tmp
   logical ::  iao_tmp,a1,a2
   element=>null();tmp=>null()
 
   !!! exchange the flavours all spins
   Element=>this%ffirst(b1,s1)%p
   do while(associated(element))
      element%orbital=b2
      element%spin=s2
      element=>element%fnext
   enddo

   !!! exchange the flavours all spins
   Element=>this%ffirst(b2,s2)%p
   do while(associated(element))
      element%orbital=b1
      element%spin=s1
      element=>element%fnext
   enddo

   !!! update nosoper array
   i=this%nosoper(b1,s1)
   this%nosoper(b1,s1)=this%nosoper(b2,s2)
   this%nosoper(b2,s2)=i
   i = this%NOSCAOper(b1, s1, 1)
   this%NOSCAOper(b1, s1, 1) = this%NOSCAOper(b2, s2, 1)
   this%NOSCAOper(b2, s2, 1) = i
   i = this%NOSCAOper(b1, s1, 2)
   this%NOSCAOper(b1, s1, 2) = this%NOSCAOper(b2, s2, 2)
   this%NOSCAOper(b2, s2, 2) = i

   !!! update pointers ffirst and fflast
   a1=associated(this%ffirst(b1,s1)%p)
   a2=associated(this%ffirst(b2,s2)%p)
   if(a1.and.a2)then
      tmp=>this%ffirst(b1,s1)%p
      this%ffirst(b1,s1)%p=>this%ffirst(b2,s2)%p
      this%ffirst(b2,s2)%p=>tmp
      tmp=>this%flast(b1,s1)%p
      this%flast(b1,s1)%p=>this%flast(b2,s2)%p
      this%flast(b2,s2)%p=>tmp
   elseif(a1.and..not.a2)then
      this%ffirst(b2,s2)%p=>this%ffirst(b1,s1)%p
      this%ffirst(b1,s1)%p=>null()
      this%flast(b2,s2)%p=>this%flast(b1,s1)%p
      this%flast(b1,s1)%p=>null()
   elseif(a2.and..not.a1)then
      this%ffirst(b1,s1)%p=>this%ffirst(b2,s2)%p
      this%ffirst(b2,s2)%p=>null()
      this%flast(b1,s1)%p=>this%flast(b2,s2)%p
      this%flast(b2,s2)%p=>null()
   endif

   !!! also exchange initial_outerstate_old
   iao_tmp=this%initial_outerstate_old(b2,s2)
   this%initial_outerstate_old(b2,s2)=this%initial_outerstate_old(b1,s1)
   this%initial_outerstate_old(b1,s1)=iao_tmp

end subroutine FlavourExchangeUpdate

!===============================================================================
!> caution, this function has bugs!
real(kindr) function get_trace_seg_flavourexchange(this,DStates,b1,s1,b2,s2)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
   type(TOperPointer)                  :: Oper(2)
   integer :: b1,s1,b2,s2,b(2),s(2)
!local
   type(TOper),pointer                 :: Element
   integer :: ib1,is1,i
   real(kindr) :: intp(2)
   real(kindr) :: overlap(2,dstates%nbands,2),overlap1(dstates%nbands,2)
   real(kindr) :: mup(2), mu_tmp

   overlap=0d0
   overlap1=0d0

   b(1)=b1
   b(2)=b2
   s(1)=s1
   s(2)=s2

   mup=0d0
   intp=0d0

   !!!!!!!!!!!!!!!!!
   do i=1,2
   !!!!!!!!!!!!!!!!!

      call get_state_at_zero(this,oper)
      get_trace_seg_flavourexchange=0d0

      if(this%nosoper(b(i),s(i)).eq.0)then

         if(this%initial_outerstate(b(i),s(i)).eqv..false.)then

            overlap=0d0

         else

            call get_overlap(this,dstates,b(i),s(i),0d0,this%beta,overlap(i,:,:))

         endif

      else

         if(this%ffirst(b(i),s(i))%p%ca.eq.1)then !!! segmente

            element=>this%ffirst(b(i),s(i))%p
            do while(associated(element))

               call get_overlap(this,dstates,b(i),s(i),element%tau,element%fnext%tau,overlap1)
               overlap(i,:,:)=overlap(i,:,:)+overlap1

               element=>element%fnext
               if(associated(element))element=>element%fnext

            enddo

         else !!! antisegmente

            if(this%nosoper(b(i),b2).eq.2)then

               element=>this%ffirst(b(i),s(i))%p

               call get_overlap(this,dstates,b(i),s(i),0d0,element%tau,overlap1)
               overlap(i,:,:)=overlap(i,:,:)+overlap1

               call get_overlap(this,dstates,b(i),s(i),element%fnext%tau,this%beta,overlap1)
               overlap(i,:,:)=overlap(i,:,:)+overlap1

            else

               element=>this%ffirst(b(i),s(i))%p

               call get_overlap(this,dstates,b(i),s(i),0d0,element%tau,overlap1)
               overlap(i,:,:)=overlap(i,:,:)+overlap1

               element=>element%fnext
               do while(associated(element%fnext))

                  call get_overlap(this,dstates,b(i),s(i),element%tau,element%fnext%tau,overlap1)
                  overlap(i,:,:)=overlap(i,:,:)+overlap1

                  element=>element%fnext
                  if(associated(element))element=>element%fnext

               enddo

               call get_overlap(this,dstates,b(i),s(i),this%flast(b(i),s(i))%p%tau,this%beta,overlap1)
               overlap(i,:,:)=overlap(i,:,:)+overlap1

            endif

         endif

      endif

      !write(*,*) " "
      !call print_Trace_flavour(this,dstates,1,1)
      !call print_Trace_flavour(this,dstates,1,2)
      !write(*,*) "b(i)", b(i)
      !write(*,*) "s(i)", s(i)
      !write(*,*) "this%initial_outerstate", this%initial_outerstate
      !write(*,*) "overlap(i,:,:)", overlap(i,:,:)
   
      if(this%nosoper(b(i),s(i)).ne.0)then

         !!! array to count length of segments
         mup(i)=0d0
         !!! array to store beginning of segments
         mu_tmp=0d0

         !!! go through configuration and count segment lengths
         Element=>this%ffirst(b(i),s(i))%p
         do while(associated(element))
            if(element%ca.eq.1)then
               mu_tmp=element%tau
            else
               mup(i)=&
               mup(i)+element%tau-mu_tmp
            endif
            Element=>Element%fnext
         enddo

         !!! add contribution from last operator to beta
         if(this%ffirst(b(i),s(i))%p%ca.eq.2)then
            mup(i)=mup(i)+this%beta-this%flast(b(i),s(i))%p%tau
         endif

      else

         if(this%initial_outerstate(b(i),s(i)).eqv..false.)then

            mup(i)=0d0

         else

            mup(i)=this%beta

         endif

      endif

      !write(*,*) "i", i
      !write(*,*) "mup(i)", mup(i)

      !intp(i)=0d0
      do ib1=1,dstates%nbands
      do is1=1,2
         intp(i)=intp(i)+overlap(i,ib1,is1)*(this%u_matrix(2*(ib1-1)+is1,2*(b(i)-1)+s(i),2*(ib1-1)+is1,2*(b(i)-1)+s(i))-this%u_matrix(2*(ib1-1)+is1,2*(b(i)-1)+s(i),2*(b(i)-1)+s(i),2*(ib1-1)+is1))
      enddo
      enddo


   !!!!!!!!!!!!!!!!!
   enddo
   !!!!!!!!!!!!!!!!!

   get_trace_seg_flavourexchange=exp(-mup(1)-intp(1)-mup(2)-intp(2))

end function get_trace_seg_flavourexchange





!===============================================================================
!===============================================================================
!===============================================================================
!===============================================================================

!!!the routines for the calculation of the full trace

!===============================================================================
!===============================================================================
!===============================================================================
!===============================================================================




!===============================================================================
type(TLogTr) function get_trace_seg_onestate(this,DStates,oper)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
   type(TOperPointer)                  :: Oper(2)
!local
   integer                             :: ib1,is1,ib2,is2
   real(kindr)                         :: tmp,tmp2

   call get_initial_amplum(this,oper)


   call treat_mu_part_onestate(this,dstates)

   
   call treat_int_part_onestate(this,dstates)


   !!! accumulate mu
   tmp=0d0
   tmp2=0d0

   do ib1=1,dstates%nbands
   do is1=1,2
      tmp=tmp+(this%muimp(ib1,is1))*this%mu_accumulator_full(ib1,is1)
   enddo
   enddo

   !!! accumulate interaction
   do ib1=1,dstates%nbands
   do is1=1,2
   do ib2=1,dstates%nbands
   do is2=1,2
      tmp2=tmp2+(this%u_matrix(2*(ib1-1)+is1,2*(ib2-1)+is2,2*(ib1-1)+is1,2*(ib2-1)+is2))*this%int_accumulator_full(ib1,is1,ib2,is2)/2.0
      tmp2=tmp2-(this%u_matrix(2*(ib1-1)+is1,2*(ib2-1)+is2,2*(ib2-1)+is2,2*(ib1-1)+is1))*this%int_accumulator_full(ib1,is1,ib2,is2)/2.0
      !if(this%int_accumulator_full(ib1,is1,ib2,is2)/2.0.ne.0)then
         !write(*,*) "this%int_accumulator_full(ib1,is1,ib2,is2) ",ib1,is1,ib2,is2,  this%int_accumulator_full(ib1,is1,ib2,is2)
         !write(*,*) "-this%int_accumulator_full(ib1,is1,ib2,is2) ", ib1,is1,ib2,is2, this%int_accumulator_full(ib1,is1,ib2,is2)
      !endif
   enddo
   enddo
   enddo
   enddo

   !get_trace_seg_onestate=exp(-tmp-tmp2)
   !!! here i just give back the exponent of the exponential function, 
   !!! since for heavy calculations (= low temperature, systems reaches highly
   !!! excited states), this exponential function may give an overflow.
   get_trace_seg_onestate = TLogTr(log = -tmp-tmp2, sign = 1.0_KINDR)

end function get_trace_seg_onestate

!===============================================================================
subroutine get_initial_amplum(this,oper)
!===============================================================================
   type(TTrace)                        :: this
   type(TOper),pointer                 :: Element
   type(TOperPointer)                  :: Oper(2)
   !!! find out which segments overlap beta-zero, i.e. are present from zero to
   !!! first operator

   this%initial_outerstate(:,:)=this%initial_outerstate_old(:,:)
   this%has_operator(:,:)=this%initial_outerstate_old(:,:)

   Element=>this%first
   do while(associated(element))
      if(element%ca.eq.1)then
         this%initial_outerstate(element%orbital,element%spin)=.true.
         this%has_operator(element%orbital,element%spin)=.true.
      else
         this%initial_outerstate(element%orbital,element%spin)=.false.
      endif
      Element=>Element%next
   enddo


   if(associated(oper(1)%p))then

      if(this%outerstate_full.eqv..true.)then
         this%initial_outerstate(oper(1)%p%orbital,oper(1)%p%spin)=.true.
         !this%outerstate_full=.false.
      endif

      if(this%outerstate_empty.eqv..true.)then
         this%initial_outerstate(oper(1)%p%orbital,oper(1)%p%spin)=.false.
         !this%outerstate_empty=.false.
      endif

   endif


end subroutine get_initial_amplum

!===============================================================================
subroutine treat_mu_part_onestate(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   type(TOper),pointer                 :: Element
   integer :: ib1,is1

   this%mu_accumulator_full(:,:)=0d0
   this%mu_accumulator_full_tmp(:,:)=0d0
   Element=>this%first

   do while(associated(element))
      if(element%ca.eq.1)then
         this%mu_accumulator_full_tmp(element%orbital,element%spin)=element%tau
      else
         this%mu_accumulator_full(element%orbital,element%spin)=&
         this%mu_accumulator_full(element%orbital,element%spin)+element%tau-this%mu_accumulator_full_tmp(element%orbital,element%spin)
      endif
      Element=>Element%next
   enddo

   do ib1=1,dstates%nbands
   do is1=1,2
   if(this%initial_outerstate(ib1,is1).eqv..true.)then
      this%mu_accumulator_full(ib1,is1)=&
      this%mu_accumulator_full(ib1,is1)+(this%beta-this%mu_accumulator_full_tmp(ib1,is1))
   endif
   enddo
   enddo


end subroutine treat_mu_part_onestate


!===============================================================================
subroutine treat_int_part_onestate(this,DStates)
!===============================================================================
   type(TTrace)                        :: this
!input
   type(TStates)                       :: DStates
!local
   type(TOper),pointer                 :: Element
   integer :: ib1,is1,ib2,is2

   !!! create starting states:
   this%int_accumulator_full(:,:,:,:)=0d0

   this%actual_state_full=this%initial_outerstate
   this%int_accumulator_full_tmp(:,:,:,:)=0d0
   !!! treat interaction part
   Element=>this%first
   do while(associated(element))

      if(element%ca.eq.1)then
         do ib1=1,dstates%nbands
         do is1=1,2
         if(this%actual_state_full(ib1,is1))then
            this%int_accumulator_full_tmp(element%orbital,element%spin,ib1,is1)=element%tau
            this%int_accumulator_full_tmp(ib1,is1,element%orbital,element%spin)=element%tau
         endif
         enddo
         enddo

         this%actual_state_full(element%orbital,element%spin)=.true.
      endif

      if(element%ca.eq.2)then
         this%actual_state_full(element%orbital,element%spin)=.false.
         do ib1=1,dstates%nbands
         do is1=1,2
         if(this%actual_state_full(ib1,is1))then
            this%int_accumulator_full(element%orbital,element%spin,ib1,is1)=&
            this%int_accumulator_full(element%orbital,element%spin,ib1,is1)+&
            element%tau-this%int_accumulator_full_tmp(element%orbital,element%spin,ib1,is1)
            this%int_accumulator_full(ib1,is1,element%orbital,element%spin)=&
            this%int_accumulator_full(ib1,is1,element%orbital,element%spin)+&
            element%tau-this%int_accumulator_full_tmp(ib1,is1,element%orbital,element%spin)
         endif
         enddo
         enddo

      endif
         
      Element=>Element%next
   enddo

   do ib1=1,dstates%nbands
   do is1=1,2
   do ib2=1,dstates%nbands
   do is2=1,2
   if((ib1.eq.ib2).and.(is1.eq.is2))cycle
   if(this%initial_outerstate(ib1,is1).and.this%initial_outerstate(ib2,is2))then
      !!! there has to be a double-loop over band & spin, thus to avoid double
      !!! counting of the last segment, there has to be a factor 1/2.0
      this%int_accumulator_full(ib1,is1,ib2,is2)=this%int_accumulator_full(ib1,is1,ib2,is2)+(this%beta-this%int_accumulator_full_tmp(ib1,is1,ib2,is2))/2.0
      this%int_accumulator_full(ib2,is2,ib1,is1)=this%int_accumulator_full(ib2,is2,ib1,is1)+(this%beta-this%int_accumulator_full_tmp(ib2,is2,ib1,is1))/2.0
   endif
   enddo
   enddo
   enddo
   enddo

end subroutine treat_int_part_onestate


!===============================================================================

end module MTrace
!===============================================================================


#ifdef Trace_Test

!===============================================================================
!> The test program for this the trace module.
program Prog_Trace
!===============================================================================
use MParameters
use MStates
use MOperator
use MRandom
use MTrace
use MLanczos


integer,parameter       :: NHybrF=1000
type(TStates)           :: DStates
type(TTrace)            :: DTrace
real(KINDR),pointer     :: HybrFunction(:)
real(KINDR)             :: kindofmove
integer                 :: warmup,iSSt,iSt,CA,iB,iS,i,j,FPos(2),Pos(2)
real(KINDR)             :: TEmin ,Emin,tau,DetRatNew,DetRatFull,DetFullNew,DetFullOld
real(KINDR),pointer     :: Coeff(:),temp(:,:)
type(TOperPointer)      :: Oper(2)
real(KINDR),pointer     :: Q(:),R(:)
real(KINDR)             :: S,sig
logical                 :: foo
type(TSubMatrix),pointer:: fullhybr(:,:)
integer, parameter :: nftau=2
real(KINDR), allocatable :: FTau(:,:,:)
real(KINDR),allocatable                :: HEValues(:), muimp(:,:)
type(TOperator)                        :: HEVectors
type(TPsis)                            :: DPsis
type(TOperator)                        :: DHamiltonian
real(KINDR) :: spur,spur_EB
integer                             :: nsubstates
integer, allocatable                :: states2substates(:)
integer :: anz_messungen=0
!> The discretized hybridization function diagonal in both orbit and spin.
real(KINDR)                       :: test(1,2,3)
integer                    ::  nconvct=0
integer :: lsteps=0
real(KINDR) :: norm, permsig

!Variables for testing time evolution
real(KINDR), allocatable :: state_vector(:), state_vector_evolved(:), state_vector_kopie(:)
real(KINDR), allocatable :: u_matrix(:,:,:,:)
type(TPsis) :: DTransformed_Psis
integer :: actual_superstate, state_number_in_substates, integer_rep_of_state
logical :: debug=.false.
integer, dimension(1) :: sh
integer :: iii, fails=0
logical :: write_in_file, exit_loop
!################################################################################################
integer, parameter :: N_wiederholung_eines_blocks =1
integer, parameter :: N_anzahl_messwerte_fuer_mittelung=40
integer, parameter :: N_datenpunkte=15000
integer :: wiederholung_eines_blocks
integer :: datenpunkte
real :: t1,t2,t3,t4,t5
integer :: n_zeitmessung, anzahl_messwerte_fuer_mittelung
real :: zeitmessung1=0
real :: zeitmessung2=0
real :: zeitmessung3=0
integer*8 :: zeitmessung4=0
real(kindr) :: min_trace=100d0
integer :: exponenten(2,-500:500)
character (len=3) :: string
integer :: noper
logical :: global = .true.
!################################################################################################
exponenten=0

call read_ParameterFile("Parameters.in")
call init_States(DStates)
call qns2substates(DStates, nsubstates, states2substates)
call init_SubStates(DStates, nsubstates, states2substates)


! ! ! !!!!!!!!!!!!!!!!generate u-matrix; do not use; broken!!!
! ! ! call u_allocate(DStates%NBands,u_matrix)
! ! ! call u_density(get_Real_Parameter("Udd"),get_Real_Parameter("Vdd"),get_Real_Parameter("Jdd"),u_matrix)
! ! ! call u_kanamori(get_Real_Parameter("Udd"),get_Real_Parameter("Vdd"),get_Real_Parameter("Jdd"),u_matrix)

!Read in U-Matrix from file, since it now is generated in pyhton-module
!kana_2band.dat
!kana_4band.dat
!kana_3band.dat
!umatrix_density.dat
!u_matrix_coulomb_transformed.dat

call u_allocate(DStates%NBands,u_matrix)

open(4711, file='kana_4band.dat', status='old', action='read', access='sequential')
call u_readin(4711,u_matrix)
close(4711)

!allocate vectors for time-evolve-testing
allocate(state_vector(0:DStates%NStates-1))
allocate(state_vector_evolved(0:DStates%NStates-1))
allocate(state_vector_kopie(0:DStates%NStates-1))
state_vector=0
state_vector_evolved=0

call init_Psi(DPsis,DStates)
call init_Hamiltonian(DHamiltonian,DStates,DPsis)
allocate(HEValues(0:DStates%NStates-1))
call u_hamiltonian(u_matrix,DHamiltonian,DStates,DPsis)


write(stdout,*)"================================================================"
write(stdout,*)"             make hamiltonian minimal blocksize"
allocate(HEValues(0:DStates%NStates-1))
call u_hamiltonian(u_matrix,DHamiltonian,DStates,DPsis)


write(stdout,*)"================================================================"
write(stdout,*)"             make hamiltonian minimal blocksize"
write(stdout,*)"================================================================"

! start now we analyze the hamiltonian for block diagonal form
if (index(get_String_Parameter("QuantumNumbers"),"All").ne.0) then
   call analyze_hamiltonian(DHamiltonian, DStates, nsubstates, states2substates)
   call dest_Psis(DPsis)
   call dest_TOperator(DHamiltonian)
   call dest_States(DStates)
write(unit=*,fmt=*) "hier!!!!"
   call init_States(DStates)
   call init_SubStates(DStates, nsubstates, states2substates)
!  call print_SubStates(DStates)
   call init_Psi(DPsis, DStates)
   call init_Hamiltonian(DHamiltonian, DStates)
   call u_hamiltonian(u_matrix, DHamiltonian, DStates, DPsis)
   call print_Operator(DHamiltonian,DStates) 
endif
! end of analyzing hamiltonian


call init_TOperator(HEVectors,DStates) 

call diag_Operator(DHamiltonian,DStates,HEVectors,HEValues)

!write data into DStates
call set_Psis(DPsis,DStates)
call set_Hamiltonian(DHamiltonian,DStates)
call set_HEVal(HEValues,DStates)
call set_HEVec(HEVectors,DStates)

! write(unit=*,fmt=*) "print ham in occ b"
! call print_operator(DHamiltonian,DStates)

! call test_hamiltonian(DStates,HEVectors,DHamiltonian)

call transform_Psis(DStates,DTransformed_Psis)
call set_EB_Psis(DTransformed_Psis,DStates)

! call dest_Psis(DPsis)
! call dest_TOperator(DHamiltonian)
! call dest_TOperator(HEVectors)dd

allocate(FTau(DStates%NBands,2,Nftau))
ftau=0d0

!create a start configuration
allocate(Coeff(0:DStates%NStates-1))
! write(stdout,*)"beta: ", DTrace%beta

write(stdout,*)"================================================================"
write(stdout,*)"             init_Trace"
write(stdout,*)"================================================================"

call init_Trace(DTrace,DStates,ftau,nftau)


write(stdout,*)"================================================================"
write(stdout,*)"             Test get_Trace_EB"
write(stdout,*)"================================================================"

! !################################################################################################
! open(3333, file='zeitmessung_spur.dat', status='replace', action='write')
! open(3331, file='spur.dat', status='replace', action='write')
! open(3332, file='spur_eb.dat', status='replace', action='write')
! open(3334, file='delta_trace.dat', status='replace', action='write')
! 
! 
! write(unit=3333,fmt=*) "Run       OCCB     EB"
! close(3333)
! 
! !################################################################################################
! 
! exponenten=0
! ! write(unit=string,fmt="(I2.2)") noper*2
! open(3335, file="exponenten.dat", status='replace', action='write')
! 
! zeitmessung1=0
! zeitmessung2=0
! 
! !#########################################
! do datenpunkte=1,N_datenpunkte
! open(3333, file='zeitmessung_spur.dat', status='old', action='write',position='append')
! 
!  !#########################################
! 
!  call dest_Trace(DTrace)
!  call init_Trace(DTrace,DStates,ftau,nftau)
!  
!  spur=0d0
!  iii=int(rand()*20)
! !   write(unit=*,fmt=*) "iii", iii
! !   do while(spur<1d-20)
! !      write(unit=*,fmt=*) "spur", spur
! 
!     do i=1,iii
!           
!           do CA=1,2
!              allocate(Oper(CA)%p)
!              allocate(Oper(CA)%p%normr(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%norml(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%sstl(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%sstr(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%brar(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%bral(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%ketr(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%ketl(DTrace%NTruncStates))
!              Oper(CA)%p%calc_new=.true.
!              
!              allocate(Oper(CA)%p%normr_EB(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%norml_EB(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%sstl_EB(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%sstr_EB(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%brar_EB(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%bral_EB(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%ketr_EB(DTrace%NTruncStates))
!              allocate(Oper(CA)%p%ketl_EB(DTrace%NTruncStates))
!              Oper(CA)%p%calc_new_EB=.true.             
!           enddo
!           
!           do while (.not. gen_OperAdd(DTrace,DStates,Oper,FPos,Pos,2))
!           enddo
!           
!           call insert_Oper(DTrace,Oper)
! 
!     enddo
! !      spur=get_Trace(DTrace,DStates)
! !      if (spur<1d-20) then
! !         call dest_Trace(DTrace)
! !         call init_Trace(DTrace,DStates,ftau,nftau)
! ! !          write(unit=*,fmt=*) "destroy!!!"
! !      endif
! !   enddo
! write(unit=*,fmt=*) "noper", DTrace%Noper
! ! if (DTrace%Noper.ne.2*noper) then 
! ! cycle
! ! endif
! !   call print_Trace_screen(DTrace)
! 
! ! write(unit=*,fmt=*) "========================================================================================0"
! ! write(unit=*,fmt=*) "========================================================================================0"
! ! write(unit=*,fmt=*) "========================================================================================0"
! ! write(unit=*,fmt=*) "========================================================================================0"
! ! write(unit=*,fmt=*) "========================================================================================0"
! ! write(unit=*,fmt=*) "========================================================================================0"
! ! write(unit=*,fmt=*) "========================================================================================0"
! ! write(unit=*,fmt=*) "========================================================================================0"
! 
! 
! 
!  call cpu_time(t1)
!  do wiederholung_eines_blocks=1,N_wiederholung_eines_blocks 
! 
!     spur=get_Trace(DTrace,DStates)
!     
!     if (abs(spur) < min_trace .and. spur .ne. 0d0) min_trace = abs(spur)
!     
! !      write(unit=*,fmt=*) "Trace", spur
!     ! write(unit=*,fmt=*) "anzahl versuche", iii
!     ! write(unit=*,fmt=*) "==============================="
! !      stop
!     
!  enddo
! !      if (spur .ne. 0d0) write(*,*) int(log10(abs(spur)))
! 
!  if (spur .ne. 0d0) exponenten(1,int(log10(abs(spur))))=exponenten(1,int(log10(abs(spur))))+1
! !   call update_Trace(DTrace,DStates)
!  call cpu_time(t2)
!  
!  write(unit=3331,fmt=*)  DTrace%NOper, spur
!  if (datenpunkte/1000.0==datenpunkte/1000) write(unit=*,fmt=*) "datenpunkte", datenpunkte
! 
!  do wiederholung_eines_blocks=1,N_wiederholung_eines_blocks 
! 
! !      spur_EB = get_Trace_EB2(DTrace,DStates)
!     spur_EB = get_Trace_EB(DTrace,DStates)
! !      write(unit=*,fmt=*) "Trace_EB", spur_EB
! 
!  enddo
! !   if (spur_EB .ne. 0d0) write(*,*) int(log10(abs(spur_EB)))
!  if (spur_EB .ne. 0d0) exponenten(2,int(log10(abs(spur_EB))))=exponenten(2,int(log10(abs(spur_EB))))+1
! !   call update_trace_EB_neu(DTrace)
!  
! !   call update_Trace_EB(DTrace)
!  
!     write(unit=3332,fmt=*)  DTrace%NOper, spur_EB
! !      stop
! 
!  call cpu_time(t3)
!  
!  if (abs(spur-spur_EB).ge.1d-30) then
! !      write(unit=*,fmt=*) "FAIL!!!!!", fails
! !      write(unit=*,fmt=*) "spur OCCB", spur
! !      write(unit=*,fmt=*) "spur EB  ", spur_EB
! !      write(unit=*,fmt=*) "delta", abs(spur-spur_EB)
! !      write(unit=*,fmt=*) " "
!     write(unit=3334,fmt=*) fails, abs(spur-spur_EB)/abs(spur), DTrace%Noper, spur, spur_EB
!     fails = fails +1
! !      if (fails==10) stop
! !      stop
!  endif
! 
!  zeitmessung1=zeitmessung1+(t2-t1)
!  zeitmessung2 =zeitmessung2+(t3-t2)
!  !#########################################
!  write(unit=3333,fmt="(I5)",advance="no") datenpunkte
!  write(unit=3333,fmt="(A5)",advance="no") "     "
!  write(unit=3333,fmt="(F14.11)",advance="no") zeitmessung1/datenpunkte
!  write(unit=3333,fmt="(A5)",advance="no") "     "
!  write(unit=3333,fmt="(F14.11)",advance="no") zeitmessung2/datenpunkte
!  write(unit=3333,fmt="(A5)",advance="no") "     "
!  if (datenpunkte <100) then
!     write(unit=3333,fmt="(F14.11)",advance="no") 0.0
!     write(unit=3333,fmt="(A5)",advance="no") "     "
!     write(unit=3333,fmt="(F14.11)",advance="no") 0.0
!     write(unit=3333,fmt="(A5)",advance="no") "     "
!  else
!     write(unit=3333,fmt="(F14.11)",advance="no") zeitmessung1/sqrt(real(datenpunkte))/datenpunkte
!     write(unit=3333,fmt="(A5)",advance="no") "     "
!     write(unit=3333,fmt="(F14.11)",advance="no") zeitmessung2/sqrt(real(datenpunkte))/datenpunkte
!     write(unit=3333,fmt="(A5)",advance="no") "     "
!  endif
!  write(unit=3333,fmt=*) " "
!  
!  close(3333)
! 
! enddo
! 
! 
! 
! ! do datenpunkte =-500,-1
! ! 
! ! write(unit=3335,fmt=*) datenpunkte, exponenten(1, datenpunkte), exponenten(2, datenpunkte)
! ! 
! ! enddo
! 
! close(3334)
! ! close(3335)
! 
! write(unit=*,fmt=*) "min_trace", min_trace




write(stdout,*)"================================================================"
write(stdout,*)"             begin insertion and removal"
write(stdout,*)"================================================================"

open(3331, file='spur.dat', status='replace', action='write')
! open(3332, file='spur_eb.dat', status='replace', action='write')
! open(3334, file='delta_trace.dat', status='replace', action='write')
   
   
do warmup=1,10000
   kindofmove=grnd()
   write(unit=*,fmt=*) "warmup", warmup
   write(unit=*,fmt=*) "noper", DTrace%Noper
! write(unit=*,fmt=*) "hier!!!"
!===============================================================================
!insertion
!===============================================================================
   if(kindofmove<(1/(real(DTrace%Noper))))then

      do CA=1,2
         allocate(Oper(CA)%p)
!             allocate(Oper(CA)%p%normr(DTrace%NTruncStatesMax))
!             allocate(Oper(CA)%p%norml(DTrace%NTruncStatesMax))
         Oper(CA)%p%calc_new=.true.
      enddo
      
      
      if(.not. gen_OperAdd(DTrace,DStates,Oper,FPos,Pos,2))then
         cycle
      endif
      
      anz_messungen=anz_messungen+1
      call insert_Oper(DTrace,Oper)

      spur=get_Trace_EB_for_trunc(DTrace,DStates)
      !spur=get_Trace_EB(DTrace,DStates)

!     call update_Trace(DTrace)
   !  call update_trace_EB(DTrace)

   endif 
!===============================================================================
! removal
!===============================================================================
   if(kindofmove>(1/(real(DTrace%Noper)))) then
   
      if(.not.gen_OperRemove(DTrace,Oper,FPos,2))then
         cycle
      endif
      
      anz_messungen=anz_messungen+1
      
      call remove_Oper(DTrace,Oper)

       spur=get_trace_eb_for_trunc(DTrace,DStates)
     !spur=get_Trace_EB(DTrace,DStates)

!        call update_Trace(DTrace,DStates)
   !  call update_trace_EB(DTrace,DStates)

   endif
  
   call print_Trace(DTrace,3331,spur)
   
!  if (abs(spur-spur_EB).ge.1d-30) then
!     write(unit=*,fmt=*) "FAIL!!!!!", fails
!     write(unit=*,fmt=*) "spur OCCB", spur
!     write(unit=*,fmt=*) "spur EB  ", spur_EB
!     write(unit=*,fmt=*) "delta", abs(spur-spur_EB)
!     write(unit=*,fmt=*) " "
!     write(unit=3334,fmt=*) fails, abs(spur-spur_EB)/abs(spur), DTrace%Noper, spur, spur_EB, DTrace%NOper
!     fails = fails +1
!  endif
!  
!     write(unit=*,fmt=*) "NOper", DTrace%Noper
   !write(unit=3331,fmt=*) warmup, DTrace%NOper, spur
!  if (spur.ne.0d0) write(unit=3332,fmt=*) warmup, DTrace%NOper, spur_EB
   
enddo

write(unit=*,fmt=*) "Anzahl Spuren", anz_messungen

!close(3334)
close(3331)
!close(3332)

deallocate(HEValues)
call dest_Trace(DTrace)



! write_in_file=.false.


! !################################################################################################
! open(3333, file='zeitmessung.dat', status='replace', action='write')
! write(unit=3333,fmt=*) "QN:     ", trim(get_String_Parameter("QuantumNumbers"))
! write(unit=3333,fmt=*) "Run       Transformation in EB + Time-Evolve in EB     Time-Evolve in EB    Time-Evolove in OCCB"
! close(3333)
! 
! open(3333, file='zeitmessung.dat', status='old', action='write',position='append')
! !################################################################################################
! 
! !#########################################
! do datenpunkte=1,N_datenpunkte
! 
! tau=real(datenpunkte**2)*0.01
! write(unit=*,fmt=*) "tau", tau
! !#########################################
! do anzahl_messwerte_fuer_mittelung=1,N_anzahl_messwerte_fuer_mittelung
! 
! do ib=1,DStates%NBands-1
! do is=1,2
! do ca=1,2
! 
! 
! ! write(stdout,*)"================================================================"
! ! write(stdout,*)"             test time-evolution in eigenbasis"
! ! write(stdout,*)"================================================================"
! 
! do iii=0,DStates%NStates-1
!  state_vector_kopie(iii)=rand()
! enddo
! ! call ausgabe(state_vector_kopie,shape(state_vector_kopie),"state_vector_kopie",3,.false.)
! 
! do iSSt=0,DStates%NSStates-1
! !   write(unit=*,fmt=*) "ns", DStates%SubStates(iSSt)%NStates
!  state_vector_kopie(DStates%SubStates(iSSt)%Offset:DStates%SubStates(iSSt)%Offset+DStates%Substates(iSSt)%NStates-1)=&
!  &state_vector_kopie(DStates%SubStates(iSSt)%Offset:DStates%SubStates(iSSt)%Offset+DStates%Substates(iSSt)%NStates-1)/&
!  &sqrt(sum((state_vector_kopie(DStates%SubStates(iSSt)%Offset:DStates%SubStates(iSSt)%Offset+DStates%Substates(iSSt)%NStates-1))**2))
! enddo
! 
! call cpu_time(t1)
! 
! 
! 
! !#########################################
! do wiederholung_eines_blocks=1,N_wiederholung_eines_blocks
! do integer_rep_of_state=0,DStates%NStates-1
! 
!  if (write_in_file) open(12345, file='ausgabe_trace_m.txt', status='replace', action='write')
!  if (write_in_file)write(12345,fmt=*) "asdf hier test"
! 
! !   write(unit=*,fmt=*) "integer_rep_of_state", integer_rep_of_state
! !   state_number_in_substates=DStates%StatesSStates(integer_rep_of_state,3)
! !   state_vector=0
! !   state_vector(state_number_in_substates)=1
!  state_vector=state_vector_kopie
!  actual_superstate=DStates%StatesSStates(integer_rep_of_state,2)
! !   write(unit=*,fmt=*) "actual_superstate", actual_superstate
!  state_vector_evolved=0
! 
! 
!  ! call ausgabe(DStates%StatesSStates(:,2),shape(DStates%StatesSStates(:,2)),"DStates%StatesSStates(:,2)",3,.false.)
!  ! call ausgabe(DStates%StatesSStates(:,1),shape(DStates%StatesSStates(:,1)),"DStates%StatesSStates(:,1)",3,.false.)
! 
!  call transform_state_occbasis_to_eb(HEVectors,DStates,state_vector,.false.)
! 
!  call time_evolve_in_eigenbasis(DStates,state_vector,state_vector_evolved,actual_superstate,tau)
!     if (DStates%SubStates(actual_superstate)%Connect(ib,is,ca)==-1) then
! !      write(unit=*,fmt=*) "not connected#!!!!!"
!  else
!  call apply_psi_to_ket(state_vector,actual_superstate,ib,is,ca,DStates,DTransformed_Psis,.false.)
!  endif
!  
!  call transform_state_eb_to_occbasis(HEVectors,DStates,state_vector_evolved,.false.)
! 
!  !============Begin Vektor Ausgabe
!  sh=shape(state_vector_evolved)
!  if (write_in_file) write(12345,"(A)",advance="no") "state_vector_evolved, "
!  if (write_in_file) write(12345,"(A11)",advance="no") "Dimension: "
!  if (write_in_file) write(12345,"(I3)") sh(1)
! 
!  do iii=0,sh(1)-1
!  if (write_in_file) write(12345,"(I5,F30.25,A8)",advance="no")  iii, state_vector_evolved(iii), " "
!  if (write_in_file) write(12345,fmt=*) " "
!  enddo
!  !============End Vektor Ausgabe
!  if (write_in_file) close(12345)
! 
! enddo
! enddo
! !   write(unit=*,fmt=*) "huhu!!!!"
! 
! call cpu_time(t2)
! 
! !#########################################
! do wiederholung_eines_blocks=1,N_wiederholung_eines_blocks
! do integer_rep_of_state=0,DStates%NStates-1
! !   state_number_in_substates=DStates%StatesSStates(integer_rep_of_state,3)
! !   state_vector=0
! !   state_vector(state_number_in_substates)=1
!  state_vector=state_vector_kopie
!  actual_superstate=DStates%StatesSStates(integer_rep_of_state,2)
! !   do i = 1,10
! !      write(unit=*,fmt=*) i
! !   enddo
! !      
!  actual_superstate=DStates%StatesSStates(integer_rep_of_state,2)
!  call transform_state_occbasis_to_eb(HEVectors,DStates,state_vector,.false.)
!  call transform_state_eb_to_occbasis(HEVectors,DStates,state_vector,.false.)
! enddo
! enddo
! 
! 
! !   write(stdout,*)"================================================================"
! !   write(stdout,*)"             test time-evolution in occupation basis"
! !   write(stdout,*)"================================================================"
!  
! call cpu_time(t3)
! 
! !#########################################
! do wiederholung_eines_blocks=1,N_wiederholung_eines_blocks
! do integer_rep_of_state=0,DStates%NStates-1
! 
!  if (write_in_file) open(123, file='ausgabe_trace_l.csv', status='replace', action='write')
! 
! !   write(unit=*,fmt=*) "integer_rep_of_state", integer_rep_of_state
! !   state_number_in_substates=DStates%StatesSStates(integer_rep_of_state,3)
! !   state_vector=0
! !   state_vector(state_number_in_substates)=1
!  state_vector=state_vector_kopie
!  actual_superstate=DStates%StatesSStates(integer_rep_of_state,2)
! !   write(unit=*,fmt=*) "actual_superstate", actual_superstate
!  state_vector_evolved=0
! 
! 
!  ! call ausgabe(state_vector,shape(state_vector),"state_vector",3,.false.)
! 
! !   write(unit=*,fmt=*) "blockkgroesse", DStates%Substates(actual_superstate)%NStates
!  
!  call LanczosTimeEvolve(DStates,actual_superstate,&
!        state_vector(DStates%SubStates(actual_superstate)%Offset:DStates%SubStates(actual_superstate)%Offset+DStates%SubStates(actual_superstate)%NStates-1),&
!        state_vector_evolved(DStates%SubStates(actual_superstate)%Offset:DStates%SubStates(actual_superstate)%Offset+DStates%SubStates(actual_superstate)%NStates-1),&
!        tau,norm,DTrace%Egs,nconvct,lsteps)
!  
!  zeitmessung4=zeitmessung4+lsteps
! !   write(unit=*,fmt=*) "lsteps", lsteps
!     
! 
!  if (DStates%SubStates(actual_superstate)%Connect(ib,is,ca)==-1) then
! !      write(unit=*,fmt=*) "not connected#!!!!!"
!  else
!      call SMatVecProd(DStates%SubStates(actual_superstate)%Psis(ib,is,ca),&
!     state_vector_evolved(DStates%SubStates(actual_superstate)%Offset:DStates%SubStates(actual_superstate)%Offset+DStates%SubStates(actual_superstate)%NStates-1),&
!     state_vector(DStates%SubStates(DStates%SubStates(actual_superstate)%Connect(ib,is,ca))%Offset:&
!     DStates%SubStates(DStates%SubStates(actual_superstate)%Connect(ib,is,ca))%Offset+DStates%SubStates(DStates%SubStates(actual_superstate)%Connect(ib,is,ca))%NStates-1))
!  endif
! 
!  !write(unit=*,fmt=*) DStates%SubStates(actual_superstate)%Psis(1,1,1)%dat 
!  
!   
! 
!     
!  !============Begin Vektor Ausgabe
!  sh=shape(state_vector_evolved)
!  if (write_in_file) write(123,"(A)",advance="no") "state_vector_evolved, "
!  if (write_in_file) write(123,"(A11)",advance="no") "Dimension: "
!  if (write_in_file) write(123,"(I3)") sh(1)
! 
!  do iii=0,sh(1)-1
!  if (write_in_file) write(123,"(I5,F30.25,A8)",advance="no")  iii, state_vector_evolved(iii), " "
!  if (write_in_file) write(123,fmt=*) " "
!  enddo
!  !============End Vektor Ausgabe
!  if (write_in_file) close(123)
! 
!  !     write(unit=*,fmt=*) "lsteps", lsteps
!  !     write(unit=*,fmt=*) "nconvct", nconvct
!        
!  ! call ausgabe(state_vector_evolved,shape(state_vector_evolved),"state_vector_evolved",3,.false.)
! 
! enddo
! ! write(unit=*,fmt=*) "zeitmessung4", zeitmessung4
! enddo
! 
! call cpu_time(t4)
! 
! 
! !#########################################
! ! zeitmessung_array(1,anzahl_messwerte_fuer_mittelung)= t2-t1
! ! zeitmessung_array(2,anzahl_messwerte_fuer_mittelung)= t2-t1-(t3-t2)
! ! zeitmessung_array(3,anzahl_messwerte_fuer_mittelung)= t4-t3
! ! zeitmessung_array(4,anzahl_messwerte_fuer_mittelung)= real(lsteps)
! 
! zeitmessung1=zeitmessung1+t2-t1
! zeitmessung2=zeitmessung2+t2-t1-(t3-t2)
! zeitmessung3=zeitmessung3+t4-t3
! ! write(unit=*,fmt=*) "zeitmessung4", zeitmessung4
! 
! enddo
! enddo
! enddo
! 
! 
! enddo
! 
! !#########################################
! write(unit=*,fmt=*) t2-t1,t2-t1-(t3-t2),t4-t3
! write(unit=3333,fmt="(F14.11)",advance="no") tau
! write(unit=3333,fmt="(A5)",advance="no") "     "
! write(unit=3333,fmt="(F14.11)",advance="no") zeitmessung1/real(N_anzahl_messwerte_fuer_mittelung)
! write(unit=3333,fmt="(A5)",advance="no") "     "
! write(unit=3333,fmt="(F14.11)",advance="no") zeitmessung2/real(N_anzahl_messwerte_fuer_mittelung)
! write(unit=3333,fmt="(A5)",advance="no") "     "
! write(unit=3333,fmt="(F14.11)",advance="no") zeitmessung3/real(N_anzahl_messwerte_fuer_mittelung)
! write(unit=3333,fmt="(A5)",advance="no") "     "
! write(unit=3333,fmt="(F14.11)") real(zeitmessung4)/real(N_anzahl_messwerte_fuer_mittelung)/N_wiederholung_eines_blocks/(DStates%NStates-1)/(DStates%NBands-1)/2/2
! 
! zeitmessung1=0
! zeitmessung2=0
! zeitmessung3=0
! zeitmessung4=0
! 
! enddo


end program Prog_Trace

#endif
