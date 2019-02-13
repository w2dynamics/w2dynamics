!!! Module which performs a continuous time QMC simulation
! Working
!===============================================================================
module MCTQMC
!===============================================================================
use MParameters
use MOperator
use MTrace
use MStates
use MRandom
use MAngularMomentum
use LegendrePoly, only: m1pow
use type_progress
use signals
use iso_c_binding
use MAusgabe
use ft_worm
#ifdef USE_NFFT
 use ft_fast, ft_beta => beta
#else
 use ft_naive, ft_beta => beta
#endif
implicit none

#ifdef USE_REALTIME
    interface
       subroutine wasted(secs, nanos) bind(C)
           use iso_c_binding
           integer(c_long), intent(out) :: secs, nanos
       end subroutine
    end interface
#endif

   integer(C_INT64_T)         :: Nwarmups,Nmeas
   integer                    :: NGtau,NCorr,NGiw,NBands,NStates
   integer(c_int64_t)         :: NSlide
   integer                    :: N4tau, N4leg, N4iwb, N4iwf, N2iwb, N3iwb, N3iwf
   real(KINDR),allocatable    :: Gtau(:,:,:)   ! o,sp,tau
   real(KINDR),allocatable    :: Gtau_mean_step(:,:,:), Gtau_mid_step(:,:,:) ! o,sp,iNmeas
   real(KINDR),allocatable    :: sign_step(:) ! iNMeas
   integer(c_int64_t)         :: g_inmeas
   real(KINDR),allocatable    :: Gtau_full(:,:,:,:,:)   ! o,sp,tau
   real(KINDR),allocatable    :: ntau_n0(:,:,:,:,:)   ! o,sp,tau
   real(KINDR),allocatable    :: nn_disc(:)   ! o,sp,tau
   real(KINDR),allocatable    :: G4tau(:,:,:,:,:,:,:) !o1,s1,o2,s2,t14,t24,t34
   !real(KINDR),allocatable    :: GGtau(:,:,:,:,:,:)
   complex(KINDC),allocatable :: Giw(:,:,:)
   real(KINDR),allocatable    :: Giw_lookup_tau(:,:,:)
   complex(KINDC),allocatable :: Giw_lookup_M(:,:,:)
   integer,allocatable        :: Giw_lookup_row(:,:)
   integer                    :: NLookup_nfft
   complex(KINDC),allocatable :: GSigmaiw(:,:,:) !o,s,iw
   complex(KINDC),allocatable :: G2iw(:,:,:,:)  !o,s,iw,iw'
   complex(KINDC),allocatable :: G4iw(:,:,:,:,:,:,:) !o1,s1,o2,s2,iw,iw',iW
   complex(KINDC),allocatable :: G4iw_pp(:,:,:,:,:,:,:) !o1,s1,o2,s2,iw,iw',iW
   real(KINDR)                :: meas_sign
   real(KINDR)                :: beta
   real(KINDR)                :: PercentageGlobalMove, PercentageTauShiftMove
   real(KINDR)                :: PercentageOuterMove
   real(KINDR)                :: Percentage4OperatorMove
   real(KINDR),allocatable    :: globalmove_check(:)
   integer(c_int64_t), allocatable :: AccPair(:,:,:)   ! o,sp,tau
   real(KINDR), allocatable   :: AccPairTau(:)
   integer                    :: NAccPair
   real(KINDR)                :: AccPairMax
   integer                    :: apt_index
   !matsubara array for the trivial fourier transform
   real(KINDR),allocatable    :: iw(:)
   real(KINDR),allocatable    :: iwb2(:), iwb3(:), iwf3(:)
 

   !worm quantities
   complex(KINDC),allocatable :: Giw_worm(:) ! iw
   real(KINDR),allocatable    :: Gtau_worm(:) ! tau
   complex(KINDC),allocatable :: GSigmaiw_worm(:) !iw
  
   !worm quantities (component sampling)
   complex(KINDC),allocatable :: G4iw_worm(:,:,:) ! iv,iv',iw
   complex(KINDC),allocatable :: H4iw_worm(:,:,:) ! iv,iv',iw
   complex(KINDC),allocatable :: P2iw_worm(:) ! iw
   complex(KINDC),allocatable :: P2iwPP_worm(:) ! iw
   complex(KINDC),allocatable :: P3iw_worm(:,:) ! iv, iw
   complex(KINDC),allocatable :: P3iwPP_worm(:,:) ! iv, iw

   real(KINDR), allocatable   :: P2tau_worm(:) ! tau 
   real(KINDR), allocatable   :: P2tauPP_worm(:) ! tau

   !imaginary-time buffer for buffered nfft
   real(c_double),allocatable            :: tau_worm(:)     ! (tau-argumens)*buffer_size
   complex(c_double_complex),allocatable :: val_worm(:)         ! buffer_size
   integer(c_int64_t)                    :: tau_fill

   !number of worm operators in specific sector
   integer, parameter, dimension(2:9)   :: NOperWorm = (/2,4,4,6,4,4,4,4/)

   !logical to control when to Fourier transform Worm operators
   !we need to store previous worm configuration 
   logical                    :: isNew(2:9)
   !wormTau stores time differences for fermionic and bosonic fourier transform for each sector
   real(KINDR)                :: wormSum(2:9),wormTau(2:9,3)
   integer                    :: wormO(2:9,4),wormS(2:9,4)
    
   !we make the umatrix global to qmc
   real(KINDR),allocatable    :: u_matrix(:,:,:,:)
   !outer/dangling indices (spin, orbital, general) of umatrix for improved estimator
   integer                    :: u_s, u_o, u_g

   !parameter eta to control how much greens function space contributes
   !spin-orbit eta matrix with combined spin-orbit index, order c,c^dag (c,c^dag)
   real(KINDR)                :: PercentageWormInsert,PercentageWormReplace

   !we consider 9 different spaces, which are defined in the sector variable
   !in the CTQMCsimulate subroutine
   !CntSampling: steps in each space taken
   !CntMeas: measurements in each space made
   integer(c_int64_t)                  :: CntSampling(9), CntMeas(9)
   
   real(KINDR)                :: mean_sign
   ! state dimension index of the form (sst offset + in-sst index)
   integer,allocatable        :: StatesSuperstates(:)          ! superstate of (eigen-)state
   real(KINDR),allocatable    :: EigenstatesEnergies(:)        ! energy of eigenstate
   integer,allocatable        :: OccbasisMapping(:, :, :)      ! (dm-index, orb, spin) -> occ-eigval
   real(KINDR),allocatable    :: Histo(:,:,:)
   real(KINDR),allocatable    :: OuterSuperstateHisto(:)
   real(KINDR),allocatable    :: OuterStateHisto(:)
   real(KINDR),allocatable    :: SignHistoSuperstates(:)
   real(KINDR),allocatable    :: SignHistoStates(:)
   real(KINDR),allocatable    :: TraceContribSuperstates(:)
   real(KINDR),allocatable    :: TraceContribStates(:)
   real(KINDR),allocatable    :: Histo_seg(:)
   logical, allocatable       :: nf(:,:)
!> the Lanczos histogram; steps needed for convergence in the time evolution
   real(KINDR), allocatable   :: lhisto(:)
!> the non-reversibility histogram; additional steps needed for reversible
!! trace
   real(KINDR), allocatable   :: rhisto(:)
   real(KINDR),allocatable    :: occ(:,:,:,:)  ! o1,s1,o2,s2
   real(KINDR),allocatable    :: rho2(:,:,:,:)  ! 4 contracted orbitals and spins
   real(KINDR),allocatable    :: double_occ(:,:,:,:) ! o1, s1, o2, s2  ??? we have to exchange double_occ and single_occ in the old matrix functions
   real(KINDR),allocatable    :: single_occ(:,:) ! o1, s1 
   real(KINDR),allocatable    :: rho1(:,:)  ! 2 contracted orbitals and spins
   real(KINDR),allocatable    :: DensityMatrix(:,:)
   real(KINDR),allocatable    :: ExpResDensityMatrix(:,:,:,:,:)
   
   !we count the acceptance of hyb pairs in each sector including z sector
   real(KINDR)                :: AccRem(9),AccAdd(9),AccGlob,AccShift
   integer(c_int64_t)         :: TryRem(9),TryAdd(9),TryGlob,TryShift
   integer(c_int64_t)         :: AccQNRem(9),AccQNAdd(9),AccQNGlob

   ! Separate acceptance counters for outer and inner moves (try =
   ! move attempts, acc = absolute number of accepted moves, accqn =
   ! absolute number of moves rejected due to quantum number checking)
   ! These counters are set but code for post-processing or writing
   ! them is not present by default.
   real(KINDR)                :: AccAddInner,AccAddOuter,AccRemInner,AccRemOuter
   integer(c_int64_t)         :: TryAddInner,TryAddOuter,TryRemInner,TryRemOuter
   integer(c_int64_t)         :: AccQNAddInner,AccQNAddOuter,AccQNRemInner,AccQNRemOuter

   real(KINDR)                :: AccRem4,AccAdd4
   integer(C_INT64_T)         :: TryRem4,TryAdd4
   real(KINDR)                :: AccFlavc
   integer(C_INT64_T)         :: TryFlavc
   
   !worm accepts and tries
   !1 -> 1P Worm, 2-> 2P Worm, 3-> IE Sigma, 4->IE chi 
   !(mind shift of 1 as we dont consider Z)
   real(KINDR)                :: AccWormRem(2:9),AccWormAdd(2:9),AccWormRep(2:9)
   integer(c_int64_t)         :: TryWormRem(2:9),TryWormAdd(2:9),TryWormRep(2:9)
   integer(c_int64_t)         :: AccQNWormAdd(2:9),AccQNWormRem(2:9)
   
   !FIXME: this needs to generalized for components
   real(KINDR)                :: wormEta(2:9)
   
   !array for nfft
   real(c_double),allocatable                :: taus(:)
   complex(c_double_complex),allocatable     :: vals(:)
   complex(c_double_complex),allocatable     :: matsubaras(:)
   integer(c_int),allocatable                :: d(:)

   !sign measurement counter in z
   integer(c_int64_t)         :: cnt_sign_z
   integer                    :: FourPnt,compTotal
   
   integer                    :: NLegTau,NLegMax
   real(KINDR),allocatable    :: LegendreGrid(:,:)
   real(KINDR),allocatable    :: GLeg(:,:,:)  ! o,sp,l
   real(KINDR),allocatable    :: GLeg_full(:,:,:,:,:)  ! o,sp,l
   complex(KINDC),allocatable :: G4leg(:,:,:,:,:,:,:) !o1,s1,o2,s2,l,l´,n

   integer                    :: timings_maxorder
   real(KINDR), allocatable   :: timings_giw(:), timings_g4iw_ft(:), timings_g4iw_add(:)
   integer, allocatable       :: count_giw(:), count_g4iw(:)

   real(KINDR)                :: time_calibration, time_warmup, time_sim, time_sim_steps

   logical                    :: fancy_prog = .false.

! opaquely handling the pointers to derived types; see A. Pletzer et al. 
   integer(C_INT64_T)         :: ipDStates(12)
   integer(C_INT64_T)         :: ipDTrace(12)
!   type(TStates)              :: DStates
!   type(TTrace)               :: DTrace

   integer, parameter :: GF4_IMAGTIME = 1, GF4_LEGENDRE = 2, GF4_MATSUBARA = 4, GF4_WORM = 8
   integer            :: WormphConv,ZphConv
   
   logical :: b_Eigenbasis, b_Densitymatrix,b_meas_susz,b_segment
   logical :: b_meas_susz_mat
   logical :: b_Giw_lookup
   logical :: b_statesampling
   logical :: b_offdiag
   logical :: b_full_offdiag
   logical :: b_exch

   logical :: GtauDetRat

   ! communicates to the caller if the QMC sampling was aborted by a signal.
   !  0 = no abort, 1 = aborted, some data, 2 = aborted, no data
   integer                    :: aborted = 0
   integer                    :: SimID = -1
   character(4)               :: simstr

   interface
      real(KIND=KIND(0.D0)) function DDOT(N, DX, INCX, DY, INCY)
         integer                             :: N, INCX, INCY
         real(KIND=KIND(0.D0)), dimension(*) :: DX, DY
      end function DDOT
   end interface
contains

#ifndef USE_REALTIME
   subroutine wasted(secs, nanos)
      use iso_c_binding
      integer(c_long), intent(out) :: secs, nanos

      secs = 0
      nanos = 0
      stop '[wasted] Not implemented'
   end subroutine
#endif

!===============================================================================
subroutine init_CTQMC()
!===============================================================================
   use ft_common
!input
   type(TStates),pointer         :: DStates
   type(TStates_pointer)         :: pDStates
   integer                       :: i
   real(kindr), parameter        :: pi=dacos(-1d0)
   real(kindr)                   :: dt
!output

   pDStates=transfer(ipDStates,pDStates)
   DStates => pDStates%ptr

   ! segment or matrix-vector
   if (get_Integer_Parameter("segment")==0) then
      b_segment = .false.
      b_Densitymatrix = .true.   !!! FIXME: when segment, then it is not possible to measure the dm, since the beta-half state will not exist
      write(0,*) "--> Using Matrix solver."
   else
      b_segment = .true.
      b_Densitymatrix = .false.
      write(0,*) "--> Using Segment solver."
   endif

   if (.not. b_segment) then
      if (get_Integer_Parameter("statesampling") /= 0) then
         b_statesampling = .true.
      else
         b_statesampling = .false.
      end if
   else
      b_statesampling = .false.
   end if

   Nwarmups=get_LongInteger_Parameter("Nwarmups")
   Nmeas=get_LongInteger_Parameter("Nmeas")
   NGtau=get_Integer_Parameter("Ntau")
   NBands=get_Integer_Parameter("Nd")
   NStates=4**NBands
   NCorr=get_Integer_Parameter("NCorr")
   NGiw=get_Integer_Parameter("Niw")
   NLookup_nfft=get_Integer_Parameter("NLookup_nfft")
!   NBin=get_Integer_Parameter("NBin")
   N4leg=get_Integer_Parameter("N4leg")
   N4tau=get_Integer_Parameter("N4tau")
   N4iwb=get_Integer_Parameter("N4iwb")
   N4iwf=get_Integer_Parameter("N4iwf")

   N2iwb=get_Integer_Parameter("N2iwb")
   N3iwb=get_Integer_Parameter("N3iwb")
   N3iwf=get_Integer_Parameter("N3iwf")

!  NLegTau=get_Integer_Parameter("NLegTau")
   NLegMax=get_Integer_Parameter("NLegMax")
   beta=get_Real_Parameter("beta")
   !FourPnt=get_Integer_Parameter("FourPnt")
   PercentageGlobalMove=get_Real_Parameter("PercentageGlobalMove")
   if (b_segment) then
      PercentageTauShiftMove = 0d0
   else
      PercentageTauShiftMove = get_Real_Parameter("PercentageTauShiftMove")
   end if
   PercentageOuterMove=get_Real_Parameter("PercentageOuterMove")
   Percentage4OperatorMove=get_Real_Parameter("Percentage4OperatorMove")
   
   PercentageWormInsert=get_Real_Parameter("PercentageWormInsert")
   PercentageWormReplace=get_Real_Parameter("PercentageWormReplace")
   
   ! two-particl worm measurement
   if( IAND(FourPnt,GF4_WORM) /= 0 .and. get_integer_parameter("WormMeasG4iw") /= 0) then
      if(get_Integer_Parameter("WormPHConvention") /= 0) then
         WormphConv = 1
      else
         WormphConv = 0
      endif
   
      allocate(G4iw_worm(-N4iwf:N4iwf-1,-N4iwf:N4iwf-1,-N4iwb:N4iwb))
      G4iw_worm = cmplx(0, kind=KINDC)
   endif

   if( IAND(FourPnt,GF4_WORM) /= 0 .and. get_integer_parameter("WormMeasH4iw") /= 0) then
      if(get_Integer_Parameter("WormPHConvention") /= 0) then
         WormphConv = 1
      else
         WormphConv = 0
      endif
      allocate(H4iw_worm(-N4iwf:N4iwf-1,-N4iwf:N4iwf-1,-N4iwb:N4iwb))
      H4iw_worm = cmplx(0, kind=KINDC)
   endif
   
   wormEta=get_Real_Parameter("WormEta")
   
   !! choose Eigenbasis or Krylov solver
   !! we cannot use Krylov for segment algorithm
   !if ((get_Integer_Parameter("Eigenbasis")==0) .and. (b_segment.eqv..false.)) then
       !b_Eigenbasis = .false.
       !b_Densitymatrix = get_Integer_Parameter("MeasDensityMatrix")
       !!write(0,*) "--> Using Krylov solver."
    !else
       !b_Eigenbasis = .true.
       !b_Densitymatrix = .true.
       !!write(0,*) "--> Using Eigenbasis solver."
   !endif
 
   if (get_Integer_Parameter("Eigenbasis")==0) then
      stop "You requested Krylov, but this no longer supported."
   else
      b_Eigenbasis = .true.
      b_Densitymatrix = .true.
      write(0,*) "--> Using eigenbasis of impurity."
   endif

   ! choose offdiagonal or diagonal hybridisation
   if (get_Integer_Parameter("offdiag")==0) then
      b_offdiag = .false.
      write(0,*) "--> Using diagonal F(tau)."
   else
      b_offdiag = .true.
      write(0,*) "--> Using full offdiag F(tau)."
   endif

   ! in case of offdiagonal calculation, choose if hybridisation matrix is
   ! computed always from scratch for debugging, or by fast update algorithm
   if (get_Integer_Parameter("full_offdiag")==0) then
      b_full_offdiag = .false.
      write(0,*) "--> fast offdiag."
   else
      b_full_offdiag = .true.
      write(0,*) "--> full offdiag."
   endif

   if (get_Integer_Parameter("flavourchange_moves") == 0) then
      b_exch=.false.
   else
      b_exch = .true.
   end if

   if(allocated(Gtau))then
     write(*,*)"Gtau already allocated"
   endif
   allocate(Gtau(NBands,2,NGtau))
   if (get_Integer_Parameter("Gtau_mean_step") /= 0) allocate(Gtau_mean_step(NBands,2,Nmeas))
   if (get_Integer_Parameter("Gtau_mid_step") /= 0) allocate(Gtau_mid_step(NBands,2,Nmeas))
   if (get_Integer_Parameter("sign_step") /= 0) allocate(sign_step(Nmeas))
   allocate(Gtau_full(NBands,2,NBands,2,NGtau))
   allocate(ntau_n0(NBands,2,NBands,2,0:NGtau))

   dt=beta/float(ngtau)
   allocate(nn_disc(0:NGtau))
   do i=0,ngtau
   nn_disc(i)=i*dt
   enddo

   !do i=1,ngtau
   !write(*,*) "i,nn_disc(i)", i,nn_disc(i)
   !enddo
   
   if(get_integer_parameter("MeasGtauDetRat") /= 0) then
      GtauDetRat=.true.
   else
      GtauDetRat=.false.
   endif

   b_meas_susz = .false.
   b_meas_susz_mat = .false.
   if (get_Integer_Parameter("MeasSusz").eq.1) then
      if (b_segment) then
         write(*,*) "--> meas susceptibility (segment)!"
         b_meas_susz=.true.
      else
         b_meas_susz_mat = .true.
      end if
   endif

   ! Distinguish between the modes of measuring the 4-point GF
   if( IAND(FourPnt,GF4_IMAGTIME) /= 0 ) then
      allocate(G4tau(NBands, 2, NBands, 2, N4tau, N4tau, N4tau))
      G4tau = 0
   endif
   if( IAND(FourPnt,GF4_LEGENDRE) /= 0 ) then
      allocate(G4leg(NBands, 2, NBands, 2, 0:N4tau-1, 0:N4tau-1, -N4iwb:N4iwb))
      G4leg = 0
   endif
   if( IAND(FourPnt,GF4_MATSUBARA) /= 0 ) then
      allocate(G4iw(NBands, 2, NBands, 2, -N4iwf:N4iwf-1, -N4iwf:N4iwf-1, -N4iwb:N4iwb))
      G4iw = 0
   endif
   if( IAND(FourPnt,GF4_MATSUBARA) /= 0 .and. get_integer_parameter("MeasG4iwPP") /= 0) then
      allocate(G4iw_pp(NBands, 2, NBands, 2, -N4iwf:N4iwf-1, -N4iwf:N4iwf-1, -N4iwb:N4iwb))
      G4iw_pp = 0
   endif
   
   if(FourPnt /= 0) then
      if(get_Integer_Parameter("ZPHConvention") /= 0) then
         ZphConv = 1
      else
         ZphConv = 0
      endif
   endif

   if( get_integer_parameter("MeasGiw") /= 0) then
      allocate(Giw(NBands,2,-NGiw:NGiw-1))
      Giw = cmplx(0, kind=KINDC)
   endif

   if (allocated(Giw) .and. NLookup_nfft > 0) then
      allocate(Giw_lookup_tau(NLookup_nfft, NBands, 2))
      Giw_lookup_tau = 0_KINDR

      allocate(Giw_lookup_M(NLookup_nfft, NBands, 2))
      Giw_lookup_M = cmplx(0d0, 0d0, kind=KINDC)

      allocate(Giw_lookup_row(NBands, 2))
      Giw_lookup_row = 1

      b_Giw_lookup = .true.
#ifndef USE_NFFT
      write (*, *) "WARNING: Your settings for Giw measurement require compilation with NFFT support, Giw will NOT be measured!"
#endif
   else
      NLookup_nfft = 0
      b_Giw_lookup = .false.
   end if
   
   if( get_integer_parameter("WormMeasGiw") /= 0) then
      allocate(Giw_worm(-NGiw:NGiw-1))
      Giw_worm = cmplx(0, kind=KINDC)
      
      if(.not.allocated(iw)) then
         allocate(iw(-NGiw:NGiw-1))   
         do i=-NGiw,NGiw-1,1
            iw(i)=2.d0*pi/beta*(dble(i)+0.5d0)
         end do
      endif
   endif
   
   if( get_integer_parameter("WormMeasGtau") /= 0) then
      allocate(Gtau_worm(NGtau))
      Gtau_worm = 0
   endif
   
   if( get_integer_parameter("MeasGSigmaiw") /= 0) then
      allocate(GSigmaiw(NBands,2,-NGiw:NGiw-1))
      GSigmaiw = cmplx(0, kind=KINDC)
   endif
   
   if( get_integer_parameter("WormMeasGSigmaiw") /=0) then
      allocate(GSigmaiw_worm(-NGiw:NGiw-1))
      GSigmaiw_worm = cmplx(0d0, kind=KINDC)

      if(.not.allocated(iw)) then
         allocate(iw(-NGiw:NGiw-1))   
         do i=-NGiw,NGiw-1,1
            iw(i)=2.d0*pi/beta*(dble(i)+0.5d0)
         end do
      endif
   endif

   if( get_integer_parameter("WormMeasP2iwPH") /= 0) then
      allocate(P2iw_worm(-N2iwb:N2iwb))
      P2iw_worm = cmplx(0, kind=KINDC)
      
      if(.not.allocated(iwb2)) then
         allocate(iwb2(-N2iwb:N2iwb))   
         do i=-N2iwb,N2iwb,1
            iwb2(i)=(2.d0*pi/beta)*dble(i)
         end do
      endif
   endif


   if( get_integer_parameter("WormMeasP2iwPP") /= 0) then
      allocate(P2iwPP_worm(-N2iwb:N2iwb))
      P2iwPP_worm = cmplx(0, kind=KINDC)
      
      if(.not.allocated(iwb2)) then
         allocate(iwb2(-N2iwb:N2iwb))   
         do i=-N2iwb,N2iwb,1
            iwb2(i)=(2.d0*pi/beta)*dble(i)
         end do
      endif
   endif

   if( get_integer_parameter("WormMeasP2tauPH") /= 0) then
      allocate(P2tau_worm(NGtau))
      P2tau_worm = 0d0
   endif

   if( get_integer_parameter("WormMeasP2tauPP") /= 0) then
      allocate(P2tauPP_worm(NGtau))
      P2tauPP_worm = 0d0
   endif

   if( get_integer_parameter("WormMeasP3iwPH") /= 0) then
      allocate(P3iw_worm(-N3iwf:N3iwf-1,-N3iwb:N3iwb))
      P3iw_worm = cmplx(0, kind=KINDC)
      
      if(.not.allocated(iwb3)) then
         allocate(iwb3(-N3iwb:N3iwb))   
         do i=-N3iwb,N3iwb,1
            iwb3(i)=(2.d0*pi/beta)*dble(i)
         end do
      endif

      if(.not.allocated(iwf3)) then
         allocate(iwf3(-N3iwf:N3iwf-1))   
         do i=-N3iwf,N3iwf-1,1
            iwf3(i)=(2.d0*pi/beta)*(dble(i)+0.5d0)
         end do
      endif
   endif

   if( get_integer_parameter("WormMeasP3iwPP") /= 0) then
      allocate(P3iwPP_worm(-N3iwf:N3iwf-1,-N3iwb:N3iwb))
      P3iwPP_worm = cmplx(0, kind=KINDC)
      
      if(.not.allocated(iwb3)) then
         allocate(iwb3(-N3iwb:N3iwb))   
         do i=-N3iwb,N3iwb,1
            iwb3(i)=(2.d0*pi/beta)*dble(i)
         end do
      endif

      if(.not.allocated(iwf3)) then
         allocate(iwf3(-N3iwf:N3iwf-1))   
         do i=-N3iwf,N3iwf-1,1
            iwf3(i)=(2.d0*pi/beta)*(dble(i)+0.5d0)
         end do
      endif
   endif
   
   if(allocated(occ))then
     write(*,*)"occ already allocated"
   endif
   
   allocate(occ(NBands, 2,NBands, 2))
   allocate(rho2(2*NBands,2*NBands,2*NBands,2*NBands))
   allocate(double_occ(NBands, 2, NBands, 2))
   allocate(single_occ(NBands, 2))
   allocate(rho1(2*NBands,2*NBands))
   
   if(allocated(densitymatrix))then
     write(*,*)"densitymatrix already allocated"
   endif
   
   if(b_Densitymatrix .eqv. .true.) &
    allocate(DensityMatrix(NStates,NStates))
    
   if(get_Integer_Parameter("MeasExpResDensityMatrix").ne.0)&
     allocate(ExpResDensityMatrix(NBands,2,0:get_Integer_Parameter("MaxHisto"),NStates,NStates))

   if(allocated(gleg))then
     write(*,*)"gleg already allocated"
   endif
   allocate(GLeg(NBands,2,NLegMax))
   allocate(GLeg_full(NBands,2,NBands,2,NLegMax))

   allocate(StatesSuperstates(0:DStates%NStates-1))
   allocate(EigenstatesEnergies(0:DStates%NStates-1))
   allocate(OccbasisMapping(0:DStates%NStates-1, NBands, 2))
   if(allocated(histo))then
     write(*,*)"histo already allcoated"
   endif
   allocate(Histo(NBands,2,0:get_Integer_Parameter("MaxHisto")))
   allocate(OuterSuperstateHisto(0:DStates%NSStates-1))
   allocate(OuterStateHisto(0:DStates%NStates-1))
   allocate(SignHistoSuperstates(0:DStates%NSStates-1))
   allocate(SignHistoStates(0:DStates%NStates-1))
   allocate(TraceContribSuperstates(0:DStates%NSStates-1))
   allocate(TraceContribStates(0:DStates%NStates-1))
   allocate(Histo_seg(0:NBands*2))
   allocate(nf(NBands,2))
   allocate(lhisto(PMAX))
   allocate(rhisto(PMAX))

   allocate(globalmove_check(NBands*2))
   globalmove_check(:)=0d0

   allocate(AccPair(NBands,2,NGtau))
   NAccPair = NGtau
   AccPairMax = beta

   meas_sign = 0.0_KINDR
   cnt_sign_z = 0
   Gtau=0d0
   if (allocated(Gtau_mean_step)) Gtau_mean_step=0d0
   if (allocated(Gtau_mid_step)) Gtau_mid_step=0d0
   if (allocated(sign_step)) sign_step=0d0
   Gtau_full=0d0
   ntau_n0=0d0
   GLeg=0d0
   GLeg_full=0d0
   Histo=0d0
   OuterSuperstateHisto=0d0
   OuterStateHisto = 0_KINDR
   SignHistoSuperstates = 0_KINDR
   SignHistoStates = 0_KINDR
   TraceContribSuperstates=0d0
   TraceContribStates=0d0
   histo_seg=0d0
   occ=0d0
   rho2=0d0
   double_occ=0d0
   single_occ=0d0
   rho1=0d0
   
   call init_counters()
   
   if(allocated(DensityMatrix)) DensityMatrix=0d0
   if(allocated(ExpResDensityMatrix)) ExpResDensityMatrix=0d0
   
   
   lhisto = 0d0
   rhisto = 0d0
         
   call ft_common_init(2*NGiw, 2*(N4iwf+N4iwb), beta)

   timings_maxorder = get_integer_parameter("MeasurementTiming")
   if(timings_maxorder >= 0) then
       allocate(timings_giw(0:timings_maxorder))
       allocate(count_giw(0:timings_maxorder))
       timings_giw = 0
       count_giw = 0
       allocate(timings_g4iw_ft(0:timings_maxorder))
       allocate(timings_g4iw_add(0:timings_maxorder))
       allocate(count_g4iw(0:timings_maxorder))
       timings_g4iw_ft = 0
       timings_g4iw_add = 0
       count_g4iw = 0
   endif

   !initialize sector to 1, i.e. start sampling in partition function space
   time_calibration = 0.0_KINDR
   time_warmup = 0_KINDR
   time_sim = 0_KINDR
   time_sim_steps = 0_KINDR

   if(get_Integer_Parameter("MeasG2iw") /= 0) then
       allocate(G2iw(NBands,2,2*(N4iwf+N4iwb),2*(N4iwf+N4iwb)))
       G2iw = 0
   endif

   aborted = 0
end subroutine init_CTQMC

!===============================================================================
subroutine init_counters()
!===============================================================================
   AccRem=0
   AccAdd=0
   AccPair=0
   AccGlob=0
   AccShift=0
   AccAddInner=0
   AccAddOuter=0
   AccRemInner=0
   AccRemOuter=0
   AccWormAdd=0
   AccWormRem=0
   AccWormRep=0
   AccQNAdd=0
   AccQNRem=0
   AccQNGlob=0
   AccQNAddInner=0
   AccQNAddOuter=0
   AccQNRemInner=0
   AccQNRemOuter=0
   AccQNWormAdd=0
   AccQNWormRem=0
   TryRem=0
   TryAdd=0
   TryGlob=0
   TryShift=0
   TryAddInner=0
   TryAddOuter=0
   TryRemInner=0
   TryRemOuter=0
   TryWormAdd=0
   TryWormRem=0
   TryWormRep=0
   TryFlavc=0
   AccFlavc=0

   apt_index = 1
   CntSampling = 0
   CntMeas = 0

   isNew = .false.
   wormSum = 0d0
   wormTau = 0d0

end subroutine init_counters
!===============================================================================


!===============================================================================
subroutine dest_CTQMC()
!===============================================================================
   type(TStates), pointer        :: DStates
   type(TTrace), pointer         :: DTrace
   type(TStates_pointer)         :: pDStates
   type(TTrace_pointer)          :: pDTrace

   pDStates = transfer(ipDStates, pDStates)
   pDTrace = transfer(ipDTrace, pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr
   call ft_cleanup()
   call dest_Trace(DTrace)
   call dest_States(DStates)
   deallocate(pDStates%ptr)
   deallocate(pDTrace%ptr)

   if(allocated(Gtau))deallocate(Gtau)
   if(allocated(Gtau_mean_step))deallocate(Gtau_mean_step)
   if(allocated(Gtau_mid_step))deallocate(Gtau_mid_step)
   if(allocated(sign_step))deallocate(sign_step)
   if(allocated(Gtau_full))deallocate(Gtau_full)
   if(allocated(nn_disc))deallocate(nn_disc)
   if(allocated(ntau_n0))deallocate(ntau_n0)
   if(allocated(StatesSuperstates))deallocate(StatesSuperstates)
   if(allocated(EigenstatesEnergies))deallocate(EigenstatesEnergies)
   if(allocated(OccbasisMapping))deallocate(OccbasisMapping)
   if(allocated(Histo))deallocate(Histo)
   if(allocated(OuterSuperstateHisto))deallocate(OuterSuperstateHisto)
   if(allocated(OuterStateHisto))deallocate(OuterStateHisto)
   if(allocated(SignHistoSuperstates))deallocate(SignHistoSuperstates)
   if(allocated(SignHistoStates))deallocate(SignHistoStates)
   if(allocated(TraceContribSuperstates))deallocate(TraceContribSuperstates)
   if(allocated(TraceContribStates))deallocate(TraceContribStates)
   if(allocated(globalmove_check))deallocate(globalmove_check)
   if(allocated(AccPair))deallocate(AccPair)
   if(allocated(AccPairTau))deallocate(AccPairTau)
   if(allocated(histo_seg))deallocate(histo_seg)
   if(allocated(nf))deallocate(nf)
   if(allocated(Occ))deallocate(Occ)
   if(allocated(double_occ))deallocate(double_occ)
   if(allocated(single_occ))deallocate(single_occ)
   if(allocated(rho1))deallocate(rho1)
   if(allocated(rho2))deallocate(rho2)
   if(allocated(DensityMatrix))deallocate(DensityMatrix)
   if(allocated(ExpResDensityMatrix))deallocate(ExpResDensityMatrix)
!   if(allocated(LegendreGrid))deallocate(LegendreGrid)
   if(allocated(GLeg))deallocate(GLeg)
   if(allocated(GLeg_full))deallocate(GLeg_full)
   if(allocated(iw)) deallocate(iw)
   if(allocated(iwb2)) deallocate(iwb2)
   if(allocated(iwb3)) deallocate(iwb3)
   if(allocated(iwf3)) deallocate(iwf3)

   if(allocated(G4tau)) deallocate(G4tau)
   if(allocated(G4leg)) deallocate(G4leg)
   if(allocated(G4iw)) deallocate(G4iw)
   if(allocated(G4iw_pp)) deallocate(G4iw_pp)
   if(allocated(G2iw)) deallocate(G2iw)

   if(allocated(Giw)) deallocate(Giw)
   if(allocated(Giw_worm)) deallocate(Giw_worm)
   if(allocated(Gtau_worm)) deallocate(Gtau_worm)
   if(allocated(GSigmaiw)) deallocate(GSigmaiw)
   if(allocated(GSigmaiw_worm)) deallocate(GSigmaiw_worm)
   if (allocated(Giw_lookup_tau)) deallocate(Giw_lookup_tau)
   if (allocated(Giw_lookup_M)) deallocate(Giw_lookup_M)
   if (allocated(Giw_lookup_row)) deallocate(Giw_lookup_row)

   if(allocated(P2iw_worm)) deallocate(P2iw_worm)
   if(allocated(P2iwPP_worm)) deallocate(P2iwPP_worm)
   if(allocated(P2tau_worm)) deallocate(P2tau_worm)
   if(allocated(P2tauPP_worm)) deallocate(P2tauPP_worm)
   if(allocated(P3iw_worm)) deallocate(P3iw_worm)
   if(allocated(P3iwPP_worm)) deallocate(P3iwPP_worm)

   if(allocated(timings_giw)) deallocate(timings_giw)
   if(allocated(timings_g4iw_add)) deallocate(timings_g4iw_add)
   if(allocated(timings_g4iw_ft)) deallocate(timings_g4iw_ft)
   if(allocated(count_giw)) deallocate(count_giw)
   if(allocated(count_g4iw)) deallocate(count_g4iw)
   if(allocated(lhisto)) deallocate(lhisto)
   if(allocated(rhisto)) deallocate(rhisto)
   
   if (allocated(u_matrix)) deallocate(u_matrix)
   if(allocated(g4iw_worm)) deallocate(g4iw_worm)
   if(allocated(h4iw_worm)) deallocate(h4iw_worm)
   if(allocated(tau_worm)) deallocate(tau_worm)  
   if(allocated(val_worm)) deallocate(val_worm)
 
end subroutine dest_CTQMC

!===============================================================================

subroutine measure_giw_nfft()
   use ft_common, only: nfreq
   logical, parameter :: debug = .false.
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace

   integer :: iorb, ispin
   integer(c_long) :: start_secs, start_nanos, end_secs, end_nanos
   complex(KINDC) :: this_giw(nfreq,2)

   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   if (timings_maxorder >= DTrace%noper) then
      call wasted(start_secs, start_nanos)
   endif

   do iorb = 1, NBands
      do ispin = 1, 2
         if(allocated(GSigmaiw)) then
            ! Measure improved estimators
            if(debug) write(0,"(I1,'/',I1,100ES10.2)") iorb, ispin, DTrace%urho(iorb,ispin)%v
            this_giw = get_giw( &
                size(DTrace%Minv(iorb,ispin)%Mat,1), DTrace%Minv(iorb,ispin)%Mat, &
                DTrace%tau_c(iorb,ispin)%v, DTrace%tau_a(iorb,ispin)%v, &
                m_prefactors=DTrace%urho(iorb,ispin)%v)
            GSigmaiw(iorb,ispin,:) = GSigmaiw(iorb,ispin,:) + this_giw(:,2)
         else
            this_giw = get_giw( &
                size(DTrace%Minv(iorb,ispin)%Mat,1), DTrace%Minv(iorb,ispin)%Mat, &
                DTrace%tau_c(iorb,ispin)%v, DTrace%tau_a(iorb,ispin)%v)
         endif
         Giw(iorb,ispin,:) = Giw(iorb,ispin,:) + this_giw(:,1)
      enddo
   enddo

   if (timings_maxorder >= DTrace%noper) then
      call wasted(end_secs, end_nanos)
      timings_giw(DTrace%noper) = timings_giw(DTrace%noper) + &
                    (end_secs - start_secs + (end_nanos - start_nanos)*1.0e-9)
      count_giw(DTrace%noper) = count_giw(DTrace%noper) + 1
!      write (0,'(I5,F10.6,2X,F10.6)') &
!          DTrace%noper, end_secs - start_secs + (end_nanos - start_nanos)*1.0e-9, &
!          timings_giw(DTrace%noper)/count_giw(DTrace%noper)
   endif
end subroutine

!TODO: combine all iw worm measurements and all tau worm measurements into two 
!      worm measurement routines
!measurement of giw_worm in greens function space
subroutine measure_giw_worm()
   type(TTrace), pointer       :: DTrace
   type(TTrace_pointer)        :: pDTrace
   
   real(KINDR)                :: prefact
   
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   prefact=1d0/wormEta(2)
   
   !we only fourier transform new worm configs, old ones are just accumulated
   if(isNew(2).eqv..true.) then
   
      if(tau_fill .eq. size(val_worm)) then

         allocate(matsubaras(4*NGiw),d(1))
         d(1)=size(matsubaras)

         call ft_nd(tau_worm,val_worm,matsubaras,d)
         giw_worm(:)=giw_worm(:)+matsubaras(2::2)

         deallocate(matsubaras,d)
         tau_fill=0
      endif

      tau_worm(tau_fill+1) = (DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(2)%p%tau)/(2d0*DTrace%beta)
      val_worm(tau_fill+1) = DTrace%cfgsign*prefact
      tau_fill = tau_fill + 1
      isNew(2)=.false.

   else
      val_worm(tau_fill) = val_worm(tau_fill)+DTrace%cfgsign*prefact
   endif

end subroutine

!measurement of gtau_worm in greens function space
subroutine measure_gtau_worm()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace
   
   integer                    :: itau
   real(KINDR)                :: tau, sgn, prefact
   
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr
   
   tau=DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(2)%p%tau
   sgn=1d0
   
   if(tau<0.d0) then
      !antiperiodicity includes time ordering
      tau=tau+DTrace%beta
      sgn=-sgn
   endif
   
   sgn=sgn*DTrace%cfgsign
   
   itau=int(dble(NGtau)*tau/DTrace%beta)+1
   if(itau.gt.NGtau) itau=NGtau
   

   prefact=1d0/wormEta(2)
   
   Gtau_worm(itau)=Gtau_worm(itau)+sgn*prefact 

end subroutine

!measurement of gsigma in worm space using improved estimators
subroutine measure_ie_sigma_worm()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace
   
   integer                    :: s(4), o(4), g(4), i
   real(KINDR)                :: prefact

   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   prefact = 1d0/wormEta(3)   

   do i=1,4
      o(i)=DTrace%wormContainer(i)%p%Orbital
      s(i)=DTrace%wormContainer(i)%p%Spin
      g(i)=2*(o(i)-1)+s(i)
   enddo
      
   !u_s,u_o and u_g set to dangling index
   prefact = prefact*sign(1d0,0.5d0*(u_matrix(u_g,g(2),g(1),g(3))-u_matrix(g(2),u_g,g(1),g(3))))

   !we only fourier transform new worm configs, old ones are just accumulated
   if(isNew(3).eqv..true.) then      

      if(tau_fill .eq. size(val_worm)) then

         allocate(matsubaras(4*NGiw),d(1))
         d(1)=size(matsubaras)

         call ft_nd(tau_worm,val_worm,matsubaras,d)
         GSigmaiw_worm(:)=GSigmaiw_worm(:)+matsubaras(2::2)

         deallocate(matsubaras,d)
         tau_fill=0
      endif

      tau_worm(tau_fill+1) = (DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(4)%p%tau)/(2d0*DTrace%beta)
      val_worm(tau_fill+1) = DTrace%cfgsign*prefact
      tau_fill = tau_fill + 1
      isNew(3)=.false.

   else
      val_worm(tau_fill) = val_worm(tau_fill)+DTrace%cfgsign*prefact
   endif

end subroutine

!measurement of p2iw_worm
subroutine measure_p2iw_worm()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace
   
   real(KINDR)                :: prefact
   
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   prefact=1d0/wormEta(6)
   
   !we only fourier transform new worm configs, old ones are just accumulated
   if(isNew(6).eqv..true.) then
   
      if(tau_fill .eq. size(val_worm)) then

         allocate(matsubaras(4*N2iwb+3),d(1))
         d(1)=size(matsubaras)
   
         call ft_nd(tau_worm,val_worm,matsubaras,d)
         p2iw_worm(:)=p2iw_worm(:)+matsubaras(2::2)
       
         deallocate(matsubaras,d)
         tau_fill=0
      endif

      tau_worm(tau_fill+1) = (DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(3)%p%tau)/(2d0*DTrace%beta)
      val_worm(tau_fill+1) = DTrace%cfgsign*prefact
      tau_fill = tau_fill + 1
      isNew(6)=.false.
      
   else
      val_worm(tau_fill) = val_worm(tau_fill)+DTrace%cfgsign*prefact
   endif

end subroutine

!measurement of p2tau_worm
subroutine measure_p2tau_worm()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace
   
   integer                    :: itau
   real(KINDR)                :: tau, sgn, prefact
   
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr
   
   tau=DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(3)%p%tau
   sgn=1d0
   
   if(tau<0.d0) then
      tau=tau+DTrace%beta
      !does not pick up a sign due to pairs being exchanged
      sgn=sgn
   endif
   
   sgn=sgn*DTrace%cfgsign
   
   itau=int(dble(NGtau)*tau/DTrace%beta)+1
   if(itau.gt.NGtau) itau=NGtau
   
   prefact=1d0/wormEta(6)
   
   p2tau_worm(itau)=p2tau_worm(itau)+sgn*prefact 

end subroutine

!measurement of p2iwpp_worm
subroutine measure_p2iwpp_worm()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace

   real(KINDR)                :: prefact
   
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   prefact=1d0/wormEta(7)
   
   !we only fourier transform new worm configs, old ones are just accumulated
   if(isNew(7).eqv..true.) then
  
      if(tau_fill .eq. size(val_worm)) then
         allocate(matsubaras(4*N2iwb+3),d(1))
         d(1)=size(matsubaras)

         call ft_nd(tau_worm,val_worm,matsubaras,d)
         p2iwpp_worm(:) = p2iwpp_worm(:)+matsubaras(2::2)

         deallocate(matsubaras,d)
         tau_fill=0
      endif

      tau_worm(tau_fill+1) = (DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(2)%p%tau)/(2d0*DTrace%beta)
      val_worm(tau_fill+1) = DTrace%cfgsign*prefact
      tau_fill = tau_fill + 1
      isNew(7)=.false.
  
   else
      val_worm(tau_fill) = val_worm(tau_fill)+DTrace%cfgsign*prefact
   endif

end subroutine

!measurement of p2taupp_worm
subroutine measure_p2taupp_worm()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace
   
   integer                    :: itau
   real(KINDR)                :: tau, sgn, prefact
   
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr
   
   tau=DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(2)%p%tau
   sgn=1d0
   
   if(tau<0.d0) then
      tau=tau+DTrace%beta
      !does not pick up a sign due to pairs being exchanged
      sgn=sgn
   endif
     
   sgn=sgn*DTrace%cfgsign   
   itau=int(dble(NGtau)*tau/DTrace%beta)+1
   if(itau.gt.NGtau) itau=NGtau
     
   prefact=1d0/wormEta(7)
  
   p2taupp_worm(itau)=p2taupp_worm(itau)+sgn*prefact 

end subroutine

!measurement of p3iw_worm
subroutine measure_p3iw_worm()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace
   
   integer                    :: i,j,k
   real(KINDR)                :: prefact
   
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   prefact=1d0/wormEta(8)
   
   !we only fourier transform new worm configs, old ones are just accumulated
   if(isNew(8).eqv..true.) then
   
      !fourier transform previous worm
      if(tau_fill .eq. size(val_worm)) then

         allocate(matsubaras((4*N3iwf)*(4*N3iwb+3)),d(2))
         d(1)=4*N3iwf
         d(2)=4*N3iwb+3
         
         call ft_nd(tau_worm,val_worm,matsubaras,d)
      
         do i=-N3iwf,N3iwf-1
            do j=-N3iwb,N3iwb
                k=(4*N3iwb+3)*(2*(i+N3iwf)+1)+(2*(j+N3iwb)+1)+1
                p3iw_worm(i,j)=p3iw_worm(i,j)+matsubaras(k)
            enddo
         enddo

         deallocate(matsubaras,d)
         tau_fill = 0

      endif
      
      tau_worm(2*tau_fill+1) = (DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(2)%p%tau)/(2d0*DTrace%beta)
      tau_worm(2*tau_fill+2) = (DTrace%wormContainer(2)%p%tau-DTrace%wormContainer(3)%p%tau)/(2d0*DTrace%beta)

      val_worm(tau_fill+1)=DTrace%cfgsign*prefact
      tau_fill = tau_fill + 1 
      isNew(8)=.false.
      
   else
      val_worm(tau_fill) = val_worm(tau_fill)+DTrace%cfgsign*prefact
   endif

end subroutine


!measurement of p3iwpp_worm
subroutine measure_p3iwpp_worm()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace
   
   integer                    :: i,j,k
   real(KINDR)                :: prefact
   
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   prefact=1d0/wormEta(9)
   
   !we only fourier transform new worm configs, old ones are just accumulated
   if(isNew(9).eqv..true.) then
   
      !fourier transform previous worm
      if(tau_fill .eq. size(val_worm)) then

         allocate(matsubaras((4*N3iwf)*(4*N3iwb+3)),d(2))
         d(1)=4*N3iwf
         d(2)=4*N3iwb+3
     
         call ft_nd(tau_worm,val_worm,matsubaras,d)
 
         do i=-N3iwf,N3iwf-1
            do j=-N3iwb,N3iwb
                k=(4*N3iwb+3)*(2*(i+N3iwf)+1)+(2*(j+N3iwb)+1)+1
                p3iwpp_worm(i,j)=p3iwpp_worm(i,j)+matsubaras(k)
            enddo
         enddo

         deallocate(matsubaras,d)
         tau_fill = 0

      endif

      tau_worm(2*tau_fill+1)=(DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(3)%p%tau)/(2d0*DTrace%beta)
      tau_worm(2*tau_fill+2)=(DTrace%wormContainer(3)%p%tau-DTrace%wormContainer(2)%p%tau)/(2d0*DTrace%beta)

      val_worm(tau_fill+1)=DTrace%cfgsign*prefact
      tau_fill = tau_fill + 1 
      isNew(9)=.false.
      
   else
      val_worm(tau_fill) = val_worm(tau_fill)+DTrace%cfgsign*prefact
   endif

end subroutine

subroutine measure_g4iw_worm()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace
   integer                    :: in1,in2,in3,k
   real(KINDR)                :: prefact

   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr
   
   prefact=1d0/wormEta(4)
   
   !we only fourier transform new worm configs, old ones are just accumulated
   if(isNew(4).eqv..true.) then
   
      !fourier transform if buffer is full
      if(tau_fill .eq. size(val_worm)) then

         allocate(matsubaras((4*N4iwf)*(4*N4iwf)*(4*N4iwb+3)),d(3))
         d(1)=4*N4iwf
         d(2)=4*N4iwf
         d(3)=4*N4iwb+3
         
         call ft_nd(tau_worm,val_worm,matsubaras,d)

         do in3=-N4iwb,N4iwb
             do in2=-N4iwf,N4iwf-1
                do in1=-N4iwf,N4iwf-1

                   k=(4*N4iwf)*(4*N4iwb+3)*(2*(in1+N4iwf)+1)+(4*N4iwb+3)*(2*(in2+N4iwf)+1)+(2*(in3+N4iwb)+1)+1
                   g4iw_worm(in1,in2,in3)=g4iw_worm(in1,in2,in3) + matsubaras(k)

               enddo
            enddo
         enddo

         deallocate(matsubaras,d)
         tau_fill = 0

      endif

      tau_worm(3*tau_fill+1) = (DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(2)%p%tau)/(2d0*DTrace%beta)
      tau_worm(3*tau_fill+2) = (DTrace%wormContainer(3)%p%tau-DTrace%wormContainer(4)%p%tau)/(2d0*DTrace%beta)
      if(WormphConv.eq.0) then
         tau_worm(3*tau_fill+3)= (DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(4)%p%tau)/(2d0*DTrace%beta)
      else
         tau_worm(3*tau_fill+3)= (DTrace%wormContainer(2)%p%tau-DTrace%wormContainer(3)%p%tau)/(2d0*DTrace%beta)
      endif

      val_worm(tau_fill+1)=DTrace%cfgsign*prefact
      tau_fill = tau_fill + 1
      isNew(4)=.false.
      
   else
      val_worm(tau_fill) = val_worm(tau_fill)+DTrace%cfgsign*prefact
   endif
   
end subroutine

!measurement of chi in worm space using improved estimators
subroutine measure_ie_chi_worm()
   
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace
   integer                    :: in1,in2,in3,k,i
   integer                    :: s(6), o(6), g(6)
   real(KINDR)                :: prefact

   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr
   
   do i=1,6
      o(i)=DTrace%wormContainer(i)%p%Orbital
      s(i)=DTrace%wormContainer(i)%p%Spin
      g(i)=2*(o(i)-1)+s(i)
   enddo
   
   !u_s,u_o and u_g set to dangling index
   prefact = sign(1d0,0.5d0*(u_matrix(u_g,g(2),g(1),g(3))-u_matrix(g(2),u_g,g(1),g(3))))
        
   prefact=prefact/wormEta(5)
   
   !we only fourier transform new worm configs, old ones are just accumulated
   if(isNew(5).eqv..true.) then
   
      !fourier transform previous worm
      if(tau_fill .eq. size(val_worm)) then
         
         allocate(matsubaras((4*N4iwf)*(4*N4iwf)*(4*N4iwb+3)),d(3))
         d(1)=4*N4iwf
         d(2)=4*N4iwf
         d(3)=4*N4iwb+3
         
         call ft_nd(tau_worm,val_worm,matsubaras,d)

         do in1=-N4iwf,N4iwf-1
            do in2=-N4iwf,N4iwf-1
               do in3=-N4iwb,N4iwb

                   k=(4*N4iwf)*(4*N4iwb+3)*(2*(in1+N4iwf)+1)+(4*N4iwb+3)*(2*(in2+N4iwf)+1)+(2*(in3+N4iwb)+1)+1
                   h4iw_worm(in1,in2,in3)=h4iw_worm(in1,in2,in3) + matsubaras(k)

               enddo
            enddo
         enddo

         deallocate(matsubaras,d)
         tau_fill = 0

      endif

      tau_worm(3*tau_fill+1)=(DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(4)%p%tau)/(2d0*DTrace%beta)
      tau_worm(3*tau_fill+2)=(DTrace%wormContainer(5)%p%tau-DTrace%wormContainer(6)%p%tau)/(2d0*DTrace%beta)
      if(WormphConv.eq.0) then
         tau_worm(3*tau_fill+3)=(DTrace%wormContainer(1)%p%tau-DTrace%wormContainer(6)%p%tau)/(2d0*DTrace%beta)
      else
         tau_worm(3*tau_fill+3)=(DTrace%wormContainer(4)%p%tau-DTrace%wormContainer(5)%p%tau)/(2d0*DTrace%beta)
      endif


      val_worm(tau_fill+1)=DTrace%cfgsign*prefact
      tau_fill = tau_fill + 1
      isNew(5)=.false.
      
   else
      val_worm(tau_fill)=val_worm(tau_fill)+DTrace%cfgsign*prefact
   endif
   
end subroutine

subroutine measure_g4iw_nfft()
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace

   integer :: o1, sp1, o2, sp2, w1, w2, wb
   integer(c_long) :: start_secs, start_nanos, end_secs, end_nanos

   complex(KINDC) :: this_g2iw(NBands,2,-n2freq/2:n2freq/2-1,-n2freq/2:n2freq/2-1)

   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   if (timings_maxorder >= DTrace%noper) then
      call wasted(start_secs, start_nanos)
   endif

   ! MEASURE_GIW_NFFT MUST BE CALLED BEFORE THAT
   ! Do the two-frequency Fourier-transform beforehand
   do o1 = 1, NBands; do sp1 = 1, 2
      this_g2iw(o1,sp1,:,:) = get_g2iw(size(DTrace%Minv(o1,sp1)%Mat,1), &
           DTrace%Minv(o1,sp1)%Mat, &
           DTrace%tau_c(o1,sp1)%v, DTrace%tau_a(o1,sp1)%v)
   enddo; enddo

   if (timings_maxorder >= DTrace%noper) then
      call  wasted(end_secs, end_nanos)
      timings_g4iw_ft(DTrace%noper) = timings_g4iw_ft(DTrace%noper) + &
                    (end_secs - start_secs + (end_nanos - start_nanos)*1.0e-9)
      call wasted(start_secs, start_nanos)
   endif

   if (allocated(g2iw)) g2iw = g2iw + this_g2iw

   ! Scaling (2b)^2 * nw^2 * nwb -> OK!
   ! At the level of a single measurement, the two-particle Green's function
   ! (here in the particle-hole frequency convention) is just the product of
   ! two single-particle Green's functions, but depending on two frequencies
   ! each because only the mean makes it tau-translation invariant.
   !
   !      A -->---+-----+-->--- A      /  A --->--- A    A --  ->- A  \
   !        w1    | chi | w1+wb       /                      \/        \
   !              |     |          =  \                -     /\        /
   !      B -->---+-----+-->--- B      \  B --->--- B    A --  ->- A  / MC
   !        w2+wb         w2
   !
   ! This follows by Fourier-transformation of the Wick-like formula for the
   ! measurement by removal in tau (found, e.g., in Christian's thesis).
   !
   !FIXME: replace foralls with for and column major layout
   if(allocated(g4iw)) then
      if(ZphConv.eq.0) then
         forall(w1=-N4iwf:N4iwf-1, w2=-N4iwf:N4iwf-1, wb=-N4iwb:N4iwb)
            forall(o1=1:NBands, sp1=1:2, o2=1:NBands, sp2=1:2)
               g4iw(o1,sp1,o2,sp2,w1,w2,wb) = g4iw(o1,sp1,o2,sp2,w1,w2,wb) + &
                     this_g2iw(o1,sp1,w1+wb,w1)*this_g2iw(o2,sp2,w2,w2+wb)
            end forall
            ! Crossing term in case both orbitals are equal
            forall(o1=1:NBands, sp1=1:2)
               g4iw(o1,sp1,o1,sp1,w1,w2,wb) = g4iw(o1,sp1,o1,sp1,w1,w2,wb) - &
                     this_g2iw(o1,sp1,w2,w1)*this_g2iw(o1,sp1,w1+wb,w2+wb)
            end forall
         end forall
      else
         forall(w1=-N4iwf:N4iwf-1, w2=-N4iwf:N4iwf-1, wb=-N4iwb:N4iwb)
            forall(o1=1:NBands, sp1=1:2, o2=1:NBands, sp2=1:2)
               g4iw(o1,sp1,o2,sp2,w1,w2,wb) = g4iw(o1,sp1,o2,sp2,w1,w2,wb) + &
                     this_g2iw(o1,sp1,w1,w1-wb)*this_g2iw(o2,sp2,w2-wb,w2)
            end forall
            ! Crossing term in case both orbitals are equal
            forall(o1=1:NBands, sp1=1:2)
               g4iw(o1,sp1,o1,sp1,w1,w2,wb) = g4iw(o1,sp1,o1,sp1,w1,w2,wb) - &
                     this_g2iw(o1,sp1,w1,w2)*this_g2iw(o1,sp1,w2-wb,w1-wb)
            end forall
         end forall
      endif
   endif

   ! For the particle-particle channel, we only need to perform the frequency
   ! substitution wb -> wb-w1-w2 in the Green's functions (see Rohringer et al.
   ! 2012, arXiv:1202.2796)
   if(allocated(g4iw_pp)) then
      if(ZphConv.eq.0) then
         forall(w1=-N4iwf:N4iwf-1, w2=-N4iwf:N4iwf-1, wb=-N4iwb:N4iwb)
            forall(o1=1:NBands, sp1=1:2, o2=1:NBands, sp2=1:2)
               g4iw_pp(o1,sp1,o2,sp2,w1,w2,wb) = g4iw_pp(o1,sp1,o2,sp2,w1,w2,wb) +&
                     this_g2iw(o1,sp1,wb-w2-1,w1)*this_g2iw(o2,sp2,w2,wb-w1-1)
            end forall
            ! Crossing term in case both orbitals are equal
            forall(o1=1:NBands, sp1=1:2)
               g4iw_pp(o1,sp1,o1,sp1,w1,w2,wb) = g4iw_pp(o1,sp1,o1,sp1,w1,w2,wb) -&
                     this_g2iw(o1,sp1,w2,w1)*this_g2iw(o1,sp1,wb-w2-1,wb-w1-1)
            end forall
         end forall
      else
         stop 'G4iwPP Measurement not implemented for adapted convention'
      endif
   endif


   if (timings_maxorder >= DTrace%noper) then
      call wasted(end_secs, end_nanos)
      timings_g4iw_add(DTrace%noper) = timings_g4iw_add(DTrace%noper) + &
                    (end_secs - start_secs + (end_nanos - start_nanos)*1.0e-9)
      count_g4iw(DTrace%noper) = count_g4iw(DTrace%noper) + 1
!      write (0,'(I5,F10.6,2X,F10.6)') &
!          DTrace%noper, end_secs - start_secs + (end_nanos - start_nanos)*1.0e-9, &
!          timings_g4iw(DTrace%noper)/count_g4iw(DTrace%noper)
   endif
end subroutine


!> \brief Measures the imaginary time Green´s function in Legendre polynomials.
!!
!! Measures all possible imaginary-time one-particle Green´s functions
!! \f[
!!        G(\tau) = -\langle T_\tau c(\tau)c^+(0)\rangle
!! \f]
!! from the current trace. Two optimisations are made:
!!
!!  # Instead of expensive insertion of operators in the local trace, instead
!!    two hybridisation lines in the bath are removed (Gull 2008, §7.2)
!!
!!  # Instead of measuring in tau bins or Matsubara frequencies, scaled Legendre
!!    polynomials in \f$ x = 2\tau/\beta - 1 \f$ are used, which provides better
!!    statistics and better convergence at the anti-periodic step (Boehnke 2011)
!!
subroutine MeasGtauRem(isDetRat)
!   type(TTrace), intent(in)    :: DTrace
!   type(TStates), intent(in)   :: DStates

   !isDetRat==false measure using MInv, isDetRat==true measure using determinant ratio
   logical, intent(in)           :: isDetRat
   integer                       :: TFPos(NBands,2,2)  ! Position in trace (orbital,spin,ca)
   integer                       :: idt,iLegMax,l,iB,iS
   type(TOper),pointer           :: opa,opc
   type(TSubMatrix),pointer      :: tmpFullHybr(:,:), ReducedHybr(:,:)
   real(KINDR)                   :: M,x,Plp1x,Plx,Plm1x
   type(TLogDet)                 :: Det_G, Det_Z, tDet
!   type(TOperPointer)         :: Oper(2)
   type(TStates),pointer         :: DStates
   type(TTrace),pointer          :: DTrace
   type(TStates_pointer)         :: pDStates
   type(TTrace_pointer)          :: pDTrace

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   nullify(tmpFullHybr);nullify(ReducedHybr) 

   if(isDetRat) then
      !extract the old determinant from the full hybridisation matrix
      tmpFullHybr=>get_FullHybr(DTrace,DStates,DET_Z)
      !allocating the same space for the reduced matrix
      !by calling the same function out of lazyness
      ReducedHybr=>get_FullHybr(DTrace,DStates,Det_G)
   endif
   
   ! Iterate through all the annihilation operators of the trace. TFPos tracks
   ! indices in the block-diagonal Minv matrix, which is an (orbital,spin)-indexed
   ! array of allocatable hybridisation sub-matrices. The rows (columns) of
   ! these sub-matrices are ordered by the order of occurrence of the annihilation
   ! (creation) operators of the local trace.
   opa=>DTrace%first
   TFPos=0
   do while(associated(opa))
     if(opa%CA.eq.2) then      ! 2 = Annihilation operator
       TFPos(opa%Orbital,opa%Spin,2)=TFPos(opa%Orbital,opa%Spin,2)+1
       TFPos(opa%Orbital,opa%Spin,1)=0      ! reset creation part

       ! Iterate through all the creation operators of the trace with matching
       ! quantum numbers. Update column indexing of the corresponding (orbital,spin)
       ! - submatrix
       opc => DTrace%ffirst(opa%Orbital, opa%Spin)%p
       do while(associated(opc))
         if(opc%CA.eq.1)then
           TFPos(opc%Orbital,opc%Spin,1)=TFPos(opc%Orbital,opc%Spin,1)+1

           if(.not.isDetRat) then
               ! Read out hybridisation matrix elements (obital,spin)(column(c),row(a))
               ! Find corresponding tau bin and (-1,1)-scaled tau coordinate
               M=DTrace%MInv(opc%Orbital,opc%Spin)%Mat &
                        (TFPos(opc%Orbital,opc%Spin,1),TFPos(opc%Orbital,opc%Spin,2))
           else
               !PATRIK'S DETERMINANT HACK
               !set column and row of hybr matrix to 0 except for intersect to 1
               ReducedHybr(opc%Orbital,opc%Spin)%Mat&
                        (:,TFPos(opc%Orbital,opc%Spin,1))=0d0
               ReducedHybr(opc%Orbital,opc%Spin)%Mat&
                        (TFPos(opc%Orbital,opc%Spin,2),:)=0d0         
               ReducedHybr(opc%Orbital,opc%Spin)%Mat&
                        (TFPos(opc%Orbital,opc%Spin,2),TFPos(opc%Orbital,opc%Spin,1))=1d0
               
               Det_G = TLogDet(log=0.0_KINDR, sign=1.0_KINDR)
               
               do iB=1,NBands
                     do iS=1,2
                        if(associated(ReducedHybr(iB,iS)%Mat))then
                           if(size(ReducedHybr(iB,iS)%Mat(1,:)).gt.0)then
                              call get_LogDetFull(size(ReducedHybr(iB,iS)%Mat(1,:)),ReducedHybr(iB,iS)%Mat,tDet)
                              Det_G=tDet*Det_G
                           endif
                        endif
                     enddo
               enddo
               
               !calculate sign of row, column
               !detSign=(-1)**(mod(TFPos(opc%Orbital,opc%Spin,2)+TFPos(opc%Orbital,opc%Spin,1),2))
               !detSign=1
               
               !this is the determiant ratio
               M=detval(Det_G/Det_Z)
               
               
               !reset colum and row of reduced hybr with actual values
               ReducedHybr(opc%Orbital,opc%Spin)%Mat&
                        (:,TFPos(opc%Orbital,opc%Spin,1))=&
                        tmpFullHybr(opc%Orbital,opc%Spin)%Mat(:,TFPos(opc%Orbital,opc%Spin,1))
               ReducedHybr(opc%Orbital,opc%Spin)%Mat&
                        (TFPos(opc%Orbital,opc%Spin,2),:)=&
                        tmpFullHybr(opc%Orbital,opc%Spin)%Mat(TFPos(opc%Orbital,opc%Spin,2),:)         
               ReducedHybr(opc%Orbital,opc%Spin)%Mat&
                        (TFPos(opc%Orbital,opc%Spin,2),TFPos(opc%Orbital,opc%Spin,1))=&
                        tmpFullHybr(opc%Orbital,opc%Spin)%Mat&
                        (TFPos(opc%Orbital,opc%Spin,2),TFPos(opc%Orbital,opc%Spin,1))
           endif
           
           if(opa%tau-opc%tau.gt.0d0)then
               idt=int((opa%tau-opc%tau)/DTrace%beta*dble(NGtau))+1
               x=2d0/DTrace%Beta*(opa%tau-opc%tau)-1d0
           else ! take into account beta periodicity
               idt=int((DTrace%beta+opa%tau-opc%tau)/DTrace%beta*dble(NGtau))+1
               x=2d0/DTrace%Beta*(DTrace%Beta+opa%tau-opc%tau)-1d0
               M=-M
           endif

           ! Probably better:
           !   beta = DTrace%beta
           !   dtau = opa%tau-opc%tau
           !   if(dtau<0) dtau = dtau + beta
           !   idt = int(dtau/beta*NGtau) + 1
           !   x = 2/DTrace%Beta * dtau - 1

           ! taking care of signum
           M=M*DTrace%cfgsign

           ! Measure to Legendre polynomials
           Plx=0d0
           Plm1x=0d0
           do iLegMax=0,NLegMax-1
              if(x.eq.-1d0.or.x.eq.1d0)then
                 Plp1x=(x)**iLegMax
              endif
              if(iLegMax.eq.0)then
                 Plp1x=1d0
              elseif(iLegMax.eq.1)then
                 Plp1x=x
              else
                 l=iLegMax-1
                 Plp1x=(2d0*dble(l)+1d0)/dble(l+1)*x*Plx-dble(l)/dble(l+1)*Plm1x
              endif
              GLeg(opa%Orbital,opa%Spin,iLegMax+1)=&
                 GLeg(opa%Orbital,opa%Spin,iLegMax+1)&
                 +M*Plp1x
              Plm1x=Plx
              Plx=Plp1x
           enddo

           ! Measure to tau bins
           if(idt.gt.NGtau) idt=NGtau
           Gtau(opa%Orbital,opa%Spin,idt)=&
              Gtau(opa%Orbital,opa%Spin,idt)+M
           if (allocated(Gtau_mean_step)) Gtau_mean_step(opa%Orbital,opa%Spin,g_inmeas)=&
              Gtau_mean_step(opa%Orbital,opa%Spin,g_inmeas) + M * DTrace%cfgsign
           if (allocated(Gtau_mid_step)) then
              if (abs(opa%tau - opc%tau)/DTrace%beta >= 0.4_KINDR&
                  .and. abs(opa%tau - opc%tau)/DTrace%beta <= 0.6_KINDR)&
                 Gtau_mid_step(opa%Orbital,opa%Spin,g_inmeas)=&
                    Gtau_mid_step(opa%Orbital,opa%Spin,g_inmeas) + M * DTrace%cfgsign
           end if
           if (allocated(sign_step)) sign_step(g_inmeas) = DTrace%cfgsign
#ifdef USE_NFFT
           if (b_Giw_lookup) then
              M = - DTrace%cfgsign * DTrace%MInv(opc%Orbital,opc%Spin)%Mat&
                 (TFPos(opc%Orbital,opc%Spin,1),TFPos(opc%Orbital,opc%Spin,2))
              Giw_lookup_tau(Giw_lookup_row(opa%Orbital,opa%Spin),opa%Orbital,opa%Spin)&
                 = (opa%tau - opc%tau)/(2.0 * DTrace%beta)
              Giw_lookup_M(Giw_lookup_row(opa%Orbital,opa%Spin),opa%Orbital,opa%Spin)&
                 = M * exp(cmplx(0.,(opa%tau - opc%tau) * PI / DTrace%beta, kind=KINDC))

              Giw_lookup_row(opa%Orbital,opa%Spin) = Giw_lookup_row(opa%Orbital,opa%Spin) + 1

              if (Giw_lookup_row(opa%Orbital,opa%Spin) == NLookup_nfft + 1) then
                 call MeasGiwFlushLookup(opa%Orbital, opa%Spin)
              end if
           end if
#endif
         endif
         opc => opc%fnext
       enddo
     endif
     opa=>opa%next
   enddo
   
   if(isDetRat) then
      !remove the temp hybridisation matrix again
      do iB=1,NBands
         do iS=1,2
            deallocate(tmpFullHybr(iB,iS)%Mat)
            deallocate(ReducedHybr(iB,iS)%Mat)
         enddo
      enddo
   
      deallocate(tmpFullHybr)
      deallocate(ReducedHybr)
   endif
   
end subroutine MeasGtauRem

subroutine MeasGiwFlushLookup(Orbital, Spin)
   integer, intent(in)  :: Orbital, Spin
#ifdef USE_NFFT
   complex(KINDC)       :: giw_bunch(nfreq, 2)

   giw_bunch = get_giw2(Giw_lookup_row(Orbital,Spin) - 1,&
                        Giw_lookup_M(1:Giw_lookup_row(Orbital,Spin) - 1,Orbital,Spin),&
                        Giw_lookup_tau(1:Giw_lookup_row(Orbital,Spin) - 1,Orbital,Spin))
   Giw(Orbital, Spin, :) = Giw(Orbital, Spin, :) + giw_bunch(:, 1)
   Giw_lookup_row(Orbital,Spin) = 1
#else
   stop "Unreachable: MeasGiwFlushLookup called in code compiled without USE_NFFT"
#endif
end subroutine MeasGiwFlushLookup

!> \brief Measures the imaginary time Green´s function in Legendre polynomials.
!!
!! Measures all possible imaginary-time one-particle Green´s functions
!! \f[
!!        G(\tau) = -\langle T_\tau c(\tau)c^+(0)\rangle
!! \f]
!! from the current trace. Two optimisations are made:
!!
!!  # Instead of expensive insertion of operators in the local trace, instead
!!    two hybridisation lines in the bath are removed (Gull 2008, §7.2)
!!
!!  # Instead of measuring in tau bins or Matsubara frequencies, scaled Legendre
!!    polynomials in \f$ x = 2\tau/\beta - 1 \f$ are used, which provides better
!!    statistics and better convergence at the anti-periodic step (Boehnke 2011)
!!
!===============================================================================
!> Regenerate all arrays containing the possible hybridizations after a global
!! update and returns the hybrizization matrix.
subroutine MeasGtauRem_full()
!===============================================================================
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
!local   
   integer                    :: iters, itere
   integer                    :: pos(2),b1,s1,b2,s2
   type(TOper),pointer        :: ElementS,ElementE
   real(KINDR)                :: m, x
   real(KINDR)                :: Plp1x, Plx, Plm1x
   integer                    :: iLegMax, l, idt

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   if(dtrace%noper.eq.0)then
      return
   endif


   pos=1
   do b1=1,dstates%nbands
   do s1=1,2
   !write(*,*) "b1,s1", b1,s1

   ElementS=>dtrace%first
   do iters=1,dtrace%NOper

      if(elements%orbital.eq.b1.and.elements%spin.eq.s1)then !select1

      if(ElementS%has_hyb) then
         if(ElementS%CA.eq.1)then

            pos(2)=1
            do b2=1,dstates%nbands
            do s2=1,2

            ElementE=>dtrace%first
            do itere=1,dtrace%NOper

               if(elemente%orbital.eq.b2.and.elemente%spin.eq.s2)then !select2

               if(ElementE%has_hyb) then
                  if(ElementE%CA.eq.2)then

                     if(size(DTrace%Minv_full%mat).ne.0)then

                        M=DTrace%Minv_full%mat(pos(1),pos(2))

                        if(elemente%tau-elements%tau.gt.0d0)then
                            idt=int((elemente%tau-elements%tau)/DTrace%beta*dble(NGtau))+1
                            x=2d0/DTrace%Beta*(elemente%tau-elements%tau)-1d0
                        else ! take into account beta periodicity
                            idt=int((DTrace%beta+elemente%tau-elements%tau)/DTrace%beta*dble(NGtau))+1
                            x=2d0/DTrace%Beta*(DTrace%Beta+elemente%tau-elements%tau)-1d0
                            M=-M
                        endif

                        M=M*DTrace%cfgsign

                        Gtau_full(elemente%Orbital,elemente%Spin,elements%Orbital,elements%Spin,idt)=&
                           Gtau_full(elemente%Orbital,elemente%Spin,elements%Orbital,elements%Spin,idt)+M

                        ! Measure to Legendre polynomials
                        Plx=0d0
                        Plm1x=0d0
                        do iLegMax=0,NLegMax-1
                           if(x.eq.-1d0.or.x.eq.1d0)then
                              Plp1x=(x)**iLegMax
                           endif
                           if(iLegMax.eq.0)then
                              Plp1x=1d0
                           elseif(iLegMax.eq.1)then
                              Plp1x=x
                           else
                              l=iLegMax-1
                              Plp1x=(2d0*dble(l)+1d0)/dble(l+1)*x*Plx-dble(l)/dble(l+1)*Plm1x
                           endif
                           GLeg_full(elemente%Orbital,elemente%Spin,elements%Orbital,elements%Spin,iLegMax+1)=&
                              GLeg_full(elemente%Orbital,elemente%Spin,elements%Orbital,elements%Spin,iLegMax+1)&
                              +M*Plp1x
                           Plm1x=Plx
                           Plx=Plp1x
                        enddo

                     endif


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

end subroutine MeasGtauRem_full

!> \brief Measures the two-particle imaginary time Green´s function in
!!        Legendre polynomials and one bosonic Matsubara frequency.
!!
!! Measures all possible imaginary-time two-particle Green´s functions
!! \f[
!!      G(\tau_{12}, \tau_{34}, \tau_{14}) =
!!          +\langle T_tau c(\tau_1) c^+(\tau_2) d(\tau_3) d^+(\tau_4) \rangle
!! \f]
!! from the current trace, where \f$ c \f$ and \f$ d \f$ denote fermionic
!! operators to possibly distinct quantum numbers (sometimes written as
!! \f$ G_{AABB} \f$).
!!
!! The chosen time differences \f$ \tau_ij := \tau_i - \tau_j \f$ yield two
!! beta-anti-periodic (1-2, 3-4) and one beta-periodic (1-4) time coordinate,or,
!! equivalently, two fermionic and one bosonic Matsubara frequencies. Following
!! (Boehnke 2011), the bosonic frequency is retained, while the other time
!! coordinates are mapped to scaled Legendre polynomials.
!!
!! Instead of expensive insertion of operators in the local trace, instead two
!! hybridisation lines in the bath are removed (Gull 2008, §7.2). The term
!! "remove" may be misleading here, since we can multiply invoke that procedure
!! if two operators in the two-particle Green´s function coincide.
!!
subroutine MeasGtau4PntLeg()
    use LegendrePoly, only: legendrep, PI

    integer             :: TF1Pos(NBands,2,2),TF2Pos(NBands,2,2)
    type(TOper),pointer :: opa1, opa2, opc1, opc2     ! trace pointers
    real(KINDR)         :: M
!    type(TOperPointer)  :: Operators(2)

    ! tau
    real(KINDR)         :: tau12, tau34, tau14, x12, x34, beta, sgn12, sgnfull
    integer             :: bin12, bin34, bin14

    ! Caching
    real(KINDR)         :: leg12(0:N4leg-1)
    real(KINDR)         :: leg34(0:N4leg-1)
    complex(KINDC)      :: ph14(-N4iwb:N4iwb)
    complex(KINDC)      :: leg12ph14(0:N4tau-1, -N4iwb:N4iwb)

    integer             :: l, lp, wm
    integer             :: nins
    type(TTrace),pointer       :: DTrace
    type(TStates_pointer)      :: pDStates
    type(TTrace_pointer)       :: pDTrace

 
   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr


    TF1Pos=0
    beta = DTrace%beta
    nins = 0
    M = 0

    opa1=>null(); opc1=>null(); opc2=>null(); opa2=>null()

    opa1=>DTrace%first          ! iterates over c(tau_1)
    do while(associated(opa1))
      if(opa1%CA.eq.2)then
        TF1Pos(opa1%Orbital,opa1%Spin,2)=TF1Pos(opa1%Orbital,opa1%Spin,2)+1
        TF1Pos(opa1%Orbital,opa1%Spin,1)=0

        opc1 => DTrace%ffirst(opa1%Orbital, opa1%Spin)%p      ! iterates over c^+(tau_2)
        do while(associated(opc1))
          if(opc1%CA.eq.1)then
            TF1Pos(opc1%Orbital,opc1%Spin,1)=TF1Pos(opc1%Orbital,opc1%Spin,1)+1
            TF2Pos=0

            ! determine stuff for the first pair 1 <--> 2
            ! Find out tau coordinate of the differences, use beta periodicity.
            ! Track the sign for the fermionic part for the Legendre coefficients
            ! (Boehnke 2011, C20) and the imaginary time measurement (C16)
            tau12 = opa1%tau - opc1%tau
            sgn12 = sign(1._KINDR, tau12)
            if(tau12 < 0) tau12 = tau12 + beta
            if(IAND(FourPnt, GF4_LEGENDRE) /= 0) then
                x12 = 2/DTrace%Beta * tau12 - 1
                call legendrep(x12, leg12)
            endif

            opc2=>DTrace%first  ! iterates over c^+(tau_4)
            do while(associated(opc2))
              if(opc2%CA.eq.1) then
                TF2Pos(opc2%Orbital,opc2%Spin,1)=TF2Pos(opc2%Orbital,opc2%Spin,1)+1
                TF2Pos(opc2%Orbital,opc2%Spin,2) = 0

                ! determine stuff for the bosonic pair 1 <--> 4
                tau14 = opa1%tau - opc2%tau
                ! no sign here!
                if(tau14 < 0) tau14 = tau14 + beta
                if( IAND(FourPnt, GF4_LEGENDRE) /= 0 ) then
                    forall(wm = -N4iwb:N4iwb) &
                        ph14(wm) = exp(cmplx(0.D0, 2.D0*wm*PI*tau14/beta, KINDR))
                    ! Caching this quantity for faster access in the multiplication,
                    ! which is very heavy
                    forall(l = 0:N4leg-1) &
                        leg12ph14(l,:) = leg12(l)*ph14
                endif

                opa2 => DTrace%ffirst(opc2%Orbital, opc2%Spin)%p  ! iterates over c(tau_3)
                do while(associated(opa2))
                  if(opa2%CA.eq.2)then
                      TF2Pos(opa2%Orbital,opa2%Spin,2)=TF2Pos(opa2%Orbital,opa2%Spin,2)+1

                      M=0d0
                      if(opc1%Orbital.eq.opc2%Orbital.and.opc1%Spin.eq.opc2%Spin)then
                         M=DTrace%MInv(opc1%Orbital,opc1%Spin)%Mat&
                            (TF1Pos(opc1%Orbital,opc1%Spin,1),TF1Pos(opc1%Orbital,opc1%Spin,2))*&
                           DTrace%MInv(opc2%Orbital,opc2%Spin)%Mat&
                            (TF2Pos(opc2%Orbital,opc2%Spin,1),TF2Pos(opc2%Orbital,opc2%Spin,2))-&
                           DTrace%MInv(opc1%Orbital,opc1%Spin)%Mat&
                            (TF2Pos(opc1%Orbital,opc1%Spin,1),TF1Pos(opc1%Orbital,opc1%Spin,2))*&
                           DTrace%MInv(opc2%Orbital,opc2%Spin)%Mat&
                            (TF1Pos(opc2%Orbital,opc2%Spin,1),TF2Pos(opc2%Orbital,opc2%Spin,2))
                      else !if(opc1%Orbital.ne.opc2%Orbital.or.opc1%Spin.ne.opc2%Spin)then
                         M=DTrace%MInv(opc1%Orbital,opc1%Spin)%Mat&
                            (TF1Pos(opc1%Orbital,opc1%Spin,1),TF1Pos(opc1%Orbital,opc1%Spin,2))*&
                           DTrace%MInv(opc2%Orbital,opc2%Spin)%Mat&
                            (TF2Pos(opc2%Orbital,opc2%Spin,1),TF2Pos(opc2%Orbital,opc2%Spin,2))
                      endif

                      ! determine stuff for second fermionic pair 3 <--> 4
                      tau34 = opa2%tau - opc2%tau
                      sgnfull = sgn12*sign(1._KINDR, tau34)
                      if(tau34 < 0) tau34 = tau34 + beta
                      if( IAND(FourPnt, GF4_LEGENDRE) /= 0 ) then
                          x34 = 2/DTrace%Beta * tau34 - 1
                          call legendrep(x34, leg34)
                      endif

                      ! add sign and include M into leg34
                      nins = nins + 1
                      M = -sgnfull*M

                      ! perform the mapping onto mixed Legendre/Matsubara basis
                      if( IAND(FourPnt, GF4_LEGENDRE) /= 0 ) then
                        ! Include M in leg34
                        leg34 = leg34*M
                        forall(lp = 0:N4leg-1)
                          ! includes M
                          G4leg(opc1%Orbital,opc1%Spin,opc2%Orbital,&
                                        opc2%Spin,:,lp,:) = &
                              G4leg(opc1%Orbital,opc1%Spin,opc2%Orbital,&
                                            opc2%Spin,:,lp,:) + leg34(lp)*leg12ph14
                        end forall
                      endif

                      ! Perform the direct measurement to the imaginary time axis
                      if( IAND(FourPnt, GF4_IMAGTIME) /= 0 ) then
                          ! Compute bins
                          bin12 = int(N4tau * tau12/beta) + 1
                          bin14 = int(N4tau * tau14/beta) + 1
                          bin34 = int(N4tau * tau34/beta) + 1

                          ! We removed the averaging guard, because the index are
                          ! purely fermionic or bosonic, therefore we should not
                          ! "blur" any indices.
                          G4tau(opa1%Orbital,opa1%Spin,opa2%Orbital,&
                                        opa2%Spin,bin12,bin34,bin14) = &
                              G4tau(opa1%Orbital,opa1%Spin,opa2%Orbital,&
                                            opa2%Spin,bin12,bin34,bin14) + M
                      endif

                  endif
                  opa2=>opa2%fnext
                enddo
              endif

              opc2=>opc2%next
            enddo
          endif

          opc1=>opc1%fnext
        enddo
      endif

      opa1=>opa1%next
    enddo
    !write (*,*) 'Trace size:', DTrace%NOper, ' Ins:', nins
end subroutine MeasGtau4PntLeg


subroutine MeasSingleOcc_from_Densitymatrix()

    type(tstates), pointer    :: dstates
    type(tstates_pointer)     :: pdstates
    real(kindr), allocatable  :: tempmat(:, :)
    integer                   :: iB1, iS1, iSSt, iSt, nsst, offset, NStates_sst, NStates_nsst

    pdstates = transfer(ipdstates, pdstates)
    dstates => pdstates%ptr

    do iB1=1,NBands
       do iS1=1,2
          do iSSt=0,DStates%NSStates-1
             nsst = DStates%Substates(iSSt)%Connect(ib1, is1, 2)
             if (nsst == -1) cycle

             offset = DStates%SubStates(iSSt)%offset + 1
             NStates_sst = DStates%Substates(iSSt)%NStates
             NStates_nsst = DStates%Substates(nsst)%NStates

             allocate(tempmat(NStates_sst, NStates_nsst))

             ! occ increment per state is <iSt|rho*psi^dag(ib1, is1)*psi(ib1, is1)|iSt>
             ! - <iSt|rho*psi^dag(ib1, is1) is the iSt-th row of tempmat
             ! - psi(ib1, is1)|iSt> is the iSt-th column of psi(ib1, is1)
             ! - we take their dot products as the end results in the loop
             ! analogous procedure for doubleocc, rho1, rho2 in subroutines below

             call DGEMM('N', 'N', NStates_sst, NStates_nsst, NStates_sst,&
                        1.0d0, densitymatrix(offset, offset), NStates,&
                        DStates%SubStates(nsst)%Psis_EB(ib1, is1, 1)%Op, NStates_sst,&
                        0.0d0, tempmat, NStates_sst)

             do iSt = 1, NStates_sst
                Single_occ(ib1,is1) = Single_occ(ib1,is1)&
                   + DDOT(NStates_nsst, tempmat(iSt, 1), NStates_sst,&
                          DStates%SubStates(iSSt)%Psis_EB(ib1, is1, 2)%Op(1, iSt), 1)
             enddo
             deallocate(tempmat)
          enddo
       enddo
    enddo

end subroutine MeasSingleOcc_from_Densitymatrix


subroutine MeasDoubleOcc_from_Densitymatrix()

    type(tstates),pointer         :: dstates
    type(tstates_pointer)         :: pdstates
    real(kindr), allocatable      :: tempmat1(:, :), tempmat2(:, :)
    integer                       :: iB1, iS1, iB2, iS2, iSSt, iSt, nsst1, nsst2, offset
    integer                       :: NStates_sst, NStates_nsst1, NStates_nsst2

    pdstates = transfer(ipdstates, pdstates)
    dstates => pdstates%ptr

    do iB1=1,NBands
    do iS1=1,2
    do iB2=1,NBands
    do iS2=1,2
    do iSSt=0,DStates%NSStates-1
       nsst1 = DStates%Substates(iSSt)%Connect(ib1, is1, 2)
       if (nsst1 == -1) cycle
       nsst2 = DStates%Substates(iSSt)%Connect(ib2, is2, 2)
       if (nsst2 == -1) cycle

       offset = DStates%SubStates(iSSt)%offset + 1
       NStates_sst = DStates%Substates(iSSt)%NStates
       NStates_nsst1 = DStates%Substates(nsst1)%NStates
       NStates_nsst2 = DStates%Substates(nsst2)%NStates

       allocate(tempmat1(NStates_sst, NStates_nsst2))

       call DGEMM('N', 'N', NStates_sst, NStates_nsst2, NStates_sst,&
                  1.0d0, densitymatrix(offset, offset), NStates,&
                  DStates%SubStates(nsst2)%Psis_EB(ib2, is2, 1)%Op, NStates_sst,&
                  0.0d0, tempmat1, NStates_sst)

       allocate(tempmat2(NStates_sst, NStates_sst))

       call DGEMM('N', 'N', NStates_sst, NStates_sst, NStates_nsst2,&
                  1.0d0, tempmat1, NStates_sst,&
                  DStates%SubStates(iSSt)%Psis_EB(ib2, is2, 2)%Op, NStates_nsst2,&
                  0.0d0, tempmat2, NStates_sst)

       deallocate(tempmat1)
       allocate(tempmat1(NStates_sst, NStates_nsst1))

       call DGEMM('N', 'N', NStates_sst, NStates_nsst1, NStates_sst,&
                  1.0d0, tempmat2, NStates_sst,&
                  DStates%SubStates(nsst1)%Psis_EB(ib1, is1, 1)%Op, NStates_sst,&
                  0.0d0, tempmat1, NStates_sst)

       deallocate(tempmat2)

       do iSt = 1, NStates_sst
          occ(ib1, is1, ib2, is2) = occ(ib1, is1, ib2, is2)&
             + DDOT(NStates_nsst1, tempmat1(iSt, 1), NStates_sst,&
                    DStates%SubStates(iSSt)%Psis_EB(ib1, is1, 2)%Op(1, iSt), 1)
       enddo
       deallocate(tempmat1)
    enddo
    enddo
    enddo
    enddo 
    enddo

end subroutine MeasDoubleOcc_from_Densitymatrix


subroutine Meas_rho2_from_Densitymatrix()

    type(tstates), pointer   :: dstates
    type(tstates_pointer)    :: pdstates
    real(kindr), allocatable :: tempmat1(:, :), tempmat2(:, :)
    integer                  :: iB1, iS1, iB2, iS2, iB3, iS3, iB4, iS4
    integer                  :: iSSt, iSt, nsst1, nsst2, nsst3, offset
    integer                  :: NStates_sst, NStates_nsst1, NStates_nsst2, NStates_nsst3

    pdstates = transfer(ipdstates, pdstates)
    dstates => pdstates%ptr

    do iB1=1,NBands
    do iS1=1,2
    do iB2=1,NBands
    do iS2=1,2   
    do iB3=1,NBands
    do iS3=1,2
    do iB4=1,NBands
    do iS4=1,2   
    do iSSt=0,DStates%NSStates-1
       nsst1 = DStates%Substates(iSSt)%Connect(ib1, is1, 2)
       if (nsst1 == -1) cycle
       nsst2 = DStates%Substates(nsst1)%Connect(ib2, is2, 2)
       if (nsst2 == -1) cycle
       nsst3 = DStates%Substates(nsst2)%Connect(ib3, is3, 1)
       if (nsst3 == -1) cycle
       if (DStates%Substates(nsst3)%Connect(ib4, is4, 1) /= iSSt) cycle

       offset = DStates%SubStates(iSSt)%offset + 1
       NStates_sst = DStates%Substates(iSSt)%NStates
       NStates_nsst1 = DStates%Substates(nsst1)%NStates
       NStates_nsst2 = DStates%Substates(nsst2)%NStates
       NStates_nsst3 = DStates%Substates(nsst3)%NStates

       allocate(tempmat1(NStates_sst, NStates_nsst3))

       call DGEMM('N', 'N', NStates_sst, NStates_nsst3, NStates_sst,&
                  1.0d0, densitymatrix(offset, offset), NStates,&
                  DStates%SubStates(nsst3)%Psis_EB(ib4, is4, 1)%Op, NStates_sst,&
                  0.0d0, tempmat1, NStates_sst)

       allocate(tempmat2(NStates_sst, NStates_nsst2))

       call DGEMM('N', 'N', NStates_sst, NStates_nsst2, NStates_nsst3,&
                  1.0d0, tempmat1, NStates_sst,&
                  DStates%SubStates(nsst2)%Psis_EB(ib3, is3, 1)%Op, NStates_nsst3,&
                  0.0d0, tempmat2, NStates_sst)

       deallocate(tempmat1)
       allocate(tempmat1(NStates_sst, NStates_nsst1))

       call DGEMM('N', 'N', NStates_sst, NStates_nsst1, NStates_nsst2,&
                  1.0d0, tempmat2, NStates_sst,&
                  DStates%SubStates(nsst1)%Psis_EB(ib2, is2, 2)%Op, NStates_nsst2,&
                  0.0d0, tempmat1, NStates_sst)

       deallocate(tempmat2)

       do iSt = 1, NStates_sst
          rho2(2*ib4-2+is4, 2*ib3-2+is3, 2*ib2-2+is2, 2*ib1-2+is1)=&
             rho2(2*ib4-2+is4, 2*ib3-2+is3, 2*ib2-2+is2, 2*ib1-2+is1)&
             + DDOT(NStates_nsst1, tempmat1(iSt, 1), NStates_sst,&
                    DStates%SubStates(iSSt)%Psis_EB(ib1, is1, 2)%Op(1, iSt), 1)
       enddo
       deallocate(tempmat1)
    enddo
    enddo
    enddo 
    enddo
    enddo
    enddo
    enddo 
    enddo
    enddo

end subroutine Meas_rho2_from_Densitymatrix


subroutine Meas_rho1_from_Densitymatrix()

    type(tstates), pointer        :: dstates
    type(tstates_pointer)         :: pdstates
    real(kindr), allocatable      :: tempmat(:, :)
    integer                       :: iB1, iS1, iB2, iS2, iSSt, iSt, nsst, offset
    integer                       :: NStates_sst, NStates_nsst

    pdstates = transfer(ipdstates, pdstates)
    dstates => pdstates%ptr

    do iB1=1,NBands
    do iS1=1,2
    do iB2=1,NBands
    do iS2=1,2
    do iSSt=0,DStates%NSStates-1
       nsst = DStates%Substates(iSSt)%Connect(ib1, is1, 2)
       if (nsst == -1) cycle
       if (DStates%Substates(nsst)%Connect(ib2, is2, 1) /= iSSt) cycle

       offset = DStates%SubStates(iSSt)%offset + 1
       NStates_sst = DStates%Substates(iSSt)%NStates
       NStates_nsst = DStates%Substates(nsst)%NStates

       allocate(tempmat(NStates_sst, NStates_nsst))


       call DGEMM('N', 'N', NStates_sst, NStates_nsst, NStates_sst,&
                  1.0d0, densitymatrix(offset, offset), NStates,&
                  DStates%SubStates(nsst)%Psis_EB(ib2, is2, 1)%Op, NStates_sst,&
                  0.0d0, tempmat, NStates_sst)

       do iSt = 1, NStates_sst
          rho1(2*ib2-2+is2, 2*ib1-2+is1) = rho1(2*ib2-2+is2, 2*ib1-2+is1)&
             + DDOT(NStates_nsst, tempmat(iSt, 1), NStates_sst,&
                    DStates%SubStates(iSSt)%Psis_EB(ib1, is1, 2)%Op(1, iSt), 1)
       enddo
       deallocate(tempmat)
    enddo
    enddo
    enddo
    enddo
    enddo

end subroutine Meas_rho1_from_Densitymatrix


! ! > \brief Measures the density matrix \f$ | \psi \rangle  \langle \psi | \f$ 
! ! at \f$ \beta/2 \f$.
! ! 
! ! Measures the density matrix
! ! \f[
! !     | \phi \rangle \langle \phi |  = \frac { \sum_m | \psi_m \rangle 
! !                                                     \langle \psi_m | 
! !                                        }
! !                                     { \sum_m \langle \psi_m | 
! !                                       \psi_m \rangle }
! !                             = \frac {\sum_m w_m 
! !                                             \frac {|psi_m \rangle
! !                                                    \langle \psi_m |}
! !                                                   { \langle \psi_m |
! !                                                    \psi_m \rangle} }
! !                                     { \sum_m w_m}
! ! \f]
! ! where \f$ m \f$ runs over all states considered in the local trace and  
! ! \f$ w_m \f$ is the value of one summand of the local trace. We need to 
! ! renormalise the wave vectors since they are not normalised before they
! ! are stored. The measurement is performed close to \f$ \beta/2 \f$ so that 
! ! the effect of the truncation at \f$ 0 \f$ and \f$ \beta \f$ is small. 

! Measures the Densitymatrix with the states considered in truncation, that are
! evolved to beta half
subroutine MeasDensityMatrix_beta_half()

   integer                    :: st1, st2, trst, offset, vec_outer_ind
   real(KINDR)                :: factor

   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace

   pDStates = transfer(ipDStates,pDStates)
   pDTrace = transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   call update_b2_states_eb(DTrace, DStates)

   ! if there are operators in the trace, use the states evolved to beta/2
   if(DTrace%NOper .ne. 0)then
      offset = dstates%substates(DTrace%sst_ket_b2)%offset
      do trst = DTrace%outer_state, DTrace%outer_state + DTrace%outer_sst_size - 1
         if (.not. associated(DTrace%ket_b2)&
             .or. .not. associated(dtrace%bra_b2)&
             .or. DTrace%sst_ket_b2 /= DTrace%sst_bra_b2)&
            stop "unreachable: beta/2-data invalid in MeasDensityMatrix_beta_half"

         if (b_statesampling) then
            vec_outer_ind = 1
         else
            vec_outer_ind = trst
         end if

         ! a needed normalization factor is
         !     norm = DTrace%parttrace(trst)/&
         !            DTrace%trace/&
         !            braket,
         ! where
         !     braket = dot_product(ket_b2(:, vec_outer_ind),&
         !                          bra_b2(:, vec_outer_ind)),
         ! but
         !     braket == parttrace(trst)
         ! up to numerical differences, so only
         !     norm = 1 / DTrace%Trace
         ! remains

         ! after this rearranging, braket and parttrace cancel, but it
         ! is probably still not a good idea to add contributions from
         ! configuration parts with zero weight contribution
         if (trval(dtrace%parttrace(trst)) == 0.0_KINDR) cycle

         do st1 = 1, size(dtrace%bra_b2(:, vec_outer_ind))
            factor = dtrace%bra_b2(st1, vec_outer_ind)&
                     / trval(DTrace%Trace /&
                             TLogTr(log = DTrace%ket_b2_logmax(vec_outer_ind)&
                                          + DTrace%bra_b2_logmax(vec_outer_ind),&
                                    sign = -DTrace%cfgsign))
            do st2 = 1, size(dtrace%ket_b2(:, vec_outer_ind))
               densitymatrix(offset+st1, offset+st2) =&
                  densitymatrix(offset+st1, offset+st2) +&
                  factor * dtrace%ket_b2(st2, vec_outer_ind)
            enddo
         enddo

      enddo

! otherwise the trace is empty and we use the eigenvectors of the local 
! hamiltonian
   else
      offset = dstates%substates(DTrace%outer_sst)%offset
      do trst = DTrace%outer_state, DTrace%outer_state + DTrace%outer_sst_size - 1
         densitymatrix(offset+trst, offset+trst) =&
            densitymatrix(offset+trst, offset+trst)&
            - DTrace%cfgsign * trval(DTrace%parttrace(trst)/DTrace%Trace)
      enddo 
   endif

end subroutine MeasDensityMatrix_beta_half


subroutine MeasSusz()
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace

   pDStates = transfer(ipDStates,pDStates)
   pDTrace = transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr
   call meas_ntau_n0(DTrace, DStates, ntau_n0)
end subroutine MeasSusz


!===============================================================================
!! This subroutine transforms the densitymatrix measured in eigenbasis to 
!! occupation basis.
subroutine Transform_DensityMatrix_EB_to_OCCB()
!===============================================================================
   integer                    :: isst, offset, ns
   type(TStates),pointer      :: DStates
   type(TStates_pointer)      :: pDStates

   pDStates = transfer(ipDStates,pDStates)
   DStates => pDStates%ptr

   do isst=0, DStates%NSStates-1

      offset = DStates%SubStates(isst)%offset
      ns = DStates%Substates(isst)%NStates

      densitymatrix(offset+1:offset+ns,offset+1:offset+ns)=&
         matmul(DStates%Substates(isst)%Evec,&
         matmul(densitymatrix(offset+1:offset+ns,offset+1:offset+ns),&
                transpose(DStates%Substates(isst)%Evec)))
   enddo

end subroutine Transform_DensityMatrix_EB_to_OCCB


!> Measures the expansion resolved density matrix.
!! We are going to use this to identify the local state the system is
!! in with certain expansion orders.
subroutine MeasExpResDensityMatrix()
   integer                       :: st1, st2, trst, offset, iB, iS, vec_outer_ind
   real(KINDR)                   :: factor

   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace

   pDStates = transfer(ipDStates,pDStates)
   pDTrace = transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   call update_b2_states_eb(DTrace, DStates)

   ! if there are operators in the trace, use the states evolved to beta/2
   if(DTrace%NOper .ne. 0)then
      offset = dstates%substates(DTrace%sst_ket_b2)%offset
      do trst = DTrace%outer_state, DTrace%outer_state + DTrace%outer_sst_size - 1
         if (.not. associated(DTrace%ket_b2)&
             .or. .not. associated(dtrace%bra_b2)&
             .or. DTrace%sst_ket_b2 /= DTrace%sst_bra_b2)&
            stop "unreachable: beta/2-data invalid in MeasExpResDensityMatrix"

         if (b_statesampling) then
            vec_outer_ind = 1
         else
            vec_outer_ind = trst
         end if

         ! a needed normalization factor is
         !     norm = DTrace%parttrace(trst)/&
         !            DTrace%trace/&
         !            braket,
         ! where
         !     braket = dot_product(ket_b2(:, vec_outer_ind),&
         !                          bra_b2(:, vec_outer_ind)),
         ! but
         !     braket == parttrace(trst)
         ! up to numerical differences, so only
         !     norm = 1 / DTrace%Trace
         ! remains

         ! after this rearranging, braket and parttrace cancel, but it
         ! is probably still not a good idea to add contributions from
         ! configuration parts with zero weight contribution
         if (trval(dtrace%parttrace(trst)) == 0.0_KINDR) cycle

         do st1 = 1, size(dtrace%bra_b2(:, vec_outer_ind))
            factor = dtrace%bra_b2(st1, vec_outer_ind)&
                     / trval(DTrace%Trace /&
                             TLogTr(log = DTrace%ket_b2_logmax(vec_outer_ind)&
                                          + DTrace%bra_b2_logmax(vec_outer_ind),&
                                    sign = -DTrace%cfgsign))
            do st2 = 1, size(dtrace%ket_b2(:, vec_outer_ind))
               do iB = 1, NBands
                  do iS = 1, 2
                     ExpResDensityMatrix(iB, iS, DTrace%NOSOper(iB, iS)/2, offset+st1, offset+st2) =&
                        ExpResDensityMatrix(iB, iS, DTrace%NOSOper(iB, iS)/2, offset+st1, offset+st2)&
                        + factor * dtrace%ket_b2(st2, vec_outer_ind)
                  end do
               end do
            enddo
         enddo

      enddo

! otherwise the trace is empty and we use the eigenvectors of the local 
! hamiltonian
   else
      offset = dstates%substates(DTrace%outer_sst)%offset
      do trst = DTrace%outer_state, DTrace%outer_state + DTrace%outer_sst_size - 1
         do iB = 1, NBands
            do iS = 1, 2
               ExpResDensityMatrix(iB, iS, 0, offset+trst, offset+trst) =&
                  ExpResDensityMatrix(iB, iS, 0, offset+trst, offset+trst)&
                  - DTrace%cfgsign * trval(DTrace%parttrace(trst)/DTrace%Trace)
            end do
         end do
      enddo 
   endif

end subroutine MeasExpResDensityMatrix


!===============================================================================
subroutine MeasExpOrder()
!===============================================================================
!input
   integer                       :: iB,iS
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace

 
   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   do iB=1,NBands
      do iS=1,2   
         if(DTrace%NOSOper(iB,iS)/2.lt.size(Histo(iB,iS,:))-1)then
            Histo(iB,iS,DTrace%NOSOper(iB,iS)/2)=&
               Histo(iB,iS,DTrace%NOSOper(iB,iS)/2)+1d0
         else
            Histo(iB,iS,size(Histo(iB,iS,:))-1)=&
               Histo(iB,iS,size(Histo(iB,iS,:))-1)+1d0
         endif
      enddo
   enddo
end subroutine MeasExpOrder

!===============================================================================
subroutine MeasOuterHistograms()
!===============================================================================
   integer                    :: i, offset
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   ! number of visits (as "part of configuration" where applicable)
   OuterSuperstateHisto(DTrace%outer_sst) =&
      OuterSuperstateHisto(DTrace%outer_sst) + 1d0

   if (b_statesampling) then
      OuterStateHisto(DStates%SubStates(DTrace%outer_sst)%Offset + DTrace%outer_state - 1) =&
         OuterStateHisto(DStates%SubStates(DTrace%outer_sst)%Offset + DTrace%outer_state - 1) + 1d0
   else
      offset = DStates%SubStates(DTrace%outer_sst)%Offset + DTrace%outer_state - 1
      OuterStateHisto(offset : offset + DTrace%outer_sst_size - 1) =&
         OuterStateHisto(offset : offset + DTrace%outer_sst_size - 1) + 1d0
   end if

   ! superstate signs
   SignHistoSuperstates(DTrace%outer_sst) =&
      SignHistoSuperstates(DTrace%outer_sst)&
      + DTrace%cfgsign

   ! state signs
   if (b_statesampling) then
      SignHistoStates(DStates%SubStates(DTrace%outer_sst)%Offset + DTrace%outer_state - 1) =&
         SignHistoStates(DStates%SubStates(DTrace%outer_sst)%Offset + DTrace%outer_state - 1)&
         + DTrace%cfgsign
   else
      offset = DStates%SubStates(DTrace%outer_sst)%Offset + DTrace%outer_state - 1
      SignHistoStates(offset : offset + DTrace%outer_sst_size - 1) =&
         SignHistoStates(offset : offset + DTrace%outer_sst_size - 1)&
         + DTrace%Trace%sign&
           * DTrace%parttrace(DTrace%outer_state : DTrace%outer_state + DTrace%outer_sst_size - 1)%sign * DTrace%cfgsign
   end if

   ! superstate trace contribution
   if (b_statesampling) then
      TraceContribSuperstates(DTrace%outer_sst) =&
         TraceContribSuperstates(DTrace%outer_sst)&
         + exp(DTrace%parttrace(DTrace%outer_state)%log)
   else
      TraceContribSuperstates(DTrace%outer_sst) =&
         TraceContribSuperstates(DTrace%outer_sst)&
         + sum(exp(DTrace%parttrace(1:DTrace%outer_sst_size)%log))
   end if

   ! state trace contribution
   do i = DTrace%outer_state, DTrace%outer_state + DTrace%outer_sst_size - 1
      TraceContribStates(DStates%Substates(DTrace%outer_sst)%Offset + i - 1) =&
         TraceContribStates(DStates%Substates(DTrace%outer_sst)%Offset + i - 1)&
         + exp(DTrace%parttrace(i)%log)
   end do

end subroutine MeasOuterHistograms

!===============================================================================
!> Measures the sign in Z
subroutine MeasSign()
!===============================================================================
   type(TTrace),pointer       :: DTrace
   type(TTrace_pointer)       :: pDTrace

   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   meas_sign = meas_sign + DTrace%cfgsign
   cnt_sign_z=cnt_sign_z+1     
   
end subroutine MeasSign

!===============================================================================
subroutine StepAdd(change_outer,Sector)
!===============================================================================
!input
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   logical                    :: change_outer
   integer                    :: Sector
!local
   integer                       :: CA,FPos(2),begin_sst,end_sst
   integer                       :: old_NPairs
   real(KINDR),pointer           :: temp(:,:)
   type(TLogTr)                  :: TraceNew
   real(KINDR)                   :: S,rand
   type(TLogDet)                 :: DetNew, DetRat
   real(KINDR)                   :: BosTraceNew,taudiff_factor
   real(KINDR)                   :: outer_state_factor, outertau_old
   type(TOperPointer)            :: Oper(2),begin,end
   real(KINDR),pointer           :: Q(:),R(:)
   logical                       :: adjacent,adj_qn_violate, global
   real(kindr) :: prob
   real(kindr),pointer            :: FullHybr_offdiag(:,:)
   logical :: equal_times

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   ! Fortran 90 does not guarantee nullified pointers
   nullify(temp); nullify(Q); nullify(R)
   nullify(FullHybr_offdiag)
   global = .false.

   if(DTrace%iOperPool.eq.0)then
      do CA=1,2
         allocate(Oper(CA)%p)
!        allocate(Oper(CA)%p%normr(DTrace%NTruncStatesMax))
!        allocate(Oper(CA)%p%norml(DTrace%NTruncStatesMax))
         Oper(CA)%p%calc_new=.true.
      enddo
   else
      Oper(1)%p=>DTrace%OperPool(DTrace%iOperPool)%p
      Oper(2)%p=>DTrace%OperPool(DTrace%iOperPool-1)%p
      DTrace%OperPool(DTrace%iOperPool)%p=>null()
      DTrace%OperPool(DTrace%iOperPool-1)%p=>null()
      DTrace%iOperPool=DTrace%iOperPool-2
      Oper(1)%p%calc_new=.true.
      Oper(2)%p%calc_new=.true.
      Oper(1)%p%cache_written = .false.
      Oper(2)%p%cache_written = .false.
   endif

   TryAdd(Sector)=TryAdd(Sector)+1
   if (change_outer .eqv. .true.) TryAddOuter = TryAddOuter + 1
   if (change_outer .eqv. .false.) TryAddInner = TryAddInner + 1

   if (b_statesampling .and. associated(DTrace%first)) then
      outertau_old = DTrace%beta - DTrace%last%tau + DTrace%first%tau
   else
      outertau_old = DTrace%beta
   end if
   outer_state_factor = 1_KINDR

   if (.not. b_offdiag) then
      call pair_OperAdd(DTrace,NBands,Oper,2,taudiff_factor,force_diagonal=.true.)
   else
      call pair_OperAdd(DTrace,NBands,Oper,2,taudiff_factor,force_diagonal=.false.)
   endif

   Oper(1)%p%has_hyb=.true.
   Oper(2)%p%has_hyb=.true.


   equal_times = .not. process_OperAdd(DTrace,Oper,FPos,2)

   if (equal_times) then
      do CA=1,size(Oper)
         if(associated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(associated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(associated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(associated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         Oper(CA)%p%calc_new=.true.
         Oper(CA)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+2

      AccQNAdd(Sector)=AccQNAdd(Sector)+1
      if (change_outer .eqv. .true.) AccQNAddOuter = AccQNAddOuter + 1
      if (change_outer .eqv. .false.) AccQNAddInner = AccQNAddInner + 1
      return
   endif

   ! determine substring of trace that needs quantum number checking
   adjacent = .false.
   adj_qn_violate = .false.

   if (change_outer) then
      begin%p => Oper(2)%p%prev
      if (associated(begin%p, Oper(1)%p)) then
         ! if insertions are adjacent and a change of the outer part
         ! of the sequence should be attempted, just store the sst
         ! after the last operator before the insertion ("inner sst")
         adjacent = .true.
         begin%p => begin%p%prev
         end%p => begin%p
      else
         end%p => Oper(1)%p%next
      endif
   else
      begin%p => Oper(1)%p%prev
      end%p => Oper(2)%p%next
   endif
   call get_current_ssts(DTrace, begin, end, begin_sst, end_sst)

   if (adjacent) then
      ! check sequence from Oper(2) over beta back to itself,
      ! beginning with the result of its application to the inner sst
      ! and ending with the inner sst
      begin%p => Oper(2)%p
      end%p => Oper(2)%p
      end_sst = begin_sst
      begin_sst = DStates%SubStates(end_sst)%Connect(end%p%Orbital, end%p%Spin, end%p%CA)
      if (begin_sst == -1) adj_qn_violate = .true.
   endif

   call insert_Oper(DTrace,Oper)

   call save_outer_sst(DTrace)
   !now we can check the quantum number violations
   if ((.not.check_EqualTime(DTrace,Sector)) .or. (adj_qn_violate .eqv. .true.)&
       .or. (.not. check_sst_sequence(DTrace, DStates, begin, end, begin_sst, end_sst))) then
      call restore_outer_sst(DTrace)
      !if we violate quantum numbers we need to remove operators again before exiting
      call remove_Oper(DTrace,Oper)
      ! it is not necessary to deallocate stuff here
      do CA=1,size(Oper)
         if(associated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(associated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(associated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(associated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         Oper(CA)%p%calc_new=.true.
         Oper(CA)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+2
      
      AccQNAdd(Sector)=AccQNAdd(Sector)+1
      if (change_outer .eqv. .true.) AccQNAddOuter = AccQNAddOuter + 1
      if (change_outer .eqv. .false.) AccQNAddInner = AccQNAddInner + 1
      return
   endif

   if (b_statesampling .and. change_outer) then
      global = .true.
      outer_state_factor = WeightedOuterStateChoice(DTrace, DStates,&
                              outertau_old,&
                              DTrace%beta - DTrace%last%tau + DTrace%first%tau,&
                              DTrace%outer_sst_old,&
                              DTrace%outer_state_old,&
                              (/ DTrace%outer_sst_old /),&
                              (/ DTrace%outer_sst /))
   end if

   if(b_offdiag)then

      if(b_full_offdiag)then

         FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
         call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag, DetNew)
         DetRat = DetNew/DTrace%Det

      else

         !!! fast updates
         call add_ftau_offset(DTrace, DStates, Oper, FPos)
         DetRat = get_LogDetRatPairAdd_full(dtrace,oper,DStates,Q,R,S,FPos)
         DetNew = DetRat * DTrace%Det

      endif

   else

      DetRat = get_LogDetRatPairAdd(DTrace,DTrace%MInv,Oper,Q,R,S,FPos)
      DetNew = DetRat * DTrace%Det

   endif

   TraceNew = get_trace_EB(DTrace,DStates,global=global)
   BosTraceNew = get_BosonicTrace(DTrace,DStates)
   rand=grnd()

   old_NPairs = DTrace%NPairs
   if (old_NPairs == -1) then
      call ensure_valid_NPairs(DTrace)
   else
      call update_NPairs(DTrace, Oper, 1)
   end if

   if (b_offdiag) then
      prob = abs(outer_state_factor * taudiff_factor * (2 * NBands)**2/real(DTrace%NPairs, KINDR)&
                 *wrat(TraceNew/DTrace%Trace, DetRat)*BosTraceNew/DTrace%BosonicTrace)
   else
      prob = abs(outer_state_factor * taudiff_factor * (2 * NBands)/real(DTrace%NPairs, KINDR)&
                 *wrat(TraceNew/DTrace%Trace, DetRat)*BosTraceNew/DTrace%BosonicTrace)
   end if
   
   if(rand.lt.prob)then

      if(.not.b_offdiag)then

         temp=>get_MatAdd(&
         DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat,&
         DTrace%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)/2-1,Q,R,S,FPos)
         if(associated(DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat))then
            deallocate(DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat)
         endif
         DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat=>temp

      else

         if(b_full_offdiag)then

            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>FullHybr_offdiag

         else

            if (allocated(DTrace%wormContainer)) then
               temp=>get_MatAdd(&
               DTrace%MInv_full%Mat,&
               (DTrace%Noper-size(DTrace%wormContainer))/2-1,Q,R,S,FPos)
            else
               temp=>get_MatAdd(&
               DTrace%MInv_full%Mat,&
               DTrace%Noper/2-1,Q,R,S,FPos)
            end if


            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>temp

         endif

      endif
            
      !we save the determinant of the current configuration
      DTrace%Det = DetNew
      DTrace%Trace = TraceNew
      DTrace%BosonicTrace = BosTraceNew
      
      !calculate sign of current configuration
      DTrace%cfgsign = -get_Sign(DTrace, DStates)

      call update_trace_EB(DTrace, global=global)
      
      AccAdd(Sector)=AccAdd(Sector)+1
      AccPair(Oper(1)%p%Orbital, Oper(1)%p%Spin,&
         min(NAccPair, floor((Oper(2)%p%tau - Oper(1)%p%tau)/AccPairMax * dble(NAccPair)) + 1)) =&
         AccPair(Oper(1)%p%Orbital, Oper(1)%p%Spin,&
         min(NAccPair, floor((Oper(2)%p%tau - Oper(1)%p%tau)/AccPairMax * dble(NAccPair)) + 1)) + 1
      if (allocated(AccPairTau)) then
         if (apt_index <= size(AccPairTau)) then
            AccPairTau(apt_index) = (Oper(2)%p%tau - Oper(1)%p%tau)
            apt_index = apt_index + 1
         end if
      end if
      if (change_outer .eqv. .true.) AccAddOuter = AccAddOuter + 1
      if (change_outer .eqv. .false.) AccAddInner = AccAddInner + 1
      
   else
      call restore_outer_sst(DTrace)
      call remove_Oper(DTrace,Oper)

      if (old_NPairs /= -1) then
         DTrace%NPairs = old_NPairs
      else
         call update_NPairs(DTrace, Oper, -1)
      end if

      do CA=1,size(Oper)
         if(associated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(associated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(associated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(associated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         Oper(CA)%p%calc_new=.true.
         Oper(CA)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+2

      if (b_offdiag .and. b_full_offdiag) deallocate(FullHybr_offdiag)

   endif
   
   if (.not. (b_offdiag .and. b_full_offdiag)) deallocate(Q,R)

end subroutine StepAdd

!===============================================================================
subroutine StepAdd4(Sector)
!===============================================================================
!input
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   integer                    :: Sector
!local
   integer                       :: CA,FPos(2),iB,iS
   integer                       :: begin_sst, end_sst, old_NPairs
   type(TLogTr)                  :: TraceNew
   real(KINDR)                   :: rand
   type(TLogDet)                 :: DetNew, DetRat
   real(KINDR)                   :: BosTraceNew,taudiff_factor
   real(KINDR)                   :: rempair_factor
   integer                       :: NPairs_only12, NPairs_only34
   integer                       :: NPairs_only14, NPairs_only23
   integer                       :: NPairs_all
   type(TOperPointer)            :: Oper(4), begin, end, oplist(2)
   logical                       :: equal_times
   real(kindr) :: prob
   type(TSubMatrix),pointer      :: FullHybr(:,:)
   real(kindr),pointer           :: FullHybr_offdiag(:,:)


   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   ! Fortran 90 does not guarantee nullified pointers
   nullify(FullHybr)

   if(DTrace%iOperPool<4)then
      do CA=1,size(Oper)
         allocate(Oper(CA)%p)
         Oper(CA)%p%calc_new=.true.
      enddo
   else
      do CA=1,size(Oper)
         Oper(CA)%p=>DTrace%OperPool(DTrace%iOperPool)%p
         DTrace%OperPool(DTrace%iOperPool)%p=>null()
         DTrace%iOperPool=DTrace%iOperPool-1
         Oper(CA)%p%calc_new=.true.
         Oper(CA)%p%cache_written = .false.
      enddo
   endif

   TryAdd4=TryAdd4+1

   if (b_offdiag) then
      call pair_OperAdd(DTrace,NBands,Oper,4,taudiff_factor,force_diagonal=.false.)
   else
      call pair_OperAdd(DTrace,NBands,Oper,4,taudiff_factor,force_diagonal=.true.)
   end if

   Oper(1)%p%has_hyb=.true.
   Oper(2)%p%has_hyb=.true.
   Oper(3)%p%has_hyb=.true.
   Oper(4)%p%has_hyb=.true.

   call ensure_valid_NPairs(DTrace)
   old_NPairs = DTrace%NPairs

   ! sort the generated operators to the correct position in the
   ! trace, insert them and count the number of removable pairs with
   ! each new pair individually and both at once added
   equal_times = .not. process_OperAdd(DTrace,Oper(1:2),FPos,2)
   call insert_Oper(DTrace, Oper(1:2))
   call update_NPairs(DTrace, Oper(1:2), 1)
   NPairs_only12 = DTrace%NPairs

   equal_times = equal_times .or. .not. process_OperAdd(DTrace,Oper(3:4),FPos,2)
   call insert_Oper(DTrace, Oper(3:4))
   call remove_Oper(DTrace, Oper(1:2))
   DTrace%NPairs = old_NPairs
   call update_NPairs(DTrace, Oper(3:4), 1)
   NPairs_only34 = DTrace%NPairs

   call insert_Oper(DTrace, Oper(1:2))
   call update_NPairs(DTrace, Oper(1:2), 1)

   NPairs_all = DTrace%NPairs

   if (is_rempair(DTrace, Oper(1)%p, Oper(4)%p)&
       .and. is_rempair(DTrace, Oper(2)%p, Oper(3)%p)) then

      ! in this case, there are four different choices of pairs (every
      ! choice of one of the inserted creators and one of the inserted
      ! annihilators is removable) that would all remove the four
      ! operators, so the probability of the removal proposal is a sum
      ! of four probabilities depending on the first removed pair;
      ! removable pair counts for the two other possibly remaining
      ! pairs need to be calculated as well
      if (Oper(1)%p%tau < Oper(4)%p%tau) then
         oplist = (/Oper(1), Oper(4)/)
      else
         oplist = (/Oper(4), Oper(1)/)
      end if
      call remove_Oper(DTrace, oplist)
      call update_NPairs(DTrace, oplist, -1)
      NPairs_only23 = DTrace%NPairs
      call insert_Oper(DTrace, oplist)
      DTrace%NPairs = NPairs_all

      if (Oper(2)%p%tau < Oper(3)%p%tau) then
         oplist = (/Oper(2), Oper(3)/)
      else
         oplist = (/Oper(3), Oper(2)/)
      end if
      call remove_Oper(DTrace, oplist)
      call update_NPairs(DTrace, oplist, -1)
      NPairs_only14 = DTrace%NPairs
      call insert_Oper(DTrace, oplist)
      DTrace%NPairs = NPairs_all

      rempair_factor = (1_KINDR/real(NPairs_all, KINDR))&
         * (1_KINDR/real(NPairs_only12, KINDR)&
            + 1_KINDR/real(NPairs_only34, KINDR)&
            + 1_KINDR/real(NPairs_only14, KINDR)&
            + 1_KINDR/real(NPairs_only23, KINDR))

      ! the insertion tau-choice factor depends on the position of the
      ! annihilator of each pair relative to the window and especially
      ! not on the order of the pairs, so it is the same for all four
      ! possibilities; taudiff_factor is the inverse proposal
      ! probability
      taudiff_factor = taudiff_factor / 4_KINDR
   else
      rempair_factor = (1_KINDR/real(NPairs_all, KINDR))&
         * (1_KINDR/real(NPairs_only12, KINDR)&
            + 1_KINDR/real(NPairs_only34, KINDR))

      taudiff_factor = taudiff_factor / 2_KINDR
   end if

   if (Oper(1)%p%tau < Oper(3)%p%tau) then
      begin%p => Oper(1)%p%prev
   else
      begin%p => Oper(3)%p%prev
   end if
   if (Oper(2)%p%tau > Oper(4)%p%tau) then
      end%p => Oper(2)%p%next
   else
      end%p => Oper(4)%p%next
   end if
   call get_current_ssts(DTrace, begin, end, begin_sst, end_sst)

   call save_outer_sst(DTrace)
   !now we can check the quantum number violations
   !TODO: equal_times considers regular operators,
   !check_EqualTime considers worm operators (naming conflict?)
   if (equal_times .or. (.not.check_EqualTime(DTrace,Sector))&
       .or. (.not. check_sst_sequence(DTrace, DStates, begin, end, begin_sst, end_sst))) then
      call restore_outer_sst(DTrace)
      !if we violate quantum numbers we need to remove operators again before exiting
      call remove_Oper(DTrace,Oper(1:2))
      call remove_Oper(DTrace,Oper(3:4))
      DTrace%NPairs = old_NPairs
      do CA=1,size(Oper)
         if(associated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(associated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(associated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(associated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         Oper(CA)%p%calc_new=.true.
         Oper(CA)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+4
      
      !AccQNAdd(Sector)=AccQNAdd(Sector)+1
      return
   endif

   if(b_offdiag)then
       
      FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
      call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag, DetNew)
      DetRat = DetNew/DTrace%Det

   else

      FullHybr => get_InvFullHybr(DTrace,DStates,DetNew)
      DetRat = DetNew*get_FullLogDet(DTrace,DStates)

   endif

   TraceNew = get_trace_EB(DTrace,DStates)
   BosTraceNew = get_BosonicTrace(DTrace,DStates)
   rand=grnd()

   prob = abs(taudiff_factor * (2 * NBands)**4 * rempair_factor&
              *wrat(TraceNew/DTrace%Trace, DetRat)*BosTraceNew/DTrace%BosonicTrace)

   if(rand.lt.prob)then

      if(b_offdiag)then

         if(associated(DTrace%MInv_full%Mat))then
            deallocate(DTrace%MInv_full%Mat)
         endif
         DTrace%MInv_full%Mat=>FullHybr_offdiag

      else

         do iB=1,NBands
         do iS=1,2
            deallocate(DTrace%MInv(iB,iS)%Mat)
         enddo
         enddo
         deallocate(DTrace%MInv)
         DTrace%MInv=>FullHybr

      endif
            
      !we save the determinant of the current configuration
      DTrace%Det = DetNew
      DTrace%Trace = TraceNew
      DTrace%BosonicTrace = BosTraceNew
      
      !calculate sign of current configuration
      DTrace%cfgsign = -get_Sign(DTrace, DStates)
      
      call update_trace_EB(DTrace)

      AccAdd4=AccAdd4+1
      
   else
      call restore_outer_sst(DTrace)
      call remove_Oper(DTrace,Oper(1:2))
      call remove_Oper(DTrace,Oper(3:4))
      DTrace%NPairs = old_NPairs
      do CA=1,size(Oper)
         if(associated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(associated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(associated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(associated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         Oper(CA)%p%calc_new=.true.
         Oper(CA)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+4

      if(b_offdiag)then

         deallocate(FullHybr_offdiag)

      else

         do iB=1,NBands
         do iS=1,2
            deallocate(FullHybr(iB,iS)%Mat)
         enddo
         enddo
         deallocate(FullHybr)

      endif

   endif
   
end subroutine StepAdd4

!===============================================================================
subroutine StepRem(change_outer,Sector)
!===============================================================================
!input
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   logical                    :: change_outer
   integer                    :: Sector
!local
   integer                       :: CA,FPos(2),begin_sst,end_sst
   real(KINDR),pointer           :: temp(:,:)
   type(TLogTr)                  :: TraceNew
   real(KINDR)                   :: S,rand
   type(TLogDet)                 :: DetNew, DetRat
   real(KINDR)                   :: BosTraceNew,hybpairremfactor
   real(KINDR)                   :: outer_state_factor,outertau_old,outertau_new
   type(TOperPointer)            :: Oper(2),begin,end
   type(TOper),pointer           :: El
   logical                       :: adjacent,adj_qn_violate,global
   real(kindr) :: prob
   real(kindr),pointer            :: FullHybr_offdiag(:,:)

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   nullify(temp)
   nullify(FullHybr_offdiag)
   global = .false.

   TryRem(Sector)=TryRem(Sector)+1
   if (change_outer .eqv. .true.) TryRemOuter = TryRemOuter + 1
   if (change_outer .eqv. .false.) TryRemInner = TryRemInner + 1

   if(.not.gen_OperRemove(DTrace,Oper,FPos,2,hybpairremfactor,with_hyb=.true.))then
      AccQNRem(Sector)=AccQNRem(Sector)+1
      if (change_outer .eqv. .true.) AccQNRemOuter = AccQNRemOuter + 1
      if (change_outer .eqv. .false.) AccQNRemInner = AccQNRemInner + 1
      return
   endif

   outer_state_factor = 1_KINDR
   ! if there were not >= 2 operators, gen_OperRemove would have returned false
   outertau_old = DTrace%beta - DTrace%last%tau + DTrace%first%tau

   ! determine substring of trace that needs quantum number checking
   adjacent = .false.
   adj_qn_violate = .false.

   if (change_outer) then
      begin%p => Oper(2)%p%prev
      if (associated(begin%p, Oper(1)%p)) then
         ! if removals are adjacent and a change of the outer part of
         ! the trace should be attempted, just store the sst after the
         ! first removed operator ("inner sst")
         adjacent = .true.
         end%p => begin%p
      else
         end%p => Oper(1)%p%next
      endif
   else
      begin%p => Oper(1)%p%prev
      end%p => Oper(2)%p%next
   endif
   call get_current_ssts(DTrace, begin, end, begin_sst, end_sst)

   if (adjacent) then
      ! check sequence from Oper(2)%next over beta back to itself,
      ! beginning with the result of its application to the inner sst
      ! and ending with the inner sst
      begin%p => Oper(2)%p%next
      end%p => Oper(2)%p%next
      end_sst = begin_sst
      if (associated(end%p)) begin_sst = DStates%SubStates(end_sst)&
                                                %Connect(end%p%Orbital, end%p%Spin, end%p%CA)
      if (begin_sst == -1) adj_qn_violate = .true.
   endif

   call remove_Oper(DTrace,Oper)

   call save_outer_sst(DTrace)
   !now we can check the quantum number violations
   if ((adj_qn_violate .eqv. .true.)&
       .or. (.not. check_sst_sequence(DTrace, DStates, begin, end, begin_sst, end_sst))) then
      call restore_outer_sst(DTrace)
      !if we violate quantum numbers we need to insert operators again before exiting
      call insert_Oper(DTrace,Oper)      
      AccQNRem(Sector)=AccQNRem(Sector)+1
      if (change_outer .eqv. .true.) AccQNRemOuter = AccQNRemOuter + 1
      if (change_outer .eqv. .false.) AccQNRemInner = AccQNRemInner + 1
      return
   endif

   if (b_statesampling .and. change_outer) then
      global = .true.

      if (associated(DTrace%first)) then
         outertau_new = DTrace%beta - DTrace%last%tau + DTrace%first%tau
      else
         outertau_new = DTrace%beta
      end if

      outer_state_factor = WeightedOuterStateChoice(DTrace, DStates,&
                              outertau_old,&
                              outertau_new,&
                              DTrace%outer_sst_old,&
                              DTrace%outer_state_old,&
                              (/ DTrace%outer_sst_old /),&
                              (/ DTrace%outer_sst /))
   end if

   if (associated(Oper(2)%p%prev)) then
      Oper(2)%p%prev%calc_new = .true.
   end if

   if (associated(Oper(1)%p%next)&
       .and. .not. associated(Oper(1)%p%next, Oper(2)%p)) then

      Oper(1)%p%next%calc_new = .true.

   else if (associated(Oper(1)%p%next, Oper(2)%p)) then

      if (associated(Oper(2)%p%next)) then

         Oper(2)%p%next%calc_new = .true.

      else if (associated(Oper(1)%p%prev)) then

         Oper(1)%p%prev%calc_new = .true.

      end if

   else
      stop "unreachable in StepRem"
   end if


   if(.not.b_offdiag)then

      DetRat = get_LogDetRatPairRem(DTrace,Oper,FPos,S)
      DetNew = DetRat * DTrace%Det

   else

      if(b_full_offdiag)then

         FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
         call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag, DetNew)
         DetRat = DetNew/DTrace%Det

      else

         !!! fast updates
         call add_ftau_offset(dtrace,dstates,oper,fpos)
         DetRat = get_LogDetRatPairRem_full(DTrace,Oper,FPos,S)
         DetNew = DetRat * DTrace%Det

      endif

   endif

   TraceNew = get_trace_EB(DTrace,DStates,global=global)
   BosTraceNew = get_BosonicTrace(DTrace,DStates)

   if (b_offdiag) then
      prob = abs(outer_state_factor * hybpairremfactor/(2 * NBands)**2&
                 *wrat(TraceNew/DTrace%Trace, DetRat)*BosTraceNew/DTrace%BosonicTrace)
   else
      prob = abs(outer_state_factor * hybpairremfactor/(2 * NBands)&
                 *wrat(TraceNew/DTrace%Trace, DetRat)*BosTraceNew/DTrace%BosonicTrace)
   end if
   
   rand=grnd()
   if(rand.lt.prob)then
      
      !update hybridization matrix
      if(dabs(detval(DetRat)).lt.1d-12) stop "DetRat too small in removal"
      
      if(.not.b_offdiag)then

         temp=>get_MatRem(&
         DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat,&
               DTrace%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)/2+1,FPos)
         deallocate(DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat)
         DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat=>temp
      else

         if(b_full_offdiag)then

            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>FullHybr_offdiag

         else

            !!! fast matrix generation
            if (allocated(DTrace%wormContainer)) then
               temp=>get_MatRem(&
                  DTrace%MInv_full%Mat,&
                  (DTrace%Noper-size(DTrace%wormContainer))/2+1,FPos)
            else
               temp=>get_MatRem(&
                  DTrace%MInv_full%Mat,&
                  DTrace%Noper/2+1,FPos)
            end if

            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            !deallocate(temp)
            DTrace%MInv_full%Mat=>temp

         endif

      endif

      
      !we save the determinant of the current configuration
      DTrace%Det = DetNew
      DTrace%Trace = TraceNew
      DTrace%BosonicTrace = BosTraceNew
      
      DTrace%cfgsign = -get_Sign(DTrace, DStates)

      call update_trace_EB(DTrace, global=global)
      call update_NPairs(DTrace, Oper, -1)

      if (Oper(1)%p%CA == 1) then
         El => Oper(1)%p
      else
         El => Oper(2)%p
      end if
      AccPair(El%Orbital, El%Spin,&
         min(NAccPair, floor((Oper(2)%p%tau - Oper(1)%p%tau)/AccPairMax * dble(NAccPair)) + 1)) =&
         AccPair(El%Orbital, El%Spin,&
         min(NAccPair, floor((Oper(2)%p%tau - Oper(1)%p%tau)/AccPairMax * dble(NAccPair)) + 1)) + 1
      if (allocated(AccPairTau)) then
         if (apt_index <= size(AccPairTau)) then
            AccPairTau(apt_index) = (Oper(2)%p%tau - Oper(1)%p%tau)
            apt_index = apt_index + 1
         end if
      end if

      do CA=1,size(Oper)
         if(associated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(associated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(associated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(associated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
         Oper(CA)%p%calc_new=.true.
         Oper(CA)%p%cache_written = .false.
      enddo
      
      DTrace%iOperPool=DTrace%iOperPool+2
      
      AccRem(Sector)=AccRem(Sector)+1
      if (change_outer .eqv. .true.) AccRemOuter = AccRemOuter + 1
      if (change_outer .eqv. .false.) AccRemInner = AccRemInner + 1
      
   else
      call restore_outer_sst(DTrace)
      call insert_Oper(DTrace,Oper)

      if (associated(Oper(1)%p%prev)) Oper(1)%p%prev%calc_new = .false.
      if (associated(Oper(2)%p%prev)) Oper(2)%p%prev%calc_new = .false.
      if (associated(Oper(1)%p%next)) Oper(1)%p%next%calc_new = .false.
      if (associated(Oper(2)%p%next)) Oper(2)%p%next%calc_new = .false.

      if (b_offdiag .and. b_full_offdiag) deallocate(FullHybr_offdiag)
   endif

end subroutine StepRem

!===============================================================================
subroutine StepRem4()
!===============================================================================
!input
   type(TStates),pointer         :: DStates
   type(TTrace),pointer          :: DTrace
   type(TStates_pointer)         :: pDStates
   type(TTrace_pointer)          :: pDTrace
!local
   integer                       :: CA,FPos(2),iB,iS
   integer                       :: begin_sst, end_sst, old_NPairs
   integer                       :: NPairs_only12, NPairs_only34
   integer                       :: NPairs_only14, NPairs_only23
   type(TLogTr)                  :: TraceNew
   real(KINDR)                   :: rand
   type(TLogDet)                 :: DetNew, DetRat
   real(KINDR)                   :: BosTraceNew,rempair_factor
   real(KINDR)                   :: taudiff_factor
   type(TOperPointer)            :: Oper(4), begin, end, oplist(2)
   real(kindr) :: prob
   type(TSubMatrix),pointer      :: FullHybr(:,:)
   real(kindr),pointer           :: FullHybr_offdiag(:,:)

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr
 
   nullify(FullHybr)
   nullify(FullHybr_offdiag)

   TryRem4=TryRem4+1

   ! the actual taudiff_factor will be extracted from values written
   ! by gen_OperRemove and the actual rempair_factor will later be
   ! calculated manually, here I just use them as write targets
   if(.not.gen_OperRemove(DTrace,Oper(1:2),FPos,2,rempair_factor,with_hyb=.true.))then
      !AccQNRem(Sector)=AccQNRem(Sector)+1
      return
   endif
   old_NPairs = DTrace%NPairs
   call remove_Oper(DTrace, Oper(1:2))
   call update_NPairs(DTrace, Oper(1:2), -1)
   NPairs_only34 = DTrace%NPairs

   if(.not.gen_OperRemove(DTrace,Oper(3:4),FPos,2,taudiff_factor,with_hyb=.true.))then
      DTrace%NPairs = old_NPairs
      call insert_Oper(DTrace, Oper(1:2))
      !AccQNRem(Sector)=AccQNRem(Sector)+1
      return
   endif
   taudiff_factor = real(old_NPairs, KINDR) * real(NPairs_only34, KINDR)&
                    / rempair_factor / taudiff_factor
   call insert_Oper(DTrace, Oper(1:2))
   DTrace%NPairs = old_NPairs
   call remove_Oper(DTrace, Oper(3:4))
   call update_NPairs(DTrace, Oper(3:4), -1)
   NPairs_only12 = DTrace%NPairs

   if (is_rempair(DTrace, Oper(1)%p, Oper(4)%p)&
       .and. is_rempair(DTrace, Oper(2)%p, Oper(3)%p)) then

      ! see StepAdd4 for comments
      call insert_Oper(DTrace, Oper(3:4))
      DTrace%NPairs = old_NPairs

      if (Oper(1)%p%tau < Oper(4)%p%tau) then
         oplist = (/Oper(1), Oper(4)/)
      else
         oplist = (/Oper(4), Oper(1)/)
      end if
      call remove_Oper(DTrace, oplist)
      call update_NPairs(DTrace, oplist, -1)
      NPairs_only23 = DTrace%NPairs
      call insert_Oper(DTrace, oplist)
      DTrace%NPairs = old_NPairs

      if (Oper(2)%p%tau < Oper(3)%p%tau) then
         oplist = (/Oper(2), Oper(3)/)
      else
         oplist = (/Oper(3), Oper(2)/)
      end if
      call remove_Oper(DTrace, oplist)
      call update_NPairs(DTrace, oplist, -1)
      NPairs_only14 = DTrace%NPairs
      call insert_Oper(DTrace, oplist)

      call remove_Oper(DTrace, Oper(3:4))
      DTrace%NPairs = NPairs_only12

      rempair_factor = (1_KINDR/real(old_NPairs, KINDR))&
         * (1_KINDR/real(NPairs_only12, KINDR)&
            + 1_KINDR/real(NPairs_only34, KINDR)&
            + 1_KINDR/real(NPairs_only14, KINDR)&
            + 1_KINDR/real(NPairs_only23, KINDR))

      taudiff_factor = taudiff_factor / 4_KINDR
   else
      rempair_factor = (1_KINDR/real(old_NPairs, KINDR))&
         * (1_KINDR/real(NPairs_only12, KINDR)&
            + 1_KINDR/real(NPairs_only34, KINDR))

      taudiff_factor = taudiff_factor / 2_KINDR
   end if

   call remove_Oper(DTrace,Oper(1:2))
   call update_NPairs(DTrace, Oper(1:2), -1)

   if (Oper(1)%p%tau < Oper(3)%p%tau) then
      begin%p => Oper(1)%p%prev
   else
      begin%p => Oper(3)%p%prev
   end if
   if (Oper(2)%p%tau > Oper(4)%p%tau) then
      end%p => Oper(2)%p%next
   else
      end%p => Oper(4)%p%next
   end if
   call get_current_ssts(DTrace, begin, end, begin_sst, end_sst)

   call save_outer_sst(DTrace)
   if (.not. check_sst_sequence(DTrace, DStates, begin, end, begin_sst, end_sst)) then
      call restore_outer_sst(DTrace)
      DTrace%NPairs = old_NPairs
      !if we violate quantum numbers we need to insert operators again before exiting
      call insert_Oper(DTrace,Oper(1:2))
      call insert_Oper(DTrace,Oper(3:4))
      !AccQNRem(Sector)=AccQNRem(Sector)+1
      return
   endif

   if (associated(begin%p)) then
      if (associated(begin%p%next)) then
         begin%p%next%calc_new = .true.
      else
         begin%p%calc_new = .true.
      end if
   else
      if (associated(DTrace%first)) DTrace%first%calc_new = .true.
   end if
   if (associated(end%p)) then
      if (associated(end%p%prev)) then
         end%p%prev%calc_new = .true.
      else
         end%p%calc_new = .true.
      end if
   else
      if (associated(DTrace%last)) DTrace%last%calc_new = .true.
   end if

   if(b_offdiag)then

      FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
      call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag, DetNew)
      DetRat = DetNew/DTrace%Det

   else

      FullHybr => get_InvFullHybr(DTrace,DStates,DetNew)
      DetRat = DetNew*get_FullLogDet(DTrace,DStates)

   endif
   
   TraceNew = get_trace_EB(DTrace,DStates)
   BosTraceNew = get_BosonicTrace(DTrace,DStates)

   prob = abs(1_KINDR / real((2 * NBands)**4, KINDR) / taudiff_factor / rempair_factor&
              *wrat(TraceNew/DTrace%Trace, DetRat)*BosTraceNew/DTrace%BosonicTrace)

   rand=grnd()
   if(rand.lt.prob)then

      if(b_offdiag)then

         if(associated(DTrace%MInv_full%Mat))then
            deallocate(DTrace%MInv_full%Mat)
         endif
         DTrace%MInv_full%Mat=>FullHybr_offdiag
         
      else

         do iB=1,NBands
         do iS=1,2
            deallocate(DTrace%MInv(iB,iS)%Mat)
         enddo
         enddo
         deallocate(DTrace%MInv)
         DTrace%MInv=>FullHybr

      endif
      
      !we save the determinant of the current configuration
      DTrace%Det = DetNew
      DTrace%Trace = TraceNew
      DTrace%BosonicTrace = BosTraceNew
      
      DTrace%cfgsign = -get_Sign(DTrace, DStates)
      
      call update_trace_EB(DTrace)
      do CA=1,size(Oper)
         if(associated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(associated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(associated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(associated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         Oper(CA)%p%calc_new=.true.
         Oper(CA)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
      enddo
      
      DTrace%iOperPool=DTrace%iOperPool+4
      
      AccRem4=AccRem4+1
      
   else
      call restore_outer_sst(DTrace)
      call insert_Oper(DTrace,Oper(1:2))
      call insert_Oper(DTrace,Oper(3:4))
      DTrace%NPairs = old_NPairs

      if (associated(Oper(1)%p%prev)) Oper(1)%p%prev%calc_new = .false.
      if (associated(Oper(2)%p%prev)) Oper(2)%p%prev%calc_new = .false.
      if (associated(Oper(1)%p%next)) Oper(1)%p%next%calc_new = .false.
      if (associated(Oper(2)%p%next)) Oper(2)%p%next%calc_new = .false.
      if (associated(Oper(3)%p%prev)) Oper(3)%p%prev%calc_new = .false.
      if (associated(Oper(4)%p%prev)) Oper(4)%p%prev%calc_new = .false.
      if (associated(Oper(3)%p%next)) Oper(3)%p%next%calc_new = .false.
      if (associated(Oper(4)%p%next)) Oper(4)%p%next%calc_new = .false.
      if (associated(DTrace%first)) DTrace%first%calc_new = .false.
      if (associated(DTrace%last)) DTrace%last%calc_new = .false.

      if(b_offdiag)then

         deallocate(FullHybr_offdiag)

      else

         do iB=1,NBands
         do iS=1,2
            deallocate(FullHybr(iB,iS)%Mat)
         enddo
         enddo
         deallocate(FullHybr)

      endif

   endif
   
end subroutine StepRem4

!===============================================================================
subroutine StepGlob()
!===============================================================================
!input
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
!local 
   type(TLogTr)                  :: TraceNew
   real(KINDR)                   :: rand,sst_factor
   type(TLogDet)                 :: DetNew, DetRat
   real(KINDR)                   :: BosTraceNew
   type(TSubMatrix),pointer      :: FullHybr(:,:)
   integer                       :: iB,iS,move_type
   type(TOper),pointer           :: Element
   logical                       :: Debg
   real(kindr),pointer           :: FullHybr_offdiag(:,:)

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   nullify(FullHybr); nullify(Element)
   nullify(FullHybr_offdiag)

   TryGlob=TryGlob+1
   if(DTrace%NOper.eq.0) then
      return
   endif
   
! adding spin flip for paramagnetic case
   rand=grnd()
   if(rand.lt.4d-1)then
      move_type = 0
      call gen_SpinFlipUpdate(DTrace,DStates)
   elseif(rand.lt.7d-1.and.get_Integer_Parameter("NSymMove").gt.0)then
      move_type = 1
      call gen_SymUpdate(DTrace,DStates)
   else
      move_type = 2
      !!! in the original global moves, permutations were proposed not uniformly
      !Debg=gen_GlobalUpdate(DTrace,DStates)
      Debg=gen_GlobalUpdate_new(DTrace,DStates)
      globalmove_check=globalmove_check+dtrace%gu(1:)
   endif
   call save_outer_sst(DTrace)

   if (move_type == 0) then
      if (index(get_String_Parameter("QuantumNumbers"), "Szt") /= 0) then
         sst_factor = check_spinflip(DTrace, DStates)
      else
         sst_factor = check_qn_global(DTrace, DStates)
      endif
   else
      sst_factor = check_qn_global(DTrace, DStates)
   endif

   if (sst_factor == 0._KINDR) then
      call restore_outer_sst(DTrace)
      call gen_InverseGlobalUpdate(DTrace,DStates)
      call globalUpdate(DTrace,DStates)
      AccQNGlob=AccQNGlob+1
      return
   endif
   
   if(b_offdiag)then
      FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
      call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag, DetNew)
      DetRat = DetNew/DTrace%Det
   else
      FullHybr=>get_InvFullHybr(DTrace,DStates,DetNew)
      DetRat=DetNew*get_FullLogDet(DTrace,DStates)
   endif

   TraceNew = get_trace_EB(DTrace,DStates,global=.true.)
   BosTraceNew = get_BosonicTrace(DTrace,DStates)

   if(grnd().lt.abs(wrat(TraceNew/DTrace%Trace, DetRat)*BosTraceNew/DTrace%BosonicTrace&
                    *sst_factor)) then
      DTrace%Det=DetNew
      DTrace%Trace=TraceNew
      DTrace%BosonicTrace=BosTraceNew
      DTrace%cfgsign = -get_Sign(DTrace, DStates)

      if(b_offdiag)then

         if(associated(DTrace%MInv_full%Mat))then
            deallocate(DTrace%MInv_full%Mat)
         endif
         DTrace%MInv_full%Mat=>FullHybr_offdiag

      else

         do iB=1,NBands
         do iS=1,2
            deallocate(DTrace%MInv(iB,iS)%Mat)
         enddo
         enddo
         deallocate(DTrace%MInv)
         DTrace%MInv=>FullHybr

      endif

      AccGlob=AccGlob+1

      call update_trace_EB(DTrace, global=.true.)
   else
      call restore_outer_sst(DTrace)
      call gen_InverseGlobalUpdate(DTrace,DStates)
      call globalUpdate(DTrace,DStates)

      if(b_offdiag)then

         deallocate(FullHybr_offdiag)

      else

         do iB=1,NBands
         do iS=1,2
            deallocate(FullHybr(iB,iS)%Mat)
         enddo
         enddo
         deallocate(FullHybr)

      endif

   endif

end subroutine StepGlob

!===============================================================================
subroutine StepShiftTau()
!===============================================================================
!input
   type(TTrace),pointer          :: DTrace
   type(TTrace_pointer)          :: pDTrace

   pDTrace=transfer(ipDTrace,pDTrace)
   DTrace => pDTrace%ptr

   if (DTrace%NOper == 0) then
      call StepChangeOuter()
   else
      call StepShiftTau_proper()
   end if
end subroutine StepShiftTau

!===============================================================================
subroutine StepChangeOuter()
!===============================================================================
!input
   type(TStates),pointer         :: DStates
   type(TTrace),pointer          :: DTrace
   type(TStates_pointer)         :: pDStates
   type(TTrace_pointer)          :: pDTrace

   real(KINDR), allocatable      :: empty_trace_values(:)
   integer                       :: iSt, i, j
   real(KINDR)                   :: accumulator, rand

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   call save_outer_sst(DTrace)
1  if (b_statesampling) then
      ! choose new outer state with probability proportional to the empty trace value
      ! due to this proposal weight, the Metropolis acceptance probability is 1
      allocate(empty_trace_values(DTrace%NTruncStates))
      iSt = 1
      do i = 1, size(DTrace%States(:, 1))
         do j = 1, DTrace%States(i, 2)
            DTrace%outer_sst = DTrace%States(i, 1)
            DTrace%outer_state = j
            empty_trace_values(iSt) = trval(get_Trace_EB(DTrace, DStates, global=.true.))
            iSt = iSt + 1
         end do
      end do
      empty_trace_values = empty_trace_values / sum(empty_trace_values)

      rand = grnd()
      accumulator = 0._KINDR
      iSt = 1
      sstloop: do i = 1, size(DTrace%States(:, 1))
         do j = 1, DTrace%States(i, 2)
            DTrace%outer_sst = DTrace%States(i, 1)
            DTrace%outer_state = j
            accumulator = accumulator + empty_trace_values(iSt)
            if (accumulator >= rand .and. (empty_trace_values(iSt) /= 0._KINDR))&
               exit sstloop
            iSt = iSt + 1
         end do
      end do sstloop
      deallocate(empty_trace_values)
   else
      ! choose new outer superstate with probability proportional to the empty trace value
      ! due to this proposal weight, the Metropolis acceptance probability is 1
      allocate(empty_trace_values(size(DTrace%States(:, 1))))
      do i = 1, size(DTrace%States(:, 1))
         DTrace%outer_sst = DTrace%States(i, 1)
         empty_trace_values(i) = trval(get_Trace_EB(DTrace, DStates, global=.true.))
      end do
      empty_trace_values = empty_trace_values / sum(empty_trace_values)

      rand = grnd()
      accumulator = 0._KINDR
      do i = 1, size(DTrace%States(:, 1))
         DTrace%outer_sst = DTrace%States(i, 1)
         accumulator = accumulator + empty_trace_values(i)
         if (accumulator >= rand .and. (empty_trace_values(i) /= 0._KINDR))&
            exit
      end do
      deallocate(empty_trace_values)
   end if

   DTrace%Trace = get_Trace_EB(DTrace, DStates, global=.true.)
   ! the trace should never be 0, but just in case it is
   if (DTrace%Trace%log < -huge(0.0_KINDR)) go to 1
   call update_trace_EB(DTrace)
   ! this should not really be necessary
   DTrace%BosonicTrace = get_BosonicTrace(DTrace, DStates)
   DTrace%cfgsign = -get_Sign(DTrace, DStates)

end subroutine StepChangeOuter

!===============================================================================
subroutine StepShiftTau_proper()
!===============================================================================
!input
   type(TStates),pointer         :: DStates
   type(TTrace),pointer          :: DTrace
   type(TStates_pointer)         :: pDStates
   type(TTrace_pointer)          :: pDTrace
!local
   type(TLogTr)                  :: TraceNew
   real(KINDR)                   :: DeltaTau, BosTraceNew
   real(KINDR)                   :: outertau_old, outer_state_factor
   real(KINDR)                   :: tauwin_min_old, tauwin_max_old
   real(KINDR),pointer           :: tau_old(:)
   integer                       :: i, iS, NPairs_old, window_position_old
   integer                       :: wrapped(NBands, 2, 2)
   type(TOper),pointer           :: Element, first_old, first_temp, last_old
   type(TOper),pointer           :: prewin_old, postwin_old
   type(TOperPointer)            :: ffirst_old(NBands, 2), ffirst_temp(NBands, 2)
   type(TOperPointer)            :: flast_old(NBands, 2), fprewin_old(NBands, 2)

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   DeltaTau = grnd() * DTrace%beta
   if (DeltaTau == 0_KINDR .or. DeltaTau == DTrace%beta) return
   TryShift = TryShift + 1

   ! store old list ends and make list cyclic
   outertau_old = DTrace%beta - DTrace%last%tau + DTrace%first%tau
   first_old => DTrace%first
   first_temp => DTrace%first
   DTrace%first%prev => DTrace%last
   DTrace%last%next => DTrace%first
   last_old => DTrace%last
   do i = 1, NBands
      do iS = 1, 2
         if (associated(DTrace%ffirst(i, iS)%p)) then
            ffirst_old(i, iS)%p => DTrace%ffirst(i, iS)%p
            ffirst_temp(i, iS)%p => DTrace%ffirst(i, iS)%p
            DTrace%ffirst(i, iS)%p%fprev => DTrace%flast(i, iS)%p
            DTrace%flast(i, iS)%p%fnext => DTrace%ffirst(i, iS)%p
            flast_old(i, iS)%p => DTrace%flast(i, iS)%p
         end if
      enddo
   enddo


   ! shift operators toward 0 by DeltaTau and keep track of first elements
   wrapped = 0
   Element => DTrace%last
   allocate(tau_old(DTrace%NOper))
   do i = DTrace%NOper, 1, -1
      tau_old(i) = Element%tau
      if (Element%tau >= DeltaTau) then
         first_temp => Element
         ffirst_temp(Element%Orbital, Element%Spin)%p => Element
      else
         if (Element%has_hyb) wrapped(Element%Orbital, Element%Spin, Element%CA) =&
            wrapped(Element%Orbital, Element%Spin, Element%CA) + 1
      end if
      Element%tau = modulo(Element%tau - DeltaTau, DTrace%beta)
      Element => Element%prev
   end do


   ! set new list ends and "cut" lists between them
   call save_outer_sst(DTrace)
   if (b_statesampling) then
      outer_state_factor = WeightedOuterStateChoice(DTrace, DStates,&
                              outertau_old,&
                              (DTrace%beta - first_temp%prev%tau + first_temp%tau),&
                              DTrace%outer_sst,&
                              DTrace%outer_state,&
                              (/ DTrace%outer_sst /),&
                              (/ first_temp%preceding_sst /))
   else
      ! reject if there are no states included in potential truncation for target sst
      if (DTrace%sst_to_statesindex(DTrace%outer_sst) /= -1) then
         outer_state_factor = 1_KINDR
      else
         outer_state_factor = 0_KINDR
      end if
      DTrace%outer_sst = first_temp%preceding_sst
   end if

   first_temp%prev%next => null()
   DTrace%last => first_temp%prev
   first_temp%prev => null()
   DTrace%first => first_temp
   do i = 1, NBands
      do iS = 1, 2
         if (associated(ffirst_temp(i, iS)%p)) then
            ffirst_temp(i, iS)%p%fprev%fnext => null()
            DTrace%flast(i, iS)%p => ffirst_temp(i, iS)%p%fprev
            ffirst_temp(i, iS)%p%fprev => null()
            DTrace%ffirst(i, iS)%p => ffirst_temp(i, iS)%p
         end if
      enddo
   enddo


   ! store old window state and shift the sliding window as close to
   ! its old position relative to the operators as possible
   window_position_old = DTrace%window_position
   tauwin_min_old = DTrace%tauwin_min
   tauwin_max_old = DTrace%tauwin_max
   NPairs_old = DTrace%NPairs
   prewin_old => DTrace%prewin
   postwin_old => DTrace%postwin
   DTrace%window_position = modulo(window_position_old&
      - modulo(floor((DeltaTau/DTrace%beta) * dble(DTrace%num_windows + 1)), DTrace%num_windows),&
      DTrace%num_windows)
   DTrace%prewin => null()
   do i = 1, NBands
      do iS = 1, 2
         fprewin_old(i, iS)%p => DTrace%fprewin(i, iS)%p
         DTrace%fprewin(i, iS)%p => null()
      enddo
   enddo
   DTrace%postwin => null()
   DTrace%NPairs = -1
   call set_window(DTrace)


   if (outer_state_factor > 0_KINDR) then
      ! the scalar value can only change due to limited numerical
      ! precision, but the stored vectors need to be updated anyway
      TraceNew = get_trace_EB(DTrace, DStates, global=.true.)
      BosTraceNew = get_BosonicTrace(DTrace, DStates)
   else
      TraceNew = TLogTr(log = -huge(0.0_KINDR), sign = 0.0_KINDR)
      BosTraceNew = 0_KINDR
   end if

   ! the move can only be rejected due to limited numerical precision
   ! of the relevant quantities or the potential outer truncation
   ! (i.e. implementation details), so the determinant is not needed
   ! for the probability (no effect from outer truncation and update
   ! technique ensures ratio +-1)
   if (grnd() < abs(outer_state_factor&
                    *(trval(TraceNew/DTrace%Trace))*(BosTraceNew/DTrace%BosonicTrace))) then
      DTrace%PreWinFCountValid = .false.
      DTrace%Trace = TraceNew
      DTrace%BosonicTrace = BosTraceNew
      call update_trace_EB(DTrace, global=.true.)
      call update_InvFullHybrShiftTau(DTrace, DStates, wrapped)
      DTrace%cfgsign = -get_Sign(DTrace, DStates)
      AccShift = AccShift + 1
   else
      ! make list cyclic
      DTrace%first%prev => DTrace%last
      DTrace%last%next => DTrace%first
      do i = 1, NBands
         do iS = 1, 2
            if (associated(DTrace%ffirst(i, iS)%p)) then
               DTrace%ffirst(i, iS)%p%fprev => DTrace%flast(i, iS)%p
               DTrace%flast(i, iS)%p%fnext => DTrace%ffirst(i, iS)%p
            end if
         end do
      end do

      ! restore old list ends and "cut" lists between them
      DTrace%first => first_old
      DTrace%first%prev => null()
      DTrace%last => last_old
      DTrace%last%next => null()
      do i = 1, NBands
         do iS = 1, 2
            if (associated(ffirst_old(i, iS)%p)) then
               DTrace%ffirst(i, iS)%p => ffirst_old(i, iS)%p
               DTrace%ffirst(i, iS)%p%fprev => null()
               DTrace%flast(i, iS)%p => flast_old(i, iS)%p
               DTrace%flast(i, iS)%p%fnext => null()
            end if
         end do
      end do

      ! restore old taus
      Element => DTrace%last
      do i = DTrace%NOper, 1, -1
         Element%tau = tau_old(i)
         Element => Element%prev
      end do

      ! restore old window state
      DTrace%window_position = window_position_old
      DTrace%tauwin_min = tauwin_min_old
      DTrace%tauwin_max = tauwin_max_old
      DTrace%NPairs = NPairs_old
      DTrace%prewin => prewin_old
      do i = 1, NBands
         do iS = 1, 2
            DTrace%fprewin(i, iS)%p => fprewin_old(i, iS)%p
         enddo
      enddo
      DTrace%postwin => postwin_old
      call restore_outer_sst(DTrace)
   end if

   deallocate(tau_old)
end subroutine StepShiftTau_proper

!===============================================================================
subroutine StepFlavourchange_general()
!===============================================================================
!input
   type(TStates),pointer          :: DStates
   type(TTrace),pointer           :: DTrace
   type(TStates_pointer)          :: pDStates
   type(TTrace_pointer)           :: pDTrace
!local
   integer                       :: iB,iS,begin_sst,end_sst,throwaway(2)
   type(TLogTr)                  :: TraceNew
   real(KINDR)                   :: BosTraceNew
   type(TLogDet)                 :: DetNew, DetRat
   type(TOperPointer)            :: Oper(2), begin, end
   real(kindr),pointer           :: FullHybr_offdiag(:,:)
   real(kindr)                   :: prob
   integer :: oldorb1, oldspin1
   integer :: oldorb2, oldspin2
   type(TSubMatrix),pointer      :: FullHybr(:,:)

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   ! Fortran 90 does not guarantee nullified pointers
   nullify(FullHybr_offdiag)

   TryFlavc=TryFlavc+1

   if (propose_flavourexchange_general(DTrace,oper).eqv..false.)then
      return
   endif


   oldorb1=oper(1)%p%orbital
   oldorb2=oper(2)%p%orbital
   oldspin1=oper(1)%p%spin
   oldspin2=oper(2)%p%spin

   ! remove and later reinsert to keep flavour-specific lists correct
   call remove_Oper(DTrace, Oper)

   !!! dice new flavours

   oper(1)%p%orbital=ceiling(grnd()*dble(dstates%NBands))
   oper(1)%p%spin=ceiling(grnd()*dble(2))
   if(oper(1)%p%orbital.eq.0)oper(1)%p%orbital=1
   if(oper(1)%p%spin.eq.0)oper(1)%p%spin=1

   oper(2)%p%orbital=ceiling(grnd()*dble(dstates%NBands))
   oper(2)%p%spin=ceiling(grnd()*dble(2))
   if(oper(2)%p%orbital.eq.0)oper(2)%p%orbital=1
   if(oper(2)%p%spin.eq.0)oper(2)%p%spin=1

   if (.not. process_OperAdd(DTrace, Oper, throwaway, 2))&
      stop "unreachable: reinsertion (1) failed in StepFlavourchange_general"
   call insert_Oper(DTrace, Oper)

   begin%p => Oper(1)%p%prev
   end%p => Oper(2)%p%next
   call get_current_ssts(DTrace, begin, end, begin_sst, end_sst)
   call save_outer_sst(DTrace)
   if (.not. check_sst_sequence(DTrace, DStates, begin, end, begin_sst, end_sst)) then
      call restore_outer_sst(DTrace)
      call remove_Oper(DTrace, Oper)
      oper(1)%p%orbital=oldorb1
      oper(1)%p%spin=oldspin1
      oper(2)%p%orbital=oldorb2
      oper(2)%p%spin=oldspin2
      if (.not. process_OperAdd(DTrace, Oper, throwaway, 2))&
         stop "unreachable: reinsertion (2) failed in StepFlavourchange_general"
      call insert_Oper(DTrace, Oper)
      return
   endif

   Oper(1)%p%calc_new = .true.
   Oper(2)%p%calc_new = .true.

   TraceNew = get_trace_EB(DTrace,DStates)

   if(b_offdiag)then

      FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
      call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag, DetNew)
      DetRat = DetNew/DTrace%det

   else

      FullHybr=>get_InvFullHybr(DTrace,DStates,DetNew)
      DetRat=DetNew*get_FullLogDet(DTrace,DStates)

   endif


   BosTraceNew = get_BosonicTrace(DTrace,DStates)

   prob=abs(wrat(TraceNew/DTrace%Trace, DetRat)*BosTraceNew/DTrace%BosonicTrace)

   if(grnd().lt.prob)then
      DTrace%Det = DetNew
      DTrace%Trace = TraceNew
      DTrace%BosonicTrace = BosTraceNew
      DTrace%cfgsign = -get_Sign(DTrace, DStates)

      !update hybridizaton matrix
      if(b_offdiag)then

         if(associated(DTrace%MInv_full%Mat))then
            deallocate(DTrace%MInv_full%Mat)
         endif
         DTrace%MInv_full%Mat=>FullHybr_offdiag

      else

         DTrace%NPairs = -1
         do iB=1,NBands
         do iS=1,2
            deallocate(DTrace%MInv(iB,iS)%Mat)
         enddo
         enddo
         deallocate(DTrace%MInv)
         DTrace%MInv=>FullHybr

      endif

      AccFlavc=AccFlavc+1

      call update_trace_EB(DTrace)

   else

      !!! revert changes and return !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call restore_outer_sst(DTrace)
      call remove_Oper(DTrace, Oper)
      oper(1)%p%orbital=oldorb1
      oper(1)%p%spin=oldspin1
      oper(2)%p%orbital=oldorb2
      oper(2)%p%spin=oldspin2
      if (.not. process_OperAdd(DTrace, Oper, throwaway, 2))&
         stop "unreachable: reinsertion (3) failed in StepFlavourchange_general"
      call insert_Oper(DTrace, Oper)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(b_offdiag)then

         deallocate(FullHybr_offdiag)

      else

         do iB=1,NBands
         do iS=1,2
            deallocate(FullHybr(iB,iS)%Mat)
         enddo
         enddo
         deallocate(FullHybr)

      endif
   endif

   Oper(1)%p%calc_new = .false.
   Oper(2)%p%calc_new = .false.
end subroutine StepFlavourchange_general

!===============================================================================
subroutine StepWormAdd(Sector,flavor)
!===============================================================================
!input
   integer, intent(inout)     :: Sector
   integer, optional          :: flavor
   integer                    :: tflavor
   
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
!local
   integer                       :: begin_sst,end_sst
   type(TLogTr)                  :: TraceNew
   type(TOperPointer)            :: begin,end
   integer                       :: iO,i,FPos(2)
   real(KINDR)                   :: BosTraceNew
   real(KINDR)                   :: rand,prefact
   type(TOperPointer)            :: Oper(NOperWorm(Sector))
 

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   !generate random general flavor index for band spin pattern
   if(.not.present(flavor).or.flavor.eq.0) then 
      tflavor = ceiling(grnd()*((2*NBands)**NOperWorm(Sector)))
   else
      tflavor = flavor
   endif

   if(tflavor .eq. 0) tflavor = 1
   

   !calculate band-spin compound index for umatrix in improved estimator
   !call index2component_general(Nbands, NOperWorm(Sector), tflavor, bs, b, s)
   
   TryWormAdd(Sector)=TryWormAdd(Sector)+1
   
   !allocating oper stuff, in case memory already allocated
   !we use the opers stored in OperPool
   if(DTrace%iOperPool<NOperWorm(Sector))then
      do iO=1,size(Oper)
         allocate(Oper(iO)%p)
         Oper(iO)%p%calc_new=.true.
      enddo
   else
      do iO=1,size(Oper)
         Oper(iO)%p=>DTrace%OperPool(DTrace%iOperPool)%p
         DTrace%OperPool(DTrace%iOperPool)%p=>null()
         DTrace%iOperPool=DTrace%iOperPool-1
         Oper(iO)%p%calc_new=.true.
         Oper(iO)%p%cache_written = .false.
      enddo
   endif

   !set band spin pattern and times for worm according to flavor and Sector
   !for improved estimators the 'inner' operators are chosen randomly
   !for a given tflavor
   if(Sector==3 .or. Sector==5) then
      call set_IEOperFlavor(Oper,NOperWorm(Sector),NBands,flavor,u_o,u_s)
      u_g=2*(u_o-1)+u_s
   else
     call set_OperFlavor(Oper,NOperWorm(Sector),NBands,tflavor)
   endif
   call set_OperTime(Oper,NOperWorm(Sector),DTrace%beta,DTrace%equal_time_offset,Sector)
   
   !attaching tags to the generated operators
   !make sure they have no hyb lines attached
   allocate(DTrace%wormContainer(NOperWorm(Sector)))
   do i=1,NOperWorm(Sector)
      Oper(i)%p%has_hyb=.false.
      DTrace%wormContainer(i)%p=>Oper(i)%p
   enddo
   
   !sort the generated operators to the correct position in the trace
   call process_OperAdd_global(DTrace,Oper,FPos,NOperWorm(Sector))

   ! determine substring of trace that needs quantum number checking
   begin%p => Oper(1)%p%prev
   end%p => Oper(size(Oper))%p%next
   call get_current_ssts(DTrace, begin, end, begin_sst, end_sst)

   !we now add the oper pairs
   do i=1,NOperWorm(Sector),2
      call insert_Oper(DTrace,Oper(i:i+1))
   enddo

   call get_Proposal(Sector,prefact)
   
   call save_outer_sst(DTrace)
   !now we can check the quantum number violations
   if((.not.check_EqualTime(DTrace,Sector)).or.(prefact.eq.0d0)&
       .or. (.not. check_sst_sequence(DTrace, DStates, begin, end, begin_sst, end_sst))) then
      call restore_outer_sst(DTrace)
      !if we violate quantum numbers we need to remove operators again before exiting
      do i=1,NOperWorm(Sector),2
         call remove_Oper(DTrace,Oper(i:i+1))
      enddo
      
      deallocate(DTrace%wormContainer)
      do iO=1,size(Oper)
         if(associated(Oper(iO)%p%state))deallocate(Oper(iO)%p%state)
         if(associated(Oper(iO)%p%cache))deallocate(Oper(iO)%p%cache)
         if(associated(Oper(iO)%p%slogmax))deallocate(Oper(iO)%p%slogmax)
         if(associated(Oper(iO)%p%clogmax))deallocate(Oper(iO)%p%clogmax)
         Oper(iO)%p%calc_new=.true.
         Oper(iO)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+iO)%p=>Oper(iO)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+NOperWorm(Sector)
      
      AccQNWormAdd(Sector)=AccQNWormAdd(Sector)+1
      Sector=1
      return
   endif
    
   TraceNew = get_trace_EB(DTrace,DStates,global=.true.)
   BosTraceNew = get_BosonicTrace(DTrace,DStates)
   
   !Metropolis                        
   rand=grnd()
   if(rand.lt.dabs(prefact&
      *trval(TraceNew/DTrace%Trace)*BosTraceNew/DTrace%BosonicTrace))then
     
      DTrace%Trace=TraceNew
      DTrace%BosonicTrace=BosTraceNew
      
      DTrace%cfgsign = -get_Sign(DTrace, DStates)
     
      call update_trace_EB(DTrace,global=.true.)
      
      AccWormAdd(Sector)=AccWormAdd(Sector)+1
      
      !update control for worm measurement in order to save on Fourier transform
      if(Sector>1) then
         isNew(Sector)=.true.
      endif

   else
      call restore_outer_sst(DTrace)
      do i=1,NOperWorm(Sector),2
         call remove_Oper(DTrace,Oper(i:i+1))
      enddo
            
      deallocate(DTrace%wormContainer)
      
      !put operators memory in operator pool
      do iO=1,size(Oper)
         if(associated(Oper(iO)%p%state))deallocate(Oper(iO)%p%state)
         if(associated(Oper(iO)%p%cache))deallocate(Oper(iO)%p%cache)
         if(associated(Oper(iO)%p%slogmax))deallocate(Oper(iO)%p%slogmax)
         if(associated(Oper(iO)%p%clogmax))deallocate(Oper(iO)%p%clogmax)
         Oper(iO)%p%calc_new=.true.
         Oper(iO)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+iO)%p=>Oper(iO)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+NOperWorm(Sector)
    
      !make sure we change sector after updating pool
      Sector=1
   endif

end subroutine StepWormAdd

!===============================================================================
!> Removes worm added by StepWormAdd
subroutine StepWormRem(Sector)
!===============================================================================
!input
   integer, intent(inout)        :: Sector
!local
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   
   integer                       :: begin_sst,end_sst
   type(TLogTr)                  :: TraceNew
   real(KINDR)                   :: throwaway
   type(TOperPointer)            :: begin,end
   integer                       :: iO,i,FPos(2)
   real(KINDR)                   :: rand,prefact
   real(KINDR)                   :: BosTraceNew
   type(TOperPointer)            :: Oper(NOperWorm(Sector))
   type(TOper),pointer           :: Element
      
   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   nullify(Element)

   TryWormRem(Sector)=TryWormRem(Sector)+1
  
   !TODO: might be replaced by simply taking the adress from the worm container
   if(.not.gen_OperRemove(DTrace,Oper,FPos,NOperWorm(Sector),throwaway,with_hyb=.false.))then
      AccQNWormRem(Sector)=AccQNWormRem(Sector)+1
      return
   endif

   ! determine substring of trace that needs quantum number checking
   begin%p => Oper(1)%p%prev
   end%p => Oper(size(Oper))%p%next
   call get_current_ssts(DTrace, begin, end, begin_sst, end_sst)
   
   !we now remove the oper pairs
   do i=1,NOperWorm(Sector),2
      call remove_Oper(DTrace,Oper(i:i+1))
   enddo

   call save_outer_sst(DTrace)
   !now we can check the quantum number violations
   if (.not. check_sst_sequence(DTrace, DStates, begin, end, begin_sst, end_sst)) then
      call restore_outer_sst(DTrace)
     !if we violate quantum numbers we need to insert operators again before exiting
     do i=1,NOperWorm(Sector),2
        call insert_Oper(DTrace,Oper(i:i+1))
     enddo
     AccQNWormRem(Sector)=AccQNWormRem(Sector)+1
     return
   endif
   
   ! we need to set the operators next to the newly removed
   ! operators to calc_new=true
   do iO=1,size(Oper)
      if(associated(Oper(iO)%p%next))then
         Oper(iO)%p%next%calc_new=.true.
      endif
   enddo
   do iO=size(Oper),1,-1
      if(associated(Oper(iO)%p%prev))then
        Oper(iO)%p%prev%calc_new=.true.
      endif
   enddo

   TraceNew = get_trace_EB(DTrace,DStates,global=.true.)
   BosTraceNew = get_BosonicTrace(DTrace,DStates)
   
   call get_Proposal(Sector,prefact)
                     
   rand=grnd()
   if(rand.lt.dabs(1d0/prefact&
      *trval(TraceNew/DTrace%Trace)*BosTraceNew/DTrace%BosonicTrace))then
      
      DTrace%Trace=TraceNew
      DTrace%BosonicTrace=BosTraceNew
      
      !deallocating worm container before sign measure
      deallocate(DTrace%wormContainer)
      
      DTrace%cfgsign = -get_Sign(DTrace, DStates)
      
      call update_trace_EB(DTrace,global=.true.)

      do iO=1,size(Oper)
         if(associated(Oper(iO)%p%state))deallocate(Oper(iO)%p%state)
         if(associated(Oper(iO)%p%cache))deallocate(Oper(iO)%p%cache)
         if(associated(Oper(iO)%p%slogmax))deallocate(Oper(iO)%p%slogmax)
         if(associated(Oper(iO)%p%clogmax))deallocate(Oper(iO)%p%clogmax)
         DTrace%OperPool(DTrace%iOperPool+iO)%p=>Oper(iO)%p
         Oper(iO)%p%calc_new=.true.
         Oper(iO)%p%cache_written = .false.
      enddo
      DTrace%iOperPool=DTrace%iOperPool+NOperWorm(Sector)
      
      !if operator rem accepted and without hyb lines, change sector
      AccWormRem(Sector)=AccWormRem(Sector)+1
      Sector=1      


   else
      call restore_outer_sst(DTrace)
      do i=1,NOperWorm(Sector),2
        call insert_Oper(DTrace,Oper(i:i+1))
      enddo

      ! TODO: don't loop over the entire list
      Element=>DTrace%first
      do while(associated(Element))
        Element%calc_new=.false.
        Element=>Element%next
      enddo
      
   endif

end subroutine StepWormRem

!for offdiagonal hybridization functions one needs a mixed worm-hybridization
!4-operator move. This may be bypassed by inserting a diagonal worm and four
!(offdiagonal) hybridization operators followed by offdiagonal worm replacement 
!step and removing a pair of hybridization operators 
!in terms of compoment sampling a mixed worm insert is necessary for the
!single particle Green's function
!===============================================================================
subroutine StepWormHybAdd(Sector,flavor)
!===============================================================================
!input
   integer, intent(inout)     :: Sector
   integer, optional          :: flavor
   integer                    :: wflavor,hybflavor
   
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
!local
   integer                       :: begin_sst,end_sst
   type(TLogTr)                  :: TraceNew
   type(TOperPointer)            :: begin,end
   integer                       :: iO,i,FPos(4)
   real(KINDR)                   :: BosTraceNew
   real(KINDR)                   :: rand,prefact
   type(TOperPointer)            :: Oper(4)
   real(KINDR)                   :: S
   type(TLogDet)                 :: DetNew, DetRat
   real(KINDR),pointer           :: Q(:),R(:)
   real(KINDR),pointer           :: temp(:,:)
   real(kindr),pointer           :: FullHybr_offdiag(:,:)

   nullify(FullHybr_offdiag)

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   if(Sector .ne. 2) then
      stop 'StepWormHybAdd only implemented for one-particle Greens function'
   endif
   !furthermore Sector 2 -> NOperWorm(2)=2

   !generate worm flavor
   if(.not.present(flavor).or.flavor.eq.0) then
      wflavor = ceiling(grnd()*((2*NBands)**2))
   else
      wflavor = flavor
   endif
   !hybflavor determined randomly
   hybflavor = ceiling(grnd()*((2*NBands)**2))
   
   if(wflavor .eq. 0) wflavor = 1
   if(hybflavor .eq. 0) hybflavor = 1
   
   !TODO:change to account for WormHyb Steps
   TryWormAdd(Sector)=TryWormAdd(Sector)+1
   
   !allocating oper stuff, in case memory already allocated
   !we use the opers stored in OperPool
   if(DTrace%iOperPool<4)then
      do iO=1,size(Oper)
         allocate(Oper(iO)%p)
         Oper(iO)%p%calc_new=.true.
      enddo
   else
      do iO=1,size(Oper)
         Oper(iO)%p=>DTrace%OperPool(DTrace%iOperPool)%p
         DTrace%OperPool(DTrace%iOperPool)%p=>null()
         DTrace%iOperPool=DTrace%iOperPool-1
         Oper(iO)%p%calc_new=.true.
         Oper(iO)%p%cache_written = .false.
      enddo
   endif

   !Oper(1:2) worm operators, Oper(3:4) hyb operators
   call set_OperFlavor(Oper(1:2),2,NBands,wflavor)
   call set_OperFlavor(Oper(3:4),2,NBands,hybflavor)
   call set_OperTime(Oper,4,DTrace%beta,DTrace%equal_time_offset,Sector)
   
   !attaching tags to the generated operators
   allocate(DTrace%wormContainer(2))
   do i=1,2
      Oper(i)%p%has_hyb=.false.
      Oper(i+2)%p%has_hyb=.true.
      DTrace%wormContainer(i)%p=>Oper(i)%p
   enddo
   
   !sort the generated operators to the correct position in the trace
   call process_OperAdd_global(DTrace,Oper(1:2),FPos(1:2),2)
   call insert_Oper(DTrace,Oper(1:2))
   call process_OperAdd_global(DTrace,Oper(3:4),FPos(3:4),2)
   call insert_Oper(DTrace,Oper(3:4))

   ! determine substring of trace that needs quantum number checking
   if (Oper(1)%p%tau < Oper(3)%p%tau) then
      begin%p => Oper(1)%p%prev
   else
      begin%p => Oper(3)%p%prev
   end if
   if (Oper(2)%p%tau > Oper(4)%p%tau) then
      end%p => Oper(2)%p%next
   else
      end%p => Oper(4)%p%next
   end if
   call get_current_ssts(DTrace, begin, end, begin_sst, end_sst)

   call save_outer_sst(DTrace)
   !now we can check the quantum number violations
   if((.not.check_EqualTime(DTrace,Sector))&
       .or. (.not. check_sst_sequence(DTrace, DStates, begin, end, begin_sst, end_sst))) then
      call restore_outer_sst(DTrace)
      !if we violate quantum numbers we need to remove operators again before exiting
      call remove_Oper(DTrace,Oper(1:2))
      call remove_Oper(DTrace,Oper(3:4))
      deallocate(DTrace%wormContainer)
      do iO=1,size(Oper)
         if(associated(Oper(iO)%p%state))deallocate(Oper(iO)%p%state)
         if(associated(Oper(iO)%p%cache))deallocate(Oper(iO)%p%cache)
         if(associated(Oper(iO)%p%slogmax))deallocate(Oper(iO)%p%slogmax)
         if(associated(Oper(iO)%p%clogmax))deallocate(Oper(iO)%p%clogmax)
         Oper(iO)%p%calc_new=.true.
         Oper(iO)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+iO)%p=>Oper(iO)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+4
      
      !TODO:change to account for WormHyb Steps
      AccQNWormAdd(Sector)=AccQNWormAdd(Sector)+1
      Sector=1
      return
   endif
    
   if(b_offdiag)then
   
      if(b_full_offdiag)then
         FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
         call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag,DetNew)
         DetRat = DetNew/DTrace%Det
      else
         !!! fast updates
         call add_ftau_offset(DTrace, DStates, Oper(3:4), FPos(3:4))
         DetRat = get_LogDetRatPairAdd_full(dtrace,Oper(3:4),DStates,Q,R,S,FPos(3:4))
         DetNew = DetRat * DTrace%Det
      endif

   else

      DetRat = get_LogDetRatPairAdd(DTrace,DTrace%MInv,Oper(3:4),Q,R,S,FPos(3:4))
      DetNew = DetRat * DTrace%Det

   endif

   TraceNew = get_trace_EB(DTrace,DStates,global=.true.)
   BosTraceNew = get_BosonicTrace(DTrace,DStates)

   !Metropolis                        
   prefact = 2.0*((DTrace%beta*NBands)**4)*wormEta(Sector)/((DTrace%Noper/2-1)**2)
   rand=grnd()
   if(rand.lt.dabs(prefact&
      *wrat(TraceNew/DTrace%Trace,DetRat)*BosTraceNew/DTrace%BosonicTrace))then
     
      if(.not.b_offdiag)then

         temp=>get_MatAdd(&
         DTrace%MInv(Oper(3)%p%Orbital,Oper(3)%p%Spin)%Mat,&
         DTrace%NOSOper(Oper(3)%p%Orbital,Oper(3)%p%Spin)/2-1,Q,R,S,FPos(3:4))
         if(associated(DTrace%MInv(Oper(3)%p%Orbital,Oper(3)%p%Spin)%Mat))then
            deallocate(DTrace%MInv(Oper(3)%p%Orbital,Oper(3)%p%Spin)%Mat)
         endif
         DTrace%MInv(Oper(3)%p%Orbital,Oper(3)%p%Spin)%Mat=>temp

      else


         if(b_full_offdiag)then
            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>FullHybr_offdiag
         else
            !!! fast matrix generation
            !!! (Noper-2-2)/2 (-worm-hyb)
            temp=>get_MatAdd(&
            DTrace%MInv_full%Mat,&
            DTrace%Noper/2-2,Q,R,S,FPos(3:4))

            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>temp
         endif

      endif

      DTrace%Det=DetNew
      DTrace%Trace=TraceNew
      DTrace%BosonicTrace=BosTraceNew
      
      DTrace%cfgsign = -get_Sign(DTrace, DStates)
     
      call update_trace_EB(DTrace, global=.true.)
      
      !TODO:change to account for WormHyb Steps
      AccWormAdd(Sector)=AccWormAdd(Sector)+1
      
      !update control for worm measurement in order to save on Fourier transform
      if(Sector>1) then
         isNew(Sector)=.true.
      endif

      ! TODO: implement update to do this faster
      DTrace%prewinfcountvalid = .false.
      DTrace%NPairs = -1
   else
      call restore_outer_sst(DTrace)
      call remove_Oper(DTrace,Oper(1:2))
      call remove_Oper(DTrace,Oper(3:4))

      deallocate(DTrace%wormContainer)
      
      if(associated(FullHybr_offdiag))then
         deallocate(FullHybr_offdiag)
      endif

      !put operators memory in operator pool
      do iO=1,size(Oper)
         if(associated(Oper(iO)%p%state))deallocate(Oper(iO)%p%state)
         if(associated(Oper(iO)%p%cache))deallocate(Oper(iO)%p%cache)
         if(associated(Oper(iO)%p%slogmax))deallocate(Oper(iO)%p%slogmax)
         if(associated(Oper(iO)%p%clogmax))deallocate(Oper(iO)%p%clogmax)
         Oper(iO)%p%calc_new=.true.
         Oper(iO)%p%cache_written = .false.
         DTrace%OperPool(DTrace%iOperPool+iO)%p=>Oper(iO)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+4
    
      !make sure we change sector after updating pool
      Sector=1
   endif

end subroutine StepWormHybAdd

!===============================================================================
!> Removes worm/hyb 4-pair added by StepWormHybAdd
subroutine StepWormHybRem(Sector)
!===============================================================================
!input
   integer, intent(inout)        :: Sector
!local
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   
   integer                       :: begin_sst,end_sst
   type(TLogTr)                  :: TraceNew
   type(TOperPointer)            :: begin,end
   integer                       :: iO,FPos(4)
   real(KINDR)                   :: rand,prefact
   real(KINDR)                   :: BosTraceNew
   type(TOperPointer)            :: Oper(4)
   type(TOper),pointer           :: Element
   real(KINDR)                   :: S
   type(TLogDet)                 :: DetNew, DetRat
   real(KINDR),pointer           :: temp(:,:)
   real(kindr),pointer            :: FullHybr_offdiag(:,:)
      

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   if(Sector .ne. 2) then
      stop 'StepWormHybRem only implemented for one-particle Greens function'
   endif
   !furthermore Sector 2 -> NOperWorm(2)=2

   !TODO:change to account for WormHyb Steps
   TryWormRem(Sector)=TryWormRem(Sector)+1
  
   !select 2 worm operators and two other arbitrary operators and save to FPos
   if(.not.gen_OperRemove_global(DTrace,DStates,Oper(1:2),FPos(1:2),2,.false.)) then 
      AccQNWormRem(Sector)=AccQNWormRem(Sector)+1
      return
   endif
   call remove_Oper(DTrace,Oper(1:2))
   ! deallocate worm container to get correct operator count in gen_OperRemove_global
   deallocate(DTrace%wormContainer)
   
   if(.not.gen_OperRemove_global(DTrace,DStates,Oper(3:4),FPos(3:4),2,.true.)) then
      AccQNWormRem(Sector)=AccQNWormRem(Sector)+1
      call insert_Oper(DTrace, Oper(1:2))
      call restore_wormContainer()
      return
   endif
   call remove_Oper(DTrace,Oper(3:4))

   ! determine substring of trace that needs quantum number checking
   if (Oper(1)%p%tau < Oper(3)%p%tau) then
      begin%p => Oper(1)%p%prev
   else
      begin%p => Oper(3)%p%prev
   end if
   if (Oper(2)%p%tau > Oper(4)%p%tau) then
      end%p => Oper(2)%p%next
   else
      end%p => Oper(4)%p%next
   end if
   call get_current_ssts(DTrace, begin, end, begin_sst, end_sst)
   
   call save_outer_sst(DTrace)
   !now we can check the quantum number violations
   if (.not. check_sst_sequence(DTrace, DStates, begin, end, begin_sst, end_sst)) then
      call restore_outer_sst(DTrace)
     !if we violate quantum numbers we need to insert operators again before exiting
     call insert_Oper(DTrace,Oper(3:4))
     call insert_Oper(DTrace,Oper(1:2))
     call restore_wormContainer()
     AccQNWormRem(Sector)=AccQNWormRem(Sector)+1
     return
   endif
   
   ! we need to set the operators next to the newly removed
   ! operators to calc_new=true
   do iO=1,size(Oper)
      if(associated(Oper(iO)%p%next))then
         Oper(iO)%p%next%calc_new=.true.
      endif
   enddo
   do iO=size(Oper),1,-1
      if(associated(Oper(iO)%p%prev))then
        Oper(iO)%p%prev%calc_new=.true.
      endif
   enddo

   TraceNew = get_trace_EB(DTrace,DStates,global=.true.)
   BosTraceNew = get_BosonicTrace(DTrace,DStates)

   if(.not.b_offdiag)then

      DetRat = get_LogDetRatPairRem(DTrace,Oper(3:4),FPos(3:4),S)
      DetNew = DetRat * DTrace%Det

   else

      if(b_full_offdiag)then
         FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
         call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag,DetNew)
         DetRat = DetNew/DTrace%Det
      else
         !!! fast updates
         call add_ftau_offset(dtrace,dstates,Oper(3:4),Fpos(3:4))
         DetRat = get_LogDetRatPairRem_full(DTrace,Oper(3:4),FPos(3:4),S)
         DetNew = DetRat * DTrace%Det
      endif

   endif

   
   prefact = 2.0*((DTrace%beta*NBands)**4)*wormEta(Sector)/((DTrace%Noper/2+1)**2)
   rand=grnd()
   if(rand.lt.dabs(1d0/prefact&
      *wrat(TraceNew/DTrace%Trace,DetRat)*BosTraceNew/DTrace%BosonicTrace))then
      
      if(.not.b_offdiag)then

         temp=>get_MatRem(&
         DTrace%MInv(Oper(3)%p%Orbital,Oper(3)%p%Spin)%Mat,&
               DTrace%NOSOper(Oper(3)%p%Orbital,Oper(3)%p%Spin)/2+1,FPos(3:4))
         deallocate(DTrace%MInv(Oper(3)%p%Orbital,Oper(3)%p%Spin)%Mat)
         DTrace%MInv(Oper(3)%p%Orbital,Oper(3)%p%Spin)%Mat=>temp
      else

         if(b_full_offdiag)then
            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>FullHybr_offdiag
         else
            !!! fast matrix generation
            !!! (Noper-2+2)/2 (-worm+hyb)
            temp=>get_MatRem(&
                  DTrace%MInv_full%Mat,&
                  DTrace%NOper/2+1,FPos(3:4))
         
            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>temp
         endif

      endif

      !we save the determinant of the current configuration
      DTrace%Det=DetNew
      DTrace%Trace=TraceNew
      DTrace%BosonicTrace=BosTraceNew
      
      DTrace%cfgsign = -get_Sign(DTrace, DStates)
      
      call update_trace_EB(DTrace, global=.true.)

      do iO=1,size(Oper)
         if(associated(Oper(iO)%p%state))deallocate(Oper(iO)%p%state)
         if(associated(Oper(iO)%p%cache))deallocate(Oper(iO)%p%cache)
         if(associated(Oper(iO)%p%slogmax))deallocate(Oper(iO)%p%slogmax)
         if(associated(Oper(iO)%p%clogmax))deallocate(Oper(iO)%p%clogmax)
         DTrace%OperPool(DTrace%iOperPool+iO)%p=>Oper(iO)%p
         Oper(iO)%p%calc_new=.true.
         Oper(iO)%p%cache_written = .false.
      enddo
      DTrace%iOperPool=DTrace%iOperPool+4
      
      !TODO: might be replaced by simply taking the adress from the worm container
      AccWormRem(Sector)=AccWormRem(Sector)+1
      Sector=1      

      ! TODO: implement update to do this faster
      DTrace%prewinfcountvalid = .false.
      DTrace%NPairs = -1

   else
      call restore_outer_sst(DTrace)
      call insert_Oper(DTrace,Oper(3:4))
      call insert_Oper(DTrace,Oper(1:2))

      call restore_wormContainer()

      if (b_offdiag .and. b_full_offdiag) deallocate(FullHybr_offdiag)

      ! TODO: don't loop over the entire list
      Element=>DTrace%first
      do while(associated(Element))
        Element%calc_new=.false.
        Element=>Element%next
      enddo
      
   endif

contains
   subroutine restore_wormContainer()
      allocate(DTrace%wormContainer(2))
      if(Oper(1)%p%CA.eq.2) then
         DTrace%wormContainer(1)%p=>Oper(1)%p
         DTrace%wormContainer(2)%p=>Oper(2)%p
      else
         DTrace%wormContainer(1)%p=>Oper(2)%p
         DTrace%wormContainer(2)%p=>Oper(1)%p
      endif
   end subroutine restore_wormContainer
end subroutine StepWormHybRem

!===============================================================================
!> Calculates the proposal ratio for entering different Sectors
!TODO: all flavor-independent quantities can probably be included in wormEta
!this is not the case for the U matrix of the improved estimators
subroutine get_Proposal(Sector,prefact)
!===============================================================================
!input
integer,intent(in)                ::  Sector
real(KINDR),intent(out)           ::  prefact

type(TTrace),pointer              :: DTrace
type(TStates_pointer)             :: pDStates
type(TTrace_pointer)              :: pDTrace

!local
integer                           ::  i
integer                           ::  bs(NOperWorm(Sector))


pDStates=transfer(ipDStates,pDStates)
pDTrace=transfer(ipDTrace,pDTrace)
DTrace => pDTrace%ptr
   
   !2 operators with 2 random times and one random flavor
   if(Sector .eq. 1) then

      prefact = (DTrace%beta**NOperWorm(Sector))*NBands

   elseif(Sector .eq. 2) then

      prefact = ((DTrace%beta*NBands)**NOperWorm(Sector))*wormEta(Sector)

   !5 operatos with 5 random times and 5 random flavors
   elseif(Sector .eq. 4) then 

      prefact = ((DTrace%beta*Nbands)**NOperWorm(Sector))*wormEta(Sector)

   !4/6 operators with 2/4 random times and 4/6 flavors
   elseif(Sector .eq. 3 .or. Sector .eq. 5) then
      
      do i=1,NOperWorm(Sector)
         bs(i) = 2*(DTrace%wormContainer(i)%p%orbital-1) + DTrace%wormContainer(i)%p%spin
      enddo
      
      !the umatrix prefactor goes into the sampling
      prefact = 0.5d0*(u_matrix(u_g,bs(2),bs(1),bs(3))-u_matrix(bs(2),u_g,bs(1),bs(3)))
      
      prefact = prefact*(DTrace%beta**(NOperWorm(Sector)-2))*(Nbands**(NOperWorm(Sector)+1))*wormEta(Sector)

   !4 operators with 2 random times and 4 flavors
   elseif(Sector .eq. 6 .or. Sector .eq. 7) then 
 
      prefact = (DTrace%beta**2)*(Nbands**4)*wormEta(Sector)

   !4 operators with 3 random times and 4 flavors
   elseif(Sector .eq. 8 .or. Sector .eq. 9) then 

      prefact = (DTrace%beta**3)*(Nbands**4)*wormEta(Sector)
   
   endif

end subroutine get_Proposal

!===============================================================================
!> Switches position of a worm operator with a hybridization operator 
!! of the same flavor by calling propose_WormReplace and do_WormReplace.

!! we use this move to produce worms of lengths close to beta/2, as we tend to
!! block them out with hybridization operators due to quantum number violation
!! for high beta

!! component sampling -> only allow replacements of same flavor
subroutine StepWormReplace(Sector)
!===============================================================================
!input
   integer,intent(in)         :: Sector
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   integer                    :: operPos,tmpPos, iN
   integer                    :: wb,ws,wca
   integer                    :: N,rand
   real(KINDR)                :: DetRat,sgn
   type(TOper),pointer        :: OperHyb
   type(TOper),pointer        :: Element
   real(KINDR),pointer        :: u(:),vdag(:)
   logical                    :: force_diagonal=.true.   

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   nullify(Element);nullify(vdag);nullify(u)
   
   if((.not.allocated(DTrace%wormContainer)).or.(Sector.eq.1)) then 
      stop "Cannot replace a hyb operator with a worm operator"
   endif
   
   TryWormRep(Sector)=TryWormRep(Sector)+1
   
   !select one of the two worm operators for exchange
   if(Sector .eq. 2 .or. Sector .eq. 4) then
      rand=ceiling(grnd()*size(DTrace%wormContainer))
      if(rand.eq.0) rand=1     
   elseif(Sector .eq. 3 .or. Sector .eq. 5) then
      !we only exchange the operators which are not equal time operators
      rand=ceiling(grnd()*size(DTrace%wormContainer(4:)))+3
      if(rand.eq.3) rand=4
   elseif(Sector .eq. 8) then
     !we only exchange the operators which are not equal time operators
      rand=ceiling(grnd()*2)
      if(rand.eq.0) rand=1
   elseif(Sector .eq. 9) then
      !we only exchange the operators which are not equal time operators
      rand=ceiling(grnd()*2)
      if(rand.eq.0) rand=1
      if(rand.eq.2) rand=3
   !we do not attempt replacment moves for two legged GF
   else
      return
   endif
   
   wb=DTrace%wormContainer(rand)%p%orbital
   ws=DTrace%wormContainer(rand)%p%spin
   wca=DTrace%wormContainer(rand)%p%CA
   
   if(b_offdiag) then
      N=(DTrace%NOper-size(DTrace%wormContainer))/2
   else
      N=DTrace%NOSOper(wb,ws)/2
   endif

   !if we cant find operators for exchange we immediately reject move
   if(N.eq.0) return
   
   !build vectors for the specific flavor and type
   allocate(u(N))
   allocate(vdag(N)) 
   
   !pick random operator of same type as worm operator
   operPos=int(grnd()*N)+1
   
   !fix rare cases were grand() gives exactly 1
   if(operPos.eq.N+1) then
      operPos=N
   endif 

   !find operator in list
   tmpPos=operPos
   Element=>DTrace%first
 
   do iN=1,DTrace%NOper
      if(Element%has_hyb) then
         if(b_offdiag) then 
            if(Element%CA.eq.wca) tmpPos=tmpPos-1
         else
            if(Element%CA.eq.wca.and.&
               Element%Orbital.eq.wb.and.&
               Element%Spin.eq.ws) tmpPos=tmpPos-1
         endif
         if(tmpPos.eq.0) exit
      endif
      Element=>Element%next
   enddo

   !assigning hyb operator found
   OperHyb=>Element
   
   !component sampling requires diagonal exchanges only
   if(b_offdiag.and.force_diagonal) then
     if(OperHyb%Orbital.ne.wb .or. OperHyb%Spin.ne.ws) goto 99
   endif

   call propose_WormReplace(DTrace,DStates,DTrace%wormContainer(rand)%p,OperHyb,u,vdag,detRat,sgn)

   !metropolis detRat only needs to be calculated for flavor block, not entire matrix
   if(grnd().lt.abs(detRat))then
   
      !updating matrix,determinant etc.
      call do_WormReplace(DTrace,DStates,DTrace%wormContainer(rand)%p,OperHyb,u,vdag,detRat)
      
      !calculating the sign (should be the same as sgn from propose_wormReplace)
      DTrace%cfgsign=-get_Sign(DTrace,DStates)      

      !udpade acceptance
      AccWormRep(Sector)=AccWormRep(Sector)+1
      
      !update control for worm measurement in oder to save on Fourier transform
      isNew(Sector)=.true.

   endif
   
   !clean up
99 deallocate(u,vdag)
   
end subroutine StepWormReplace

!===============================================================================
!> this subroutine tries to dermine the best wormEta for each sector, such that each
!  sector is equally visited by finding the correct reweighing factor
subroutine findEta(iSector,iComponent,tol,Nfix,maxIter)
   use type_progress

   integer                    :: iSector,iComponent
   real(KINDR)                :: tol
   integer(C_INT64_T)         :: Nfix
   integer                    :: maxIter
!local
   real(KINDR)                :: rand_num
   integer(c_int64_t)         :: i, j
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   integer                    :: Sector
   type(progress)             :: p

   real(KINDR),allocatable    :: taus(:)
   integer,allocatable        :: orbs(:), spins(:), cas(:),hashybs(:)
   integer                    :: outer_sst, outer_state

   logical                    :: worm_offdiag
   integer                    :: bs(2), b(2), s(2)

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr
   

   !save starting configuration
   allocate(taus(DTrace%NOper),orbs(DTrace%NOper),&
            spins(DTrace%NOper),cas(DTrace%NOper),hashybs(DTrace%NOper))
   call get_mc_config(DTrace%NOper, taus, orbs, spins, cas, hashybs, outer_sst, outer_state)
 
   !starting in partition function space     
   Sector = 1

   p%adaptrate = 0.1
   call pstart(p, int(maxIter,PINT), &
               title='EtaSearch' // simstr // ':', fancy=fancy_prog)

   !determine if off_diagonal one-particle worm
   worm_offdiag=.false.
   if(iSector .eq. 2) then
      call index2component_general(Nbands, 2, iComponent, bs, b, s)
      if(bs(1).ne.bs(2)) worm_offdiag=.true.
   endif

   do i=1,maxIter
      do j=1,Nfix
         if(Sector==1) then
            CntSampling(1) = CntSampling(1) + 1
            rand_num=grnd()
            !insertion of hyb pair
            if(rand_num<(1d0-PercentageTauShiftMove-PercentageGlobalMove-PercentageWormInsert)/2d0)then
               if(grnd() >= Percentage4OperatorMove)then
                  call StepAdd(grnd() < PercentageOuterMove,Sector)
               else
                  if(.not.b_exch)then
                     call StepAdd4(Sector)
                  else
                     call StepFlavourchange_general()
                  endif
               endif

            !removal of hyb pair
            elseif(rand_num<=1d0-PercentageTauShiftMove-PercentageGlobalMove-PercentageWormInsert)then
               if(grnd() >= Percentage4OperatorMove)then
                  call StepRem(grnd() < PercentageOuterMove,Sector)
               else
                  if(.not.b_exch)then
                     call StepRem4()
                  else
                     call StepFlavourchange_general()
                  endif
               endif

            elseif(rand_num<=1d0-PercentageTauShiftMove-PercentageGlobalMove)then
               
               !attempt to go into iSector
               Sector = iSector

               !offdiagonals require mixed worm/hyb for Sector 2
               if(worm_offdiag) then
                   call StepWormHybAdd(Sector,iComponent)
               else
                   call StepWormAdd(Sector,iComponent)
               endif
               
            elseif (rand_num <= 1d0-PercentageTauShiftMove) then
               call StepGlob()
            else
               call StepShiftTau()
            endif
         !elseif(Sector .eq. iSector)
         else
            CntSampling(iSector) = CntSampling(iSector) + 1
            rand_num=grnd()
            !insertion of hyb pair
            if(rand_num<(1d0-PercentageWormReplace-PercentageWormInsert)/2d0)then
               if(grnd() >= Percentage4OperatorMove)then
                  call StepAdd(grnd() < PercentageOuterMove,Sector)
               else
                  if(.not.b_exch)then
                     call StepAdd4(Sector)
                  else
                     call StepFlavourchange_general()
                  endif
               endif
            !removal of of hyb pair
            elseif(rand_num<1d0-PercentageWormReplace-PercentageWormInsert)then
               if(grnd() >= Percentage4OperatorMove)then
                  call StepRem(grnd() < PercentageOuterMove,Sector)
               else
                  if(.not.b_exch)then
                     call StepRem4()
                  else
                     call StepFlavourchange_general()
                  endif
               endif
            elseif(rand_num<1d0-PercentageWormInsert) then
               if(Sector .ne. 6 .and. Sector .ne. 7) then
                  call StepWormReplace(Sector)
               endif
            !worm removal
            else
               !offdiagonals require mixed worm/hyb for Sector 2
               if(worm_offdiag) then
                   call StepWormHybRem(Sector)
               else
                   call StepWormRem(Sector)
               endif
            endif
         endif

         if (.not. b_segment) then
            if (modulo(sum(TryAdd) + sum(TryRem), NSlide) == 0) call shift_window(DTrace)
         end if

      enddo
   
      !after a sucessful runthrough we adjust eta
      if(CntSampling(iSector).eq.0) then 
         wormEta(iSector)=wormEta(iSector)*100d0
      elseif(CntSampling(1).eq.0) then
         wormEta(iSector)=wormEta(iSector)*0.01d0
      else
         wormEta(iSector)=wormEta(iSector)*CntSampling(1)/dble(CntSampling(iSector))
      endif
      
      
      !convergence criteria
      if(dble(abs(CntSampling(iSector)-CntSampling(1)))/dble((CntSampling(iSector)+CntSampling(1))).lt.tol) then
         Sector = 1
         call init_counters()
         !reset to initial configuraiton
         call clear_Trace(DTrace,DStates)
         call set_mc_config(size(taus), taus, orbs, spins, cas, hashybs, outer_sst, outer_state)
         deallocate(taus,orbs,spins,cas,hashybs)
         return
      endif
       
      !only on last iteration set wormEta to 0 if worm space is not entered
      if(i.eq.maxIter .and. CntSampling(iSector).eq.0) then
         wormEta(iSector)=0d0
      endif 

      Sector = 1
      call init_counters()
      !reset to initial configuraiton
      call clear_Trace(DTrace,DStates)
      call set_mc_config(size(taus), taus, orbs, spins, cas, hashybs, outer_sst, outer_state)
      
      call ptick(p)
   enddo
   
   deallocate(taus,orbs,spins,cas,hashybs)
   
end subroutine findEta


!!===============================================================================
!!===============================================================================
!!===============================================================================
!!! here come the segment moves
!!===============================================================================
!!===============================================================================
!!===============================================================================


!===============================================================================
!> compares the first significant figures of two large numbers
logical function equal_large(x1,x2)
!===============================================================================
!input
   real(KIND=KIND(0.D0)) :: x1,x2
!local
   integer :: e1,e2
   real(KIND=KIND(0.D0)) :: eps

   eps=0.000001_KINDR

   !write(*,*) " "
   !write(*,*) "equal_large"
   !!! exponents e1 from x1=a1*10^e1
   e1=floor(log10(x1))
   e2=floor(log10(x2))
   !write(*,*) "x1", x1
   !write(*,*) "x2", x2
   !write(*,*) "e1", e1
   !write(*,*) "e2", e2

   if(e1.eq.0.or.e2.eq.0)then

      if(abs(x1-x2).gt.eps)then
         equal_large=.false.
      else
         equal_large=.true.
      endif
      return

   endif
         
   if(abs(x1/10**e1-x2/10**e2).gt.eps)then
      equal_large=.false.
   else
      equal_large=.true.
   endif

end function equal_large 


!===============================================================================
!> This subroutine adds a pair of operators in the segment picture.
subroutine StepAdd_mine()
!===============================================================================
!input
   type(TStates),pointer         :: DStates
   type(TTrace),pointer          :: DTrace
   type(TStates_pointer)         :: pDStates
   type(TTrace_pointer)          :: pDTrace
!local
   integer                       :: CA,FPos(2),Ntmp
   real(KINDR),pointer           :: temp(:,:)
   type(TLogTr)                  :: TraceRat
   real(KINDR)                   :: S,rand,lmax
   real(KINDR)                   :: BosTraceNew,prob,signum
   type(TLogDet)                 :: DetNew, DetRat
   type(TOperPointer)            :: Oper(2),opertmp
   real(KINDR),pointer           :: Q(:),R(:)
   logical                       :: qn,seg,overbeta,check_qn_additionally
   real(kindr),pointer           :: FullHybr_offdiag(:,:)
   logical                       :: equal_times

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   ! Fortran 90 does not guarantee nullified pointers
   nullify(temp); nullify(Q); nullify(R)

   if(DTrace%iOperPool.eq.0)then
      do CA=1,2
         allocate(Oper(CA)%p)
!        allocate(Oper(CA)%p%normr(DTrace%NTruncStatesMax))
!        allocate(Oper(CA)%p%norml(DTrace%NTruncStatesMax))
         Oper(CA)%p%calc_new=.true.
      enddo
   else
      Oper(1)%p=>DTrace%OperPool(DTrace%iOperPool)%p
      Oper(2)%p=>DTrace%OperPool(DTrace%iOperPool-1)%p
      DTrace%OperPool(DTrace%iOperPool)%p=>null()
      DTrace%OperPool(DTrace%iOperPool-1)%p=>null()
      DTrace%iOperPool=DTrace%iOperPool-2
   endif
   Oper(1)%p%calc_new=.true.
   Oper(2)%p%calc_new=.true.
   Oper(1)%p%has_hyb=.true.
   Oper(2)%p%has_hyb=.true.
   oper(1)%p%next=>null()
   oper(1)%p%prev=>null()
   oper(2)%p%next=>null()
   oper(2)%p%prev=>null()
   oper(1)%p%fnext=>null()
   oper(1)%p%fprev=>null()
   oper(2)%p%fnext=>null()
   oper(2)%p%fprev=>null()
  
   !!! Propose a pair of operators to be inserted. lmax is the distance of the
   !!! creator to the next creator in trace. seg tells if segment or
   !!! antisegment, overbeta if the segment last over beta or not.
   !!! When a dangerous configuration gets proposed (one operator has been proposed to 
   !!! be on top of another one in the trace, and been slightly shifted), one 
   !!! has to check the validity of the configuration once more!
   qn=propose_insert_seg(DTrace,DStates,oper,lmax,seg,overbeta,check_qn_additionally)

   !!! check if proposed operators have equal time
   if(oper(1)%p%tau.eq.oper(2)%p%tau)then
      qn=.false.
   endif

   !!! if not allowed, exit
   if(.not.qn)then
      AccQNAdd=AccQNAdd+1
      do CA=1,2
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+2
      return
   endif

   !!! check if proposed operators have time equal to any of the other operators
   !!! in the configuration
   equal_times=check_equal_times(dtrace,oper(1)%p%tau,oper(2)%p%tau)
   !equal_times=.false.

   !!! if not allowed, exit
   if(equal_times)then
      write(*,*) "equal times!!"
      write(*,*) "oper(1)%p%tau", oper(1)%p%tau
      write(*,*) "oper(2)%p%tau", oper(2)%p%tau
      call print_trace_screen(dtrace)
      AccQNAdd=AccQNAdd+1
      do CA=1,2
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+2
      return
   endif

   !!! time order operators
   if(oper(2)%p%tau.lt.oper(1)%p%tau)then
      opertmp%p=>oper(1)%p
      oper(1)%p=>oper(2)%p
      oper(2)%p=>opertmp%p
   endif

   !if(check_qn_additionally)then
      !call my_insert_oper(dtrace,oper)
      !qn=check_qn_add_eachflavour(dtrace,oper(1)%p%orbital,oper(2)%p%spin)
      !call my_remove_oper(dtrace,oper)
      
      !if(.not.qn)then
         !AccQNAdd=AccQNAdd+1
         !do CA=1,2
            !DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
         !enddo
         !DTrace%iOperPool=DTrace%iOperPool+2
         !return
      !endif
   !endif

   !!! for test purposes: calculate full trace
   !!! we need no sign for the trace before adding, because in the weight only
   !!! the relative sign matters
   !tmp3=get_trace_seg_onestate(DTrace,DStates,oper)!*signum

   !!! insert pair
   call my_insert_oper(dtrace,oper)
   TryAdd(1)=TryAdd(1)+1
    
   !!! get position of operators in hybridisation matrix
   call get_fpos_insert(dtrace,oper,fpos)

   Oper(1)%p%has_hyb=.true.
   Oper(2)%p%has_hyb=.true.
   
   !!! calculate determinant ratio, and the column/row to be inserted in the
   !!! hybridisation matrix with Sherman-Morrison algorithm
   !DetRat=sig*get_DetRatPairAdd(DTrace,Oper,DStates,Q,R,S,FPos,sig)
   if(b_offdiag)then

      if(b_full_offdiag)then

         FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
         call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag, DetNew)
         DetRat = DetNew/DTrace%Det

      else

         !!! fast updates
         call add_ftau_offset(DTrace, DStates, Oper, FPos)
         DetRat = get_LogDetRatPairAdd_full(dtrace,oper,DStates,Q,R,S,FPos)
         DetNew = DetRat * DTrace%Det

      endif

   else

      DetRat = get_LogDetRatPairAdd(DTrace,DTrace%MInv,Oper,Q,R,S,FPos)
      DetNew = DetRat * DTrace%Det

   endif

   !!! calculate sign coming from the ordering of operators in a state:
   !!! c^dag_up c^dag_down |0> x |0> = c^dag_up |0> x |1> = |1> x |1>
   !!! whereas:
   !!! c^dag_down c^dag_up |0> x |0> = c^dag_down |1> x |0> = - |1> x |1>
   !!! In the Matrix-Vector Algorithm this sign is stored in the matrices 
   !!! of the Operators.

   signum=get_sign_seg(dtrace,oper)
   !!! for test purposes: calculate full trace
   !tmp2=get_trace_seg_onestate(DTrace,DStates,oper)*signum
   TraceRat = get_trace_seg_onestate_add(dtrace,DStates,oper,seg,overbeta)
   TraceRat%sign = TraceRat%sign * signum

   !if(.not.equal_large(tmp,tmp2/tmp3))then
      !write(*,*) "falsch!!!"
      !write(*,*) "tmp", tmp
      !write(*,*) "tmp2/tmp3", tmp2/tmp3
   !endif

   BosTraceNew=get_BosonicTrace(DTrace,DStates)
   rand=grnd()
   ntmp=dtrace%nosoper(oper(1)%p%orbital,oper(1)%p%spin)

   !!! the weight
   !if(b_full_trace)then
      !prob=dabs(DTrace%beta*lmax/(&
         !(dble(Ntmp)/2.0))*prefact&
         !*DetRat*TraceRat/DTrace%Trace*BosTraceNew/DTrace%BosonicTrace)
   !else
      prob=dabs(DTrace%beta*lmax/(&
         (dble(Ntmp)/2.0))&
         *wrat(TraceRat, DetRat)*BosTraceNew/DTrace%BosonicTrace)
   !endif

   if(rand.lt.prob)then
      !!! accept

      !!! Generate new hybridisation matrix and store it.
      !temp=>get_MatAdd(&
      !DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat,&
      !DTrace%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)/2-1,Q,R,S,FPos)
      !if(associated(DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat))then
         !deallocate(DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat)
      !endif
      !DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat=>temp

      if(.not.b_offdiag)then

         temp=>get_MatAdd(&
         DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat,&
         DTrace%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)/2-1,Q,R,S,FPos)
         if(associated(DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat))then
            deallocate(DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat)
         endif
         DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat=>temp

      else

         if(b_full_offdiag)then

            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>FullHybr_offdiag
            
         else

            temp=>get_MatAdd(&
            DTrace%MInv_full%Mat,&
            DTrace%Noper/2-1,Q,R,S,FPos)

            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>temp

         endif

      endif   
            
      !we save trace, determinant and sign of the current configuration
      DTrace%Det=DetNew
      DTrace%Trace = TraceRat * dtrace%trace
      !!! the outer state of the current segment configuration (segment or antisegment)
      dtrace%initial_outerstate_old(:,:)=dtrace%initial_outerstate(:,:)

      DTrace%cfgsign = -get_Sign(DTrace, DStates)
     
      DTrace%BosonicTrace=BosTraceNew
      AccAdd(1)=AccAdd(1)+1
   else
      !!! reject
      call my_remove_Oper(DTrace,Oper)
      !!! put operators back in operator pool
      do CA=1,size(Oper)
         if(associated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(associated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(associated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(associated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         !Oper(CA)%p%calc_new=.true.
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
      enddo
      DTrace%iOperPool=DTrace%iOperPool+2

      if(b_offdiag)then
      if(b_full_offdiag)then
         deallocate(FullHybr_offdiag)
      endif
      endif

   endif
   
   if(.not.b_full_offdiag)then
      deallocate(Q,R)
   endif
   
end subroutine StepAdd_mine


!===============================================================================
!> This subroutine removes a pair of operators in the segment picture.
subroutine StepRem_mine()
!===============================================================================
!input
   type(TStates),pointer         :: DStates
   type(TTrace),pointer          :: DTrace
   type(TStates_pointer)         :: pDStates
   type(TTrace_pointer)          :: pDTrace
!local
   integer                       :: CA, FPos(2), Ntmp
   real(KINDR),pointer           :: temp(:,:)
   type(TLogTr)                  :: TraceRat
   real(KINDR)                   :: S,rand
   real(KINDR)                   :: BosTraceNew,prob,signum,lmax
   type(TLogDet)                 :: DetNew, DetRat
   type(TOperPointer)            :: Oper(2),opertmp
!   type(TOper),pointer           :: Element
   logical                       :: qn,seg,overbeta
   real(kindr),pointer           :: FullHybr_offdiag(:,:)

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr
 
   nullify(temp)!; nullify(Element)

   TryRem(1)=TryRem(1)+1

   if(dtrace%noper.eq.0)then
      return
   endif

   !!! for test purposes: calculate full trace
   !tmp3=get_trace_seg_onestate(DTrace,DStates,oper)!*signum

   !!! Selects a segment or antisegment (with probability 0.5) from the configuration. Also gives the
   !!! lmax one would get when inserting this segment / antisegment.
   qn=propose_remove_seg(DTrace,DStates,Oper,2,lmax,seg,overbeta)

   if(qn.eqv..false.)then
      AccQNRem=AccQNRem+1
      return
   endif

   !!! time order operators
   if(oper(2)%p%tau.lt.oper(1)%p%tau)then
      opertmp%p=>oper(1)%p
      oper(1)%p=>oper(2)%p
      oper(2)%p=>opertmp%p
   endif

   !!! calculate sign (see StepAdd_mine for details).
   !!! before my_remove_oper, since the operators have still to be in the linked
   !!! list, in order to get their sign
   signum=get_sign_seg(dtrace,oper)
   call get_fpos_insert(dtrace,oper,fpos)
   call my_remove_Oper(dtrace,oper)

   !DetRat=sig*get_DetRatPairRem(DTrace,Oper,DStates,FPos,sig,S)
   if(.not.b_offdiag)then

      DetRat = get_LogDetRatPairRem(DTrace,Oper,FPos,S)
      DetNew = DetRat * DTrace%Det

   else

      if(b_full_offdiag)then

         FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
         call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag, DetNew)
         DetRat = DetNew/DTrace%Det

      else

         !!! fast updates
         call add_ftau_offset(dtrace,dstates,oper,fpos)
         DetRat = get_LogDetRatPairRem_full(DTrace,Oper,FPos,S)
         DetNew = DetRat * DTrace%Det

      endif

   endif

   !!! for test purposes: calculate full trace
   !tmp2=get_trace_seg_onestate(DTrace,DStates,oper)*signum
   TraceRat = get_trace_seg_onestate_rem(dtrace,DStates,oper,seg,overbeta)
   TraceRat%sign = TraceRat%sign * signum

   !if(.not.equal_large(tmp,tmp2/tmp3))then
      !write(*,*) "falsch!!!"
      !write(*,*) "tmp", tmp
      !write(*,*) "tmp2/tmp3", tmp2/tmp3
   !endif

   BosTraceNew=get_BosonicTrace(DTrace,DStates)
   rand=grnd()

   ntmp=dtrace%nosoper(oper(1)%p%orbital,oper(1)%p%spin)

   !if(b_full_trace)then
      !prob=dabs(&
         !dble((Ntmp+2.0)/2.0)&
         !/(DTrace%beta*lmax)*prefact*&
         !DetRat*TraceRat/DTrace%Trace*BosTraceNew/DTrace%BosonicTrace)
   !else
      prob=dabs(&
         dble((Ntmp+2.0)/2.0)&
         /(DTrace%beta*lmax)*&
         wrat(TraceRat, DetRat)*BosTraceNew/DTrace%BosonicTrace)
   !endif

   if(rand.lt.prob)then
      !!! accept

      if(dabs(detval(DetRat)).lt.1d-12) stop "DetRat too small in removal"
      !!! generate inverse hybridisation matrix of that flavour and store it
      !temp=>get_MatRem(&
         !DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat,&
         !DTrace%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)/2+1,FPos)
      !deallocate(DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat)
      !DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat=>temp

      if(.not.b_offdiag)then

         temp=>get_MatRem(&
         DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat,&
               DTrace%NOSOper(Oper(1)%p%Orbital,Oper(1)%p%Spin)/2+1,FPos)
         deallocate(DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat)
         DTrace%MInv(Oper(1)%p%Orbital,Oper(1)%p%Spin)%Mat=>temp
      else

         if(b_full_offdiag)then

            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            DTrace%MInv_full%Mat=>FullHybr_offdiag

         else

            !!! fast matrix generation
            temp=>get_MatRem(&
               DTrace%MInv_full%Mat,&
               DTrace%Noper/2+1,FPos)
            if(associated(DTrace%MInv_full%Mat))then
               deallocate(DTrace%MInv_full%Mat)
            endif
            !deallocate(temp)
            DTrace%MInv_full%Mat=>temp

         endif

      endif
      
      !!! we save the determinant and trace of the current configuration
      DTrace%Det=DetNew
      DTrace%Trace = TraceRat * dtrace%trace
      !!! the outer state of the current segment configuration (segment or antisegment)
      dtrace%initial_outerstate_old(:,:)=dtrace%initial_outerstate(:,:)

      !calculate sign of current configuration
      DTrace%cfgsign = -get_Sign(DTrace, DStates)
      
      DTrace%BosonicTrace=BosTraceNew
      
      !!! put operators back in operator pool
      do CA=1,size(Oper)
         if(associated(Oper(CA)%p%state))deallocate(Oper(CA)%p%state)
         if(associated(Oper(CA)%p%cache))deallocate(Oper(CA)%p%cache)
         if(associated(Oper(CA)%p%slogmax))deallocate(Oper(CA)%p%slogmax)
         if(associated(Oper(CA)%p%clogmax))deallocate(Oper(CA)%p%clogmax)
         DTrace%OperPool(DTrace%iOperPool+CA)%p=>Oper(CA)%p
         !Oper(CA)%p%calc_new=.true.
      enddo
      
      DTrace%iOperPool=DTrace%iOperPool+2
      AccRem(1)=AccRem(1)+1
   else
      !!! reject
      !!! put operator back in configuration
      call my_insert_Oper(DTrace,Oper)
!      Element=>DTrace%first
      !do while(associated(Element))
        !Element%calc_new=.false.
        !Element=>Element%next
      !enddo

      if(b_offdiag)then
      if(b_full_offdiag)then
         deallocate(FullHybr_offdiag)
      endif
      endif

   endif
   
end subroutine StepRem_mine


!===============================================================================
!> This subroutine measures the occupations
subroutine MeasOcc_seg()
!===============================================================================

   type(toper),pointer           :: element
   integer                       :: ib1, is1, ib2, is2
   type(tstates),pointer         :: dstates
   type(ttrace),pointer          :: dtrace
   type(tstates_pointer)         :: pdstates
   type(ttrace_pointer)          :: pdtrace
 
   pdstates = transfer(ipdstates, pdstates)
   pdtrace = transfer(ipdtrace, pdtrace)
   dstates => pdstates%ptr
   dtrace => pdtrace%ptr

   !===============================================================================
   ! singleoccs
   !===============================================================================

   !!! array to count length of segments
   dtrace%mu_accumulator(:,:)=0d0
   !!! array to store beginning of segments
   dtrace%mu_accumulator_tmp(:,:)=0d0

   !!! go through configuration and count segment lengths
   Element=>dtrace%first
   do while(associated(element))
      if(element%ca.eq.1)then
         dtrace%mu_accumulator_tmp(element%orbital,element%spin)=element%tau
      else
         dtrace%mu_accumulator(element%orbital,element%spin)=&
         dtrace%mu_accumulator(element%orbital,element%spin)+element%tau-dtrace%mu_accumulator_tmp(element%orbital,element%spin)
      endif
      Element=>Element%next
   enddo

   !!! add contribution from last operator to beta
   do ib1=1,dstates%nbands
   do is1=1,2
   if(dtrace%initial_outerstate_old(ib1,is1).eqv..true.)then
      dtrace%mu_accumulator(ib1,is1)=&
      dtrace%mu_accumulator(ib1,is1)+dtrace%beta-dtrace%mu_accumulator_tmp(ib1,is1)
   endif
   enddo
   enddo

   !!! save single occupancies
   single_occ=single_occ+dtrace%mu_accumulator
   do ib1=1,dstates%nbands
   do is1=1,2
   occ(ib1,is1,ib1,is1)=occ(ib1,is1,ib1,is1)+dtrace%mu_accumulator(ib1,is1)
   enddo
   enddo

!===============================================================================
! doubleoccs
!===============================================================================

   !!! array to count length of double occupations
   dtrace%int_accumulator(:,:,:,:)=0d0
   !!! array to store beginning of double occupations
   dtrace%int_accumulator_tmp(:,:,:,:)=0d0
   !!! store if at operator position is a segment in a specific flavour or antisegment
   dtrace%outerstate_tmp=dtrace%initial_outerstate_old


   Element=>dtrace%first
   do while(associated(element))

      !!! when creator is found, update position of beginning of double occupations
      if(element%ca.eq.1)then
         do ib1=1,dstates%nbands
         do is1=1,2
         if(dtrace%outerstate_tmp(ib1,is1))then
            !TODO: these objects are bosonic, but this symmetry is not used yet.
            !Does it matter?!
            dtrace%int_accumulator_tmp(element%orbital,element%spin,ib1,is1)=element%tau
            dtrace%int_accumulator_tmp(ib1,is1,element%orbital,element%spin)=element%tau
         endif
         enddo
         enddo
         dtrace%outerstate_tmp(element%orbital,element%spin)=.true.
      endif

      !!! when a annihilator is found
      if(element%ca.eq.2)then
         dtrace%outerstate_tmp(element%orbital,element%spin)=.false.
         do ib1=1,dstates%nbands
         do is1=1,2
         !!! and the state of an other flavour is occupied
         if(dtrace%outerstate_tmp(ib1,is1))then
            !!! then accumulate the double occupataion that ends here
            dtrace%int_accumulator(element%orbital,element%spin,ib1,is1)=&
            dtrace%int_accumulator(element%orbital,element%spin,ib1,is1)+&
            element%tau-dtrace%int_accumulator_tmp(element%orbital,element%spin,ib1,is1)
            dtrace%int_accumulator(ib1,is1,element%orbital,element%spin)=&
            dtrace%int_accumulator(ib1,is1,element%orbital,element%spin)+&
            element%tau-dtrace%int_accumulator_tmp(ib1,is1,element%orbital,element%spin)
         endif
         enddo
         enddo

      endif
         
      Element=>Element%next
   enddo

   !!! treat intervall from last operator to beta
   do ib1=1,dstates%nbands
   do is1=1,2
   do ib2=1,dstates%nbands
   do is2=1,2
   if((ib1.eq.ib2).and.(is1.eq.is2))cycle
   if(dtrace%initial_outerstate_old(ib1,is1).and.dtrace%initial_outerstate_old(ib2,is2))then
      dtrace%int_accumulator(ib1,is1,ib2,is2)=dtrace%int_accumulator(ib1,is1,ib2,is2)+(dtrace%beta-dtrace%int_accumulator_tmp(ib1,is1,ib2,is2))/2.0
      dtrace%int_accumulator(ib2,is2,ib1,is1)=dtrace%int_accumulator(ib2,is2,ib1,is1)+(dtrace%beta-dtrace%int_accumulator_tmp(ib2,is2,ib1,is1))/2.0
   endif
   enddo
   enddo
   enddo
   enddo

   !!! save double occupations
   double_occ=double_occ+dtrace%int_accumulator
   do ib1=1,dstates%nbands
   do is1=1,2
   do ib2=1,dstates%nbands
   do is2=1,2
   occ(ib1,is1,ib2,is2)=occ(ib1,is1,ib2,is2)+dtrace%int_accumulator(ib1,is1,ib2,is2)
   enddo
   enddo
   enddo
   enddo

end subroutine MeasOcc_seg


!===============================================================================
!> Measuring the density-density correlation function <n(tau)n(0)> in a
!primitive way.
!!! CAUTIION: do only take a few hundred tau points, since this is a doubly
!!! nested loop!
subroutine MeasSusz_seg()
!===============================================================================
   type(tstates),pointer         :: dstates
   type(ttrace),pointer          :: dtrace
   type(tstates_pointer)         :: pdstates
   type(ttrace_pointer)          :: pdtrace

   type(toper),pointer           :: e1, e2
   integer                       :: ib1, is1, ib2, is2, ig1, ig2, di
   real(KINDR)                   :: dt, dt2, tau, tau2 !, taup, taup2
   logical :: seg1,seg2,full1,full2,normal,fullboth

   !www=.true.
   !www=.false.

   pdstates = transfer(ipdstates, pdstates)
   pdtrace = transfer(ipdtrace, pdtrace)
   dstates => pdstates%ptr
   dtrace => pdtrace%ptr
   
   !www=.true.

   !if(www)write(*,*) " "
   !if(www)write(*,*) "meassusz"

   do ib1=1,dstates%nbands
   do is1=1,2
   do ib2=1,dstates%nbands
   do is2=1,2
   !if(www)write(*,*) "ib1", ib1
   !if(www)write(*,*) "is1", is1
   !if(www)write(*,*) "ib2", ib2
   !if(www)write(*,*) "is2", is2

      !!! define which case
      !if((ib1.eq.ib2).and.(is1.eq.is2))cycle
      if((dtrace%nosoper(ib1,is1).eq.0).and.(dtrace%initial_outerstate_old(ib1,is1).eqv..false.).and.&
         (dtrace%nosoper(ib2,is2).eq.0).and.(dtrace%initial_outerstate_old(ib2,is2).eqv..false.))cycle

      full1=.false.
      full2=.false.
      normal=.false.
      fullboth=.false.

      if((dtrace%nosoper(ib1,is1).eq.0).and.(dtrace%initial_outerstate_old(ib1,is1).eqv..true.).and.&
         (dtrace%nosoper(ib2,is2).ne.0))then
         full1=.true.
      endif
      if((dtrace%nosoper(ib2,is2).eq.0).and.(dtrace%initial_outerstate_old(ib2,is2).eqv..true.).and.&
         (dtrace%nosoper(ib1,is1).ne.0))then
         full2=.true.
      endif
      if((dtrace%nosoper(ib1,is1).eq.0).and.(dtrace%initial_outerstate_old(ib1,is1).eqv..true.).and.&
         (dtrace%nosoper(ib2,is2).eq.0).and.(dtrace%initial_outerstate_old(ib2,is2).eqv..true.))then
         fullboth=.true.
      endif
      if((dtrace%nosoper(ib1,is1).ne.0).and.(dtrace%nosoper(ib2,is2).ne.0))then
         normal=.true.
      endif
      !if(www)write(*,*) "full1", full1
      !if(www)write(*,*) "full2", full2
      !if(www)write(*,*) "fullboth", fullboth
      !if(www)write(*,*) "normal", normal

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!! normal case !!!!!!!!!!
      if(normal)then

         !if(www)write(*,*) "normal"

         !if(www)call print_trace_flavour(dtrace,dstates,ib1,is1)
         !if(www)call print_trace_flavour(dtrace,dstates,ib2,is2)
         e1=>dtrace%ffirst(ib1,is1)%p

         if(e1%ca.eq.1)then
            seg1=.false.
         else
            seg1=.true.
         endif
!         taup=0d0
         tau=e1%tau
         !if(www)write(*,*) "--> taup,tau", taup,tau
         !if(www)write(*,*) "seg1", seg1

         do ig1=0,ngtau-1

            dt=nn_disc(ig1)
            if(dt.gt.tau)then
               dl1: do while(dt.gt.tau)
               !if(www)write(*,*) "dowhile: dt, tau", dt, tau
               if(associated(e1%fnext))then
                  e1=>e1%fnext
!                  taup=tau
                  tau=e1%tau
                  if(seg1.eqv..true.)then
                     seg1=.false.
                  else
                     seg1=.true.
                  endif
                  !if(www)write(*,*) "--> taup,tau", taup,tau
                  !if(www)write(*,*) "seg1", seg1
               else
!                  taup=tau
                  tau=dtrace%beta
                  e1=>null()
                  if(seg1.eqv..true.)then
                     seg1=.false.
                  else
                     seg1=.true.
                  endif
                  exit dl1
                  !if(www)write(*,*) "--> taup,tau", taup,tau
                  !if(www)write(*,*) "seg1", seg1
               endif
               enddo dl1
            endif
            !if(www)write(*,*) "ig1,dt", ig1,dt

            !!! meas correlation funciton

            if(seg1.eqv..false.)cycle

            !if(www)write(*,*) "-----------------------------------"
            !if(www)write(*,*) "filled at seg1, run:"

            e2=>dtrace%ffirst(ib2,is2)%p

            if(e2%ca.eq.1)then
               seg2=.false.
            else
               seg2=.true.
            endif

!            taup2=0d0
            tau2=e2%tau
            !if(www)write(*,*) "--> taup2,tau2", taup2,tau2
            !if(www)write(*,*) "seg2", seg2

            do ig2=0,ngtau-1
               !if(abs(ig2-ig1).gt.50.and.abs(ig2-ig1+dtrace%beta).gt.50)cycle

               dt2=nn_disc(ig2)
               if(dt2.gt.tau2)then
                  dl2: do while(dt2.gt.tau2)
                  !if(www)write(*,*) "dowhile: dt2, tau", dt2, tau2
                  if(associated(e2%fnext))then
                     e2=>e2%fnext
!                     taup2=tau2
                     tau2=e2%tau
                     if(seg2.eqv..true.)then
                        seg2=.false.
                     else
                        seg2=.true.
                     endif
                     !if(www)write(*,*) "--> taup2,tau2", taup2,tau2
                     !if(www)write(*,*) "seg2", seg2
                  else
!                     taup2=tau2
                     tau2=dtrace%beta
                     e2=>null()
                     if(seg2.eqv..true.)then
                        seg2=.false.
                     else
                        seg2=.true.
                     endif
                     exit dl2
                     !if(www)write(*,*) "--> taup2,tau2", taup2,tau2
                     !if(www)write(*,*) "seg2", seg2
                  endif
                  enddo dl2
               endif
               !if(www)write(*,*) "ig2,dt2", ig2,dt2

               if(seg2)then
                  if(ig1.le.ig2)then
                     di=ig2-ig1
                  elseif(ig1.gt.ig2)then
                     di=ngtau+ig2-ig1
                  endif
                  ntau_n0(ib1,is1,ib2,is2,di)=ntau_n0(ib1,is1,ib2,is2,di)+1
                  if((di.eq.ngtau))then
                     write(*,*) "1eq!!!"
                     write(*,*) "ig1", ig1
                     write(*,*) "ig2", ig2
                  endif
                  !if((di.lt.1))then
                     !write(*,*) "lt!!!"
                  !endif
               endif

            enddo

            !!! !!!!!

         enddo

      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(fullboth)then
         !if(www)write(*,*)"fullboth"

         do ig1=0,ngtau-1
         !!! meas correlation funciton

         do ig2=0,ngtau-1
         !if(abs(ig2-ig1).gt.50.and.abs(ig2-ig1+dtrace%beta).gt.50)cycle

         dt2=nn_disc(ig2)

            if(ig1.le.ig2)then
               di=ig2-ig1
            elseif(ig1.gt.ig2)then
               di=ngtau+ig2-ig1
            endif
            ntau_n0(ib1,is1,ib2,is2,di)=ntau_n0(ib1,is1,ib2,is2,di)+1
            if((di.gt.ngtau))then
               write(*,*) "4eq!!!"
               write(*,*) "ig1", ig1
               write(*,*) "ig2", ig2
            endif
            !if((di.lt.1))then
               !write(*,*) "lt!!!"
            !endif

         enddo
         enddo

         !!! !!!!!
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!! one full !!!!!!!!!!
      if(full1)then
         !if(www)write(*,*)"full1"

         do ig1=0,ngtau-1
         !!! meas correlation funciton

         !if(www)write(*,*) "-----------------------------------"
         !if(www)write(*,*) "filled at seg1, run:"

         e2=>dtrace%ffirst(ib2,is2)%p

         if(e2%ca.eq.1)then
            seg2=.false.
         else
            seg2=.true.
         endif

!         taup2=0d0
         tau2=e2%tau
         !if(www)write(*,*) "--> taup2,tau2", taup2,tau2
         !if(www)write(*,*) "seg2", seg2

         do ig2=0,ngtau-1
         !if(abs(ig2-ig1).gt.50.and.abs(ig2-ig1+dtrace%beta).gt.50)cycle

            dt2=nn_disc(ig2)
            if(dt2.gt.tau2)then
               dl3: do while(dt2.gt.tau2)
               if(associated(e2%fnext))then
                  e2=>e2%fnext
!                  taup2=tau2
                  tau2=e2%tau
                  if(seg2.eqv..true.)then
                     seg2=.false.
                  else
                     seg2=.true.
                  endif
                  !if(www)write(*,*) "--> taup2,tau2", taup2,tau2
                  !if(www)write(*,*) "seg2", seg2
               else
!                  taup2=tau2
                  tau2=dtrace%beta
                  e2=>null()
                  if(seg2.eqv..true.)then
                     seg2=.false.
                  else
                     seg2=.true.
                  endif
                  !if(www)write(*,*) "--> taup2,tau2", taup2,tau2
                  !if(www)write(*,*) "seg2", seg2
                  exit dl3
               endif
               enddo dl3
            endif
            !if(www)write(*,*) "ig2,dt2", ig2,dt2

            if(seg2)then
               if(ig1.le.ig2)then
                  di=ig2-ig1
               elseif(ig1.gt.ig2)then
                  di=ngtau+ig2-ig1
               endif
               ntau_n0(ib1,is1,ib2,is2,di)=ntau_n0(ib1,is1,ib2,is2,di)+1
               if((di.gt.ngtau))then
                  write(*,*) "2eq!!!"
                  write(*,*) "ig1", ig1
                  write(*,*) "ig2", ig2
               endif
               !if((di.lt.1))then
                  !write(*,*) "lt!!!"
               !endif
            endif

         enddo
         enddo

         !!! !!!!!

      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!! one full !!!!!!!!!!
      if(full2)then
         !if(www)write(*,*)"full1"

         do ig1=0,ngtau-1
         !!! meas correlation funciton

         !if(www)write(*,*) "-----------------------------------"
         !if(www)write(*,*) "filled at seg1, run:"

         e1=>dtrace%ffirst(ib1,is1)%p

         if(e1%ca.eq.1)then
            seg1=.false.
         else
            seg1=.true.
         endif

!         taup=0d0
         tau=e1%tau
         !if(www)write(*,*) "--> taup,tau", taup,tau
         !if(www)write(*,*) "seg1", seg1

         do ig2=0,ngtau-1

            dt2=nn_disc(ig2)
            if(dt2.gt.tau)then
               dl4: do while(dt2.gt.tau)
               if(associated(e1%fnext))then
                  e1=>e1%fnext
!                  taup=tau
                  tau=e1%tau
                  if(seg1.eqv..true.)then
                     seg1=.false.
                  else
                     seg1=.true.
                  endif
                  !if(www)write(*,*) "--> taup,tau", taup,tau
                  !if(www)write(*,*) "seg1", seg1
               else
!                  taup=tau
                  tau=dtrace%beta
                  e1=>null()
                  if(seg1.eqv..true.)then
                     seg1=.false.
                  else
                     seg1=.true.
                  endif
                  !if(www)write(*,*) "--> taup,tau", taup,tau
                  !if(www)write(*,*) "seg1", seg1
                  exit dl4
               endif
               enddo dl4
            endif
            !if(www)write(*,*) "ig2,dt2", ig2,dt2

            if(seg1)then
               if(ig1.le.ig2)then
                  di=ig2-ig1
               elseif(ig1.gt.ig2)then
                  di=ngtau+ig2-ig1
               endif
               ntau_n0(ib1,is1,ib2,is2,di)=ntau_n0(ib1,is1,ib2,is2,di)+1
               if((di.gt.ngtau))then
                  write(*,*) "3eq!!!"
                  write(*,*) "ig1", ig1
                  write(*,*) "ig2", ig2
               endif
               !if((di.lt.1))then
                  !write(*,*) "lt!!!"
               !endif
            endif

         enddo
         enddo

         !!! !!!!!

      endif
   enddo
   enddo
   enddo
   enddo

   !www=.false.

end subroutine MeasSusz_seg

!===============================================================================
!> Globalmoves for segment algorithm
subroutine StepGlob_mine()
!===============================================================================
!input
   type(TStates),pointer      :: DStates
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
!local
   type(TLogTr)                  :: tr1, tr2
   real(KINDR)                   :: rand, prob
   real(KINDR)                   :: BosTraceNew
   type(TLogDet)                 :: DetNew, DetRat
   type(TSubMatrix),pointer      :: FullHybr(:,:)
   integer                       :: iB, iS
   integer :: update,b1,s1,b2,s2
   type(TOperPointer)            :: Oper(2)

   real(kindr),pointer            :: FullHybr_offdiag(:,:)

   oper(1)%p=>null()
   oper(2)%p=>null()

   rand=grnd()
   if(rand.lt.0.9D0)then
      update=1 !!! spinflip   
   else
      update=2 !! exchange two flavours
   endif
   !update=1

   pDStates=transfer(ipDStates,pDStates)
   pDTrace=transfer(ipDTrace,pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   nullify(FullHybr)

   TryGlob=TryGlob+1
   if(DTrace%NOper.eq.0) then
      return
   endif
   
   if(update.eq.1)then

      !!! perform spin-flip in all orbitals
      !tr1=get_trace_seg_onestate(DTrace,DStates,oper)*sign(1d0,dtrace%trace)
      tr1=get_trace_seg_onestate(DTrace,DStates,oper)
      call SpinFlipUpdate_mine(DTrace,DStates)
      tr2=get_trace_seg_onestate(DTrace,DStates,oper)

   elseif(update.eq.2)then

      !! permute 2 flavours randomly
      b1=ceiling(grnd()*dble(NBands))
      s1=ceiling(grnd()*dble(2))
      if(b1.eq.0)b1=1
      if(s1.eq.0)s1=1
      b2=ceiling(grnd()*dble(NBands))
      s2=ceiling(grnd()*dble(2))
      if(b2.eq.0)b2=1
      if(s2.eq.0)s2=1
      !b1=1
      !b2=1
      !s1=1
      !s2=2
      if((b1.eq.b2).and.(s1.eq.s2))then
         return
      endif

      tr1=get_trace_seg_onestate(DTrace,DStates,oper)
      call FlavourExchangeUpdate(dtrace,b1,s1,b2,s2)
      tr2=get_trace_seg_onestate(DTrace,DStates,oper)

   endif

   !!! calculate new hybridisation matrix and determinant
   if(b_offdiag)then
      FullHybr_offdiag => get_FullHybr_offdiag(DTrace, DStates)
      call get_MatLogDetFull(size(FullHybr_offdiag(1,:)), FullHybr_offdiag, DetNew)
      DetRat = DetNew/DTrace%Det
   else
      FullHybr=>get_InvFullHybr(DTrace,DStates,DetNew)
      DetRat=DetNew*get_FullLogDet(DTrace,DStates)
   endif

   BosTraceNew=get_BosonicTrace(DTrace,DStates)

   prob=abs(wrat(tr2/tr1, DetRat)*BosTraceNew/DTrace%BosonicTrace)
   if(grnd().lt.prob)then
      !!! accepted

      DTrace%Det=DetNew
      !DTrace%Trace=tr2
      DTrace%cfgsign = -get_Sign(DTrace, DStates)
      !!! this gets updated by the function, that applies the global update to
      !!! the linked list
      !dtrace%initial_outerstate_old(:,:)=dtrace%initial_outerstate(:,:)

      if(b_offdiag)then

         if(associated(DTrace%MInv_full%Mat))then
            deallocate(DTrace%MInv_full%Mat)
         endif
         DTrace%MInv_full%Mat=>FullHybr_offdiag

      else

         do iB=1,NBands
         do iS=1,2
            deallocate(DTrace%MInv(iB,iS)%Mat)
         enddo
         enddo
         deallocate(DTrace%MInv)
         DTrace%MInv=>FullHybr

      endif

      AccGlob=AccGlob+1

   else
      !!! rejected
      !!! revert
      if(update.eq.1)then
         call SpinFlipUpdate_mine(DTrace,DStates)
      else
         call FlavourExchangeUpdate(dtrace,b1,s1,b2,s2)
      endif

      if(b_offdiag)then

         deallocate(FullHybr_offdiag)

      else

         do iB=1,NBands
         do iS=1,2
            deallocate(FullHybr(iB,iS)%Mat)
         enddo
         enddo
         deallocate(FullHybr)

      endif

   endif

end subroutine StepGlob_mine


!!===============================================================================
!real(KINDR) function get_GLin(AGtau,Atau)
!!===============================================================================
!! This function takes an arbitrary tau point and the Gtau array given on the
!! tau-grid defined by beta and NTau and returns the linearized Gtau at this
!! tau point.
!!input
!   real(KINDR)                :: AGtau(:)
!   real(KINDR)                :: Atau
!!local
!   real(KINDR)                :: taum,taup,eps
!   integer                    :: idtm,idtp
!
!!  We introduce the small number eps to avoid divisions by zero
!!  and numerical problems in generating idtm.
!   eps=1d-12
!   idtm=int((Atau+eps)/beta*dble(NGTau))
!   idtp=idtm+1
!   if(idtm.eq.0)idtm=1
!   taum=beta/dble(NGTau)*(dble(idtm)-0.5)
!   taup=beta/dble(NGTau)*(dble(idtp)-0.5)
!!  write(*,*)tau,(tau+eps)/beta*dble(NTau),idtm,taum,idtp,taup
!   get_GLin=AGtau(idtm)+(AGtau(idtp)-AGtau(idtm))/(taup-taum+eps)*(Atau-taum)
!end function get_GLin

!===============================================================================
!subroutine buildGGTau()
!   !===============================================================================
!   ! Here we generate the not fully connected part (=GG) of the 2 particle GF.
!   ! We use the function get_GLin since we build expressions like G(tau1-tau2),
!   ! where the argument tau1-tau2 is off the tau-grid given by beta and NGtau.
!   !local
!   integer                       :: i,j,k,idt1,idt2,idt3
!   integer                       :: iB,iS,iB1,iB2,iS1,iS2,NBin,indx
!   real(KINDR)                   :: tau,tau1,tau2,tau3,beta,A,B,C,D
!   integer                       :: NTau,NBands
!
!   beta=beta
!   NTau=NGtau
!   NBands=NBands
!
!   !  A,B,C,D are the different contractions:
!   !  A=G(tau1),B=G(tau3-tau2),C=G(tau3),D=G(tau1-tau2)
!   A=0d0
!   B=0d0
!   C=0d0
!   D=0d0
!   do iB1=1,NBands
!   do iB2=1,NBands
!   do iS1=1,2
!   do iS2=1,2
!   ! The maximum rank for Fortran arrays is 7, so here we use for (Spin1,Spin2)=(i,j) the
!   ! spin index convention s=2*int(s1/2)+s2: (1,1)=1, (1,2)=2, (2,1)=3, (2,2)=4
!      iS=2*int(iS1/2)+iS2
!      indx=0
!   do i=1,NTau
!      tau1=beta/dble(Ntau)*(dble(i)-0.5)
!      idt1=int(tau1/beta*dble(Ntau))+1
!   do j=1,NTau
!      tau2=beta/dble(Ntau)*(dble(j)-0.5)
!      idt2=int(tau2/beta*dble(Ntau))+1
!   !  write(40000+iB1*1000+iB2*100+iS1*10+iS2,'(A,I6)')'index=',indx
!   do k=1,NTau
!      tau3=beta/dble(Ntau)*(dble(k)-0.5)
!      idt3=int(tau3/beta*dble(Ntau))+1
!
!      ! We mirrow the GF about beta/2 to meet the Toschi convention.
!      ! This is done by writing G(beta-tau) instead of G(tau).
!      A=get_GLin(Gtau(iB1,iS1,:,1),NTau,beta,beta-tau1)
!
!      if(tau3-tau2.gt.0)then
!         B=get_GLin(Gtau(iB2,iS2,:,1),NTau,beta,beta-(tau3-tau2))
!      elseif(tau3-tau2.lt.0)then
!         B=-get_GLin(Gtau(iB2,iS2,:,1),NTau,beta,beta-(tau3-tau2+beta))
!      elseif(tau3-tau2.eq.0)then
!         B=occ(iB2+(iS2-1)*NBands,iB2+(iS2-1)*NBands,1)
!      endif
!
!      if(tau1-tau2.gt.0)then
!         C=get_GLin(Gtau(iB1,iS1,:,1),NTau,beta,beta-(tau1-tau2))
!      elseif(tau1-tau2.lt.0)then
!         C=-get_GLin(Gtau(iB1,iS1,:,1),NTau,beta,beta-(tau1-tau2+beta))
!      elseif(tau1-tau2.eq.0)then
!         C=-occ(iB1+(iS1-1)*NBands,iB1+(iS1-1)*NBands,1)
!      endif
!
!      D=get_GLin(Gtau(iB2,iS2,:,1),NTau,beta,beta-tau3)
!
!      if(iB1.eq.iB2.and.iS1.eq.iS2)then
!         ! In the spin and band diagonal case we have:
!         ! GG(tau1,tau2,tau3)=G(tau1)G(tau3-tau2)-G(tau3)G(tau1-tau2)
!         GGtau(iB1,iB2,iS,idt1,idt2,idt3,1)=A*B-C*D
!      else
!         ! In the spin and band off-diagonal case we have:
!         ! GG(tau1,tau2,tau3)=-G(tau3)G(tau1-tau2)
!         GGtau(iB1,iB2,iS,idt1,idt2,idt3,1)=-C*D
!      endif
!
!      ! write(40000+iB1*1000+iB2*100+iS1*10+iS2,'(3(F8.3),2(E15.6))')&
!      !   tau1,tau2,tau3,A,B,C,D,GGtau(iB1,iB2,iS,idt1,idt2,idt3,1)
!      !   tau1,tau2,tau3,GGtau(iB1,iB2,iS,idt1,idt2,idt3,1),G4tau(iB1,iB2,iS,idt1,idt2,idt3,1)
!   enddo
!      indx=indx+1
!      ! write(40000+iB1*1000+iB2*100+iS1*10+iS2,*)
!      ! write(40000+iB1*1000+iB2*100+iS1*10+iS2,*)
!   enddo
!   enddo
!   enddo
!   enddo
!   enddo
!   enddo
!end subroutine  buildGGtau

!===============================================================================
subroutine init_solver(u_matrix_in,Ftau_full,muimp_full,&
                       screening_function,iter_no,Nftau,NBands)
!===============================================================================
   use signals

!input
   integer                    :: iter_no
   integer                    :: Nftau,NBands
   real(KINDR)                :: u_matrix_in(2*NBands,2*NBands,2*NBands,2*NBands)
   real(KINDR)                :: Ftau(NBands,2,Nftau)
   real(KINDR)                :: Ftau_full(NBands,2,NBands,2,Nftau)
   real(KINDR)                :: muimp(NBands,2)
   real(KINDR)                :: muimp_full(NBands,2,NBands,2)
   real(KINDR)                :: screening_function(NBands,2,NBands,2,Nftau)
!local
   type(TStates), pointer     :: DStates
   type(TTrace), pointer      :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   real(KINDR), allocatable   :: HEValues(:)
   type(TPsis)                :: DPsis, DTransformed_Psis
   type(TOperator)            :: DH
   type(TOperator)            :: HEVectors
   integer                    :: i
   integer                    :: ib, is
   integer, allocatable       :: states2substates(:)
   integer                    :: nsubstates
   
   allocate(pDStates%ptr)
   allocate(pDTrace%ptr)
   ipDStates=transfer(pDStates,ipDStates)
   ipDTrace=transfer(pDTrace,ipDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr
   
   allocate(u_matrix(2*NBands,2*NBands,2*NBands,2*NBands))
   u_matrix = u_matrix_in

   do ib=1,NBands
   do is=1,2
      Ftau(ib,is,:)=Ftau_full(ib,is,ib,is,:)
      muimp(ib,is)=muimp_full(ib,is,ib,is)
   enddo
   enddo

   call init_States(DStates)
   call qns2substates(DStates, nsubstates, states2substates)
   call init_SubStates(DStates, nsubstates, states2substates)
   call init_Psi(DPsis,DStates)
   !call init_Hamiltonian(DH,DStates,muimp)
   call oneparticle_hamiltonian(DH,DStates,DPsis,muimp_full)
   !call print_Operator(DH,DStates) 
   call u_hamiltonian(u_matrix, DH, DStates, DPsis)
   call force_hermitian(DH, DStates)
   !call print_Operator(DH,DStates) 
   !start now we analyze the hamiltonian for block diagonal form
   if (index(get_String_Parameter("QuantumNumbers"),"All").ne.0) then
    call analyze_hamiltonian(DH, DStates, nsubstates, states2substates)
    call dest_Psis(DPsis)
    call dest_TOperator(DH)
    call dest_States(DStates)
    call init_States(DStates)
    call init_SubStates(DStates, nsubstates, states2substates)
!   call print_SubStates(DStates)
    call init_Psi(DPsis, DStates)
    !call init_Hamiltonian(DH, DStates, muimp)
    call oneparticle_hamiltonian(DH,DStates,DPsis,muimp_full)
    call u_hamiltonian(u_matrix, DH, DStates, DPsis)
    call force_hermitian(DH, DStates)
!    call print_Operator(DH,DStates,102) 
   endif
! end of analyzing hamiltonian

   call init_TOperator(HEVectors,DStates)
   allocate(HEValues(0:DStates%NStates-1))
   call diag_Operator(DH,DStates,HEVectors,HEValues)
   call set_Psis(DPsis,DStates)
   call set_Hamiltonian(DH,DStates)
   !call print_Operator(DH,Dstates)
   !call print_SubStates(DStates)
   call set_HEVal(HEValues,DStates)
   !write(*,*) "hevalues", HEValues
   call set_HEVec(HEVectors,DStates)
   !write(*,*) "eigenvektoren"
   !call print_Operator(hevectors,Dstates)

   !write(*,*) "states"
   !call print_states(DStates)
   !write(*,*) "substates"
   !call print_substates(DStates)

   if(simid.eq.0)then
   if(get_Integer_Parameter("PrintDensityMatrixBasis").ne.0)then
      call print_densitymatrix_basis(Dstates,iter_no)
   endif
   endif

   !do i = 0,DStates%NStates-1
      !write (*, *)  Hevalues(i)
   !enddo
   !call print_ev(DSTates)

   call transform_Psis(DStates,DTransformed_Psis)
   call set_EB_Psis(DTransformed_Psis,DStates)   
   !call print_psi(DTransformed_Psis,Dstates)
   !stop

   deallocate(HEValues)
   call dest_Psis(DPsis)
   call dest_Psis(DTransformed_Psis)
   call dest_TOperator(DH)
   call dest_TOperator(HEVectors)

   call init_CTQMC()

   do i = 0, DStates%NSStates - 1
      StatesSuperstates(DStates%SubStates(i)%Offset:&
                        DStates%SubStates(i)%Offset + DStates%SubStates(i)%NStates - 1) = i
      EigenstatesEnergies(DStates%SubStates(i)%Offset:&
                          DStates%SubStates(i)%Offset + DStates%SubStates(i)%NStates - 1) =&
         DStates%SubStates(i)%Eval(:)
      do ib = 0, NBands - 1
         do is = 0, 1
            OccbasisMapping(DStates%SubStates(i)%Offset:&
                            DStates%SubStates(i)%Offset + DStates%SubStates(i)%NStates - 1, ib + 1, is + 1) =&
               ibits(DStates%SubStates(i)%States(:), ib + is * NBands, 1)
         end do
      end do
   end do

   call init_Trace(DTrace,DStates,Ftau,Ftau_full,screening_function,Nftau,muimp,u_matrix)
   if (.not. b_segment) NSlide = max(2 * NBands, floor(dble(NCorr) / dble(DTrace%num_windows)))
   call ft_init()
  
   
   !write (0,"('FourPoint: ',I3)") FourPnt
   write(simstr,'(I4)') SimID

   call register_trap(SIGNAL_TERM, trap_notify)
   call register_trap(SIGNAL_INT,  trap_notify)
   call register_trap(SIGNAL_USR1, trap_notify)

   !FIXME: functionality
   if((b_segment.eqv..true.).and.(PercentageWormInsert.ne.0))then
      stop "Worm and Segment together not implemented!"
   endif

end subroutine init_solver

!===============================================================================
subroutine ctqmc_calibrate(Ncalibrate, list)
!===============================================================================
   use type_progress

!local
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   type(progress)             :: p
   integer(c_int64_t)         :: i, Ncalibrate
   logical                    :: list
   real(KINDR)                :: rand_num
   real(KINDR)                :: time_start, time_end

   pDStates = transfer(ipDStates, pDStates)
   pDTrace = transfer(ipDTrace, pDTrace)
   DTrace => pDTrace%ptr

   if (list) then
      if (.not. allocated(AccPairTau)) allocate(AccPairTau(1000))
   else
      if (allocated(AccPairTau)) deallocate(AccPairTau)
      if (allocated(AccPair)) deallocate(AccPair)
      allocate(AccPair(NBands, 2, 1000))
      AccPair = 0
      NAccPair = 1000
      AccPairMax = DTrace%taudiff_max
   end if

   call cpu_time(time_start)

   p%adaptrate = 0.1
   call pstart(p, int(Ncalibrate,PINT), &
               title='WinCbr.' // simstr // ':', fancy=fancy_prog)

   do i = 1, Ncalibrate
      if (any_signal_fired) exit
      rand_num = grnd()
      if (rand_num < (1.0_KINDR - PercentageTauShiftMove-PercentageGlobalMove)/2.0_KINDR) then
         ! add pair(s)
         if (grnd() >= Percentage4OperatorMove) then
            call StepAdd(grnd() < PercentageOuterMove, 1)
         else
            if (.not. b_exch) then
               call StepAdd4(1)
            else
               call StepFlavourchange_general()
            end if
         end if
      else if (rand_num <= 1.0_KINDR - PercentageTauShiftMove - PercentageGlobalMove)then
         ! remove pair(s)
         if (grnd() >= Percentage4OperatorMove) then
            call StepRem(grnd() < PercentageOuterMove, 1)
         else
            if (.not. b_exch) then
               call StepRem4()
            else
               call StepFlavourchange_general()
            end if
         end if
      else if (rand_num <= 1.0_KINDR - PercentageTauShiftMove) then
         call StepGlob()
      else
         call StepShiftTau()
      end if
      call ptick(p)
   end do

   call cpu_time(time_end)
   time_calibration = time_calibration + (time_end - time_start)
end subroutine ctqmc_calibrate

subroutine set_taudiffmax(taudiffmax)
   real(KINDR), intent(in) :: taudiffmax
   type(TTrace), pointer   :: DTrace
   type(TTrace_pointer)    :: pDTrace

   pDTrace = transfer(ipDTrace, pDTrace)
   DTrace => pDTrace%ptr

   call init_window(DTrace, taudiffmax)
end subroutine set_taudiffmax

!===============================================================================
subroutine ctqmc_warmup(Nwarmups)
!===============================================================================
   use type_progress

!local
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   type(progress)             :: p
   integer(c_int64_t)         :: i, Nwarmups
   real(KINDR)                :: rand_num
   real(KINDR)                :: time_start

   pDStates = transfer(ipDStates, pDStates)
   pDTrace = transfer(ipDTrace, pDTrace)
   DTrace => pDTrace%ptr

   if (allocated(AccPair)) then
      deallocate(AccPair)
      allocate(AccPair(NBands, 2, NGtau))
      AccPair = 0
      NAccPair = NGtau
      AccPairMax = DTrace%beta
   end if
   call init_counters()
   call cpu_time(time_start)

   p%adaptrate = 0.1
   call pstart(p, int(NWarmups,PINT), &
               title='Warmup ' // simstr // ':', fancy=fancy_prog)

   do i=1,NWarmups
      if(any_signal_fired) exit
      rand_num=grnd()
      if(rand_num<(1d0-PercentageTauShiftMove-PercentageGlobalMove)/2d0)then
         if(b_segment)then
            call StepAdd_mine()
         else
            if(grnd() >= Percentage4OperatorMove)then
               call StepAdd(grnd() < PercentageOuterMove,1)
            else
               if(.not.b_exch)then
                  call StepAdd4(1)
               else
                  call StepFlavourchange_general()
               endif
            endif
         endif
      elseif(rand_num<=1d0-PercentageTauShiftMove-PercentageGlobalMove)then
         if(b_segment)then
            call StepRem_mine()
         else
            if(grnd() >= Percentage4OperatorMove)then
               call StepRem(grnd() < PercentageOuterMove,1)
            else
               if(.not.b_exch)then
                  call StepRem4()
               else
                  call StepFlavourchange_general()
               endif
            endif
         endif
      elseif (rand_num <= 1d0-PercentageTauShiftMove) then
         if(b_segment)then
            call StepGlob_mine()
         else
            !write(*,*) "stepglob"
            call StepGlob()
         endif
      else
         call StepShiftTau()
      endif
      if (.not. b_segment) then
         if (modulo(sum(TryAdd) + sum(TryRem), NSlide) == 0) call shift_window(DTrace)
      end if
      call ptick(p)
   enddo

   call cpu_time(time_warmup)
   time_warmup = time_warmup - time_start
end subroutine ctqmc_warmup

!estimating eta, starting from a thermalized configuration in Z
!ctqmc_warmup should be called before
subroutine ctqmc_worm_warmup(Nwarmups2,iSector,iComponent)
   integer                    :: iSector,iComponent
!local
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   integer(c_int64_t)         :: Nwarmups2
   integer                    :: maxIter
   !sector: 1 == partition function, 2 == 1P GF worm
   !sector: 3 == IE Sigma worm,      4 == 2P GF worm 
   !sector: 5 == IE Chi worm
   !sector: 6 == 2 legged 2P GF PH
   !sector: 7 == 2 legged 2P GF PP
   !sector: 8 == 3 legged 2P GF PH
   !sector: 9 == 3 legged 2P GF PP
 
   pDStates = transfer(ipDStates, pDStates)
   pDTrace = transfer(ipDTrace, pDTrace)

   if((b_offdiag.eqv..true.).and.(iSector>2)) then
       stop "Worm Sampling with offdiagonal hybridizations not implemented."
   endif

   !for component sampling, just add one additional sector
   if(iComponent.ne.0) then
      !bufferd transform for g4iw and h4iw
      allocate(val_worm(10000000))
      tau_fill = 0

      if(iSector .eq. 2 .or. iSector .eq. 3) then
         allocate(tau_worm(size(val_worm)))
      elseif(iSector .eq. 4 .or. iSector .eq. 5) then
         allocate(tau_worm(3*size(val_worm)))
      elseif(iSector .eq. 6 .or. iSector .eq. 7) then
         allocate(tau_worm(size(val_worm)))
      elseif(iSector .eq. 8 .or. iSector .eq. 9) then
         allocate(tau_worm(2*size(val_worm)))
      endif

   endif
   
   !generalize this to components in sector
   if(get_Integer_Parameter("WormSearchEta")/=0) then
      !we determine the reweighing matrix eta in n_component preruns
      !using a root finding algorithm, currently with a 10% tolerance
      maxIter=15
      call findEta(iSector,iComponent,0.1d0,Nwarmups2,maxIter)
   endif
   

end subroutine ctqmc_worm_warmup

!ctqmc sampling with measurements
!Z-sampling/conventional CT-HYB: iSector=1,iComponent=1
!Worm-sampling/component samling: iSector set for worm
!estimator and iComponent set for flavor component
!===============================================================================
subroutine ctqmc_measure(iSector,iComponent)
!===============================================================================
   use LegendrePoly, only: m1pow
   use type_progress
   use signals

   integer                    :: iSector,iComponent
!local
   type(TTrace),pointer       :: DTrace
   type(TStates_pointer)      :: pDStates
   type(TTrace_pointer)       :: pDTrace
   type(progress)             :: p
   integer                    :: wm,l,lp
   integer(c_int64_t)         :: i,j,iNMeas,iNCorr, k
   integer                    :: ib, is
   integer                    :: in1, in2, in3
   real(KINDR)                :: rand_num
   real(KINDR)                :: time_start, time_temp
   real(KINDR)                :: wormScaling(9),zScaling
   integer                    :: Sector
   logical                    :: do_sampling
   real(KINDR)                :: wormNorm
   logical                    :: worm_offdiag
   integer                    :: bs(2), b(2), s(2)

   pDStates = transfer(ipDStates, pDStates)
   pDTrace = transfer(ipDTrace, pDTrace)
   DTrace => pDTrace%ptr

   if (allocated(AccPair)) then
      deallocate(AccPair)
      allocate(AccPair(NBands, 2, NGtau))
      AccPair = 0
      NAccPair = NGtau
      AccPairMax = DTrace%beta
   end if
   call init_counters()
   call pstart(p, int(NMeas,PINT), &
               title='Simulation ' // simstr, fancy=fancy_prog)

   call cpu_time(time_sim)

   Sector = 1

   !determine if off_diagonal one-particle worm
   worm_offdiag=.false.
   if(iSector .eq. 2) then
      call index2component_general(Nbands, 2, iComponent, bs, b, s)
      if(bs(1).ne.bs(2)) worm_offdiag=.true.
   endif

   !only sample if necessary
   do_sampling = .false.
   if (iSector.eq.1) then
      do_sampling = .true.
   else if (iSector > 1) then
      if (wormEta(iSector).ne.0d0) do_sampling = .true.
   end if

   if (do_sampling) then

   OUTER_LOOP: &
   do iNMeas=1,NMeas
      g_inmeas=inmeas
      !correlation steps in z and worm space
      !we measure sign with NCorr = 1

      call cpu_time(time_start)
      do iNCorr=1, NCorr
         !Z space
         if(any_signal_fired) exit OUTER_LOOP
         if(Sector==1) then
            ! This is set by trap_notify in case of a signal
            !increase sampling counter for partition function space
            CntSampling(Sector) = CntSampling(Sector) + 1
            rand_num=grnd()
            !insertion of hyb pair
            if(rand_num<(1d0-PercentageTauShiftMove-PercentageGlobalMove-PercentageWormInsert)/2d0)then
               if(b_segment)then
                  call StepAdd_mine()
               else
                  if(grnd() >= Percentage4OperatorMove)then
                     call StepAdd(grnd() < PercentageOuterMove,Sector)
                  else
                     if(.not.b_exch)then
                        call StepAdd4(Sector)
                     else
                        call StepFlavourchange_general()
                     endif
                  endif
               endif
            !removal of hyb pair
            elseif(rand_num<=1d0-PercentageTauShiftMove-PercentageGlobalMove-PercentageWormInsert)then
               if(b_segment)then
                  call StepRem_mine()
               else
                  if(grnd() >= Percentage4OperatorMove)then
                     call StepRem(grnd() < PercentageOuterMove,Sector)
                  else
                     if(.not.b_exch)then
                        call StepRem4()
                     else
                        call StepFlavourchange_general()
                     endif
                  endif
               endif
            !worm insertion
            elseif(rand_num<=1d0-PercentageTauShiftMove-PercentageGlobalMove)then

               Sector=iSector
               
               !offdiagonals require mixed worm/hyb for Sector 2
               if(worm_offdiag) then
                   call StepWormHybAdd(Sector,iComponent)
               else
                   call StepWormAdd(Sector,iComponent)
               endif
               
            !global move
            elseif (rand_num <= 1d0-PercentageTauShiftMove) then
               if(b_segment)then
                  call StepGlob_mine()
               else
                  !write(*,*) "stepglob"
                  call StepGlob()
               endif
            else
               call StepShiftTau()
            endif
         !worm spaces   
         else
            ! This is set by trap_notify in case of a signal
            !increase sampling counter for greens function space
            CntSampling(Sector) = CntSampling(Sector) + 1
            rand_num=grnd()
            !insertion of hyb pair
            if(rand_num<(1d0-PercentageWormReplace-PercentageWormInsert)/2d0)then
               if(grnd() >= Percentage4OperatorMove)then
                  call StepAdd(grnd() < PercentageOuterMove,Sector)
               else
                  if(.not.b_exch)then
                     call StepAdd4(Sector)
                  else
                     call StepFlavourchange_general()
                  endif
               endif
            !removal of of hyb pair
            elseif(rand_num<1d0-PercentageWormReplace-PercentageWormInsert)then
               if(grnd() >= Percentage4OperatorMove)then
                  call StepRem(grnd() < PercentageOuterMove,Sector)
               else
                  if(.not.b_exch)then
                     call StepRem4()
                  else
                     call StepFlavourchange_general()
                  endif
               endif
            elseif(rand_num<1d0-PercentageWormInsert) then
               if(Sector .ne. 6 .and. Sector .ne. 7) then
                  call StepWormReplace(Sector)
               endif
            !worm removal
            else
               !offdiagonals require mixed worm/hyb for Sector 2
               if(worm_offdiag) then
                   call StepWormHybRem(Sector)
               else
                   call StepWormRem(Sector)
               endif
            endif
         endif

         if (.not. b_segment) then
            if (modulo(sum(TryAdd) + sum(TryRem), NSlide) == 0) call shift_window(DTrace)
         end if
      enddo

      call cpu_time(time_temp)
      time_sim_steps = time_sim_steps + (time_temp - time_start)

      call ptick(p)
      
      !partition function space measurements
      if(Sector==1) then
         !increase measurement counter for partition function space
         CntMeas(Sector) = CntMeas(Sector) + 1

         call MeasSign()
         if(b_segment) call MeasOcc_seg()
         if (.not. b_segment) call MeasOuterHistograms()
         call MeasExpOrder()
         
         !only measure if component samling not enabled
         if(iSector==1) then
            if(b_meas_susz)then
               call MeasSusz_seg()
            else if (b_meas_susz_mat) then
               call MeasSusz()
            end if
            
            if(b_offdiag)then
               call MeasGtauRem_full()
            else
               call MeasGtauRem(GtauDetRat)
            endif

            if(allocated(giw) .or. allocated(gsigmaiw) .or. &
               allocated(g4iw) .or. allocated(g4iw_pp)) then
               call minv_matching_taus(DTrace,DTrace%MInv)
            endif
!           if(allocated(gsigmaiw)) call urho_contracted(u_matrix, dtrace, dstates)
            if(.not. b_Giw_lookup .and. (allocated(giw) .or. allocated(gsigmaiw))) then
               call measure_giw_nfft()
            endif
            if ((b_Eigenbasis .eqv. .true.) .and. .not. b_segment) then
               if(allocated(densitymatrix)) call MeasDensityMatrix_beta_half()
            endif
            if(allocated(expresdensitymatrix)) call MeasExpResDensityMatrix()

            if(IAND(FourPnt,GF4_IMAGTIME+GF4_LEGENDRE) /= 0) call MeasGtau4PntLeg()
            if(allocated(g4iw) .or. allocated(g4iw_pp)) call measure_g4iw_nfft()
            !call ptick(p)
         endif
       !worm space measurements
       elseif(Sector==2) then
            CntMeas(Sector) = CntMeas(Sector) + 1            
            if(allocated(Giw_worm)) call measure_giw_worm()
            if(allocated(Gtau_worm)) call measure_gtau_worm()
       elseif(Sector==4 .and. (allocated(g4iw_worm))) then
            CntMeas(Sector) = CntMeas(Sector) + 1
            call measure_g4iw_worm()
       elseif(Sector==3 .and. (allocated(gsigmaiw_worm))) then
            CntMeas(Sector) = CntMeas(Sector) + 1
            call measure_ie_sigma_worm()
       elseif(Sector==5 .and. (allocated(h4iw_worm))) then
            CntMeas(Sector) = CntMeas(Sector) + 1
            call measure_ie_chi_worm()
       elseif(Sector==6) then
            CntMeas(Sector) = CntMeas(Sector) + 1
            if (allocated(p2iw_worm)) call measure_p2iw_worm()
            if (allocated(p2tau_worm))  call measure_p2tau_worm()
       elseif(Sector==7) then
            CntMeas(Sector) = CntMeas(Sector) + 1
            if (allocated(p2iwpp_worm)) call measure_p2iwpp_worm()
            if (allocated(p2taupp_worm)) call measure_p2taupp_worm()
       elseif(Sector==8 .and. (allocated(p3iw_worm))) then
            CntMeas(Sector) = CntMeas(Sector) + 1
            call measure_p3iw_worm()
       elseif(Sector==9 .and. (allocated(p3iwpp_worm))) then
            CntMeas(Sector) = CntMeas(Sector) + 1
            call measure_p3iwpp_worm()
       endif
       
   enddo OUTER_LOOP

   call cpu_time(time_temp)
   time_sim = time_temp - time_sim

   if(any_signal_fired) then
      write(0,"('Rank',I3,'> WARNING: Process caught signal')") SimID
      write(0,"('Rank',I3,'; measurement',I10,' of',I10)")&
         SimID, iNMeas, NMeas
      NMeas = iNMeas-1
      if (Nmeas == 0) then
         aborted = 2   ! some data
      else
         aborted = 1
      endif
   endif

   !endif for if(do_sampling)
   endif

   if(NMeas == 0) then
      write (0,"('Rank', I3, '> ERROR: No measurements')") SimID
   else
      !CntMeas(1): measurements in partition function space
      !CntMeas(2): measurements in greens function space
      mean_sign = meas_sign/real(cnt_sign_z, KINDR)
      
      write(0,*) "Z Mean Sign", mean_sign
      if (iSector /= 1) then
         write(0,*) "Steps in Z and wormspaces", CntSampling(1),CntSampling(iSector)
      else
         write(0,*) "Steps in Z", CntSampling(1)
      endif
      !acceptance in partition function space

      AccGlob=AccGlob/TryGlob
      AccShift=AccShift/TryShift
      
      AccAdd(1)=AccAdd(1)/TryAdd(1)
      AccRem(1)=AccRem(1)/TryRem(1)
      if (iSector /= 1) then
         AccAdd(iSector)=AccAdd(iSector)/TryAdd(iSector)
         AccRem(iSector)=AccRem(iSector)/TryRem(iSector)
         AccWormAdd(iSector)=AccWormAdd(iSector)/TryWormAdd(iSector)
         AccWormRem(iSector)=AccWormRem(iSector)/TryWormRem(iSector)
         AccWormRep(iSector)=AccWormRep(iSector)/TryWormRep(iSector)
      end if

      !write(*,*) "total:"
      !write(*,*) "TryAdd4", TryAdd4
      !write(*,*) "TryRem4", TryRem4
      !write(*,*) "AccAdd4", AccAdd4
      !write(*,*) "AccRem4", AccRem4

      AccAdd4=AccAdd4/TryAdd4
      AccRem4=AccRem4/TryRem4
      !write(*,*) "relative:"
      !write(*,*) "AccAdd4", AccAdd4
      !write(*,*) "AccRem4", AccRem4

      !write(*,*) "TryFlavc ", TryFlavc 
      !write(*,*) "AccFlavc ", AccFlavc 
      AccFlavc=AccFlavc/dble(TryFlavc) 
      !write(*,*) "relative: ", AccFlavc

      zScaling = -1./(dble(CntMeas(1)) * mean_sign)
      lhisto(:) = dble(dtrace%lhisto(:))/sum(dtrace%lhisto(:)) !FIXME
      rhisto(:) = dble(dtrace%rhisto(:))/(AccAdd(1)+AccRem(1)+AccGlob)
      if(.not.b_segment)then
         occ(:,:,:,:)=occ(:,:,:,:) * zScaling
      else
         occ(:,:,:,:)=occ(:,:,:,:) * zScaling / dtrace%beta
      endif

      if(allocated(densitymatrix)) then
         densitymatrix(:,:)=densitymatrix(:,:) * zScaling
      endif

      if ((b_Eigenbasis .eqv. .true.) .and. (b_segment .eqv. .false.)) then
         call Meas_rho1_from_Densitymatrix()
         call Meas_rho2_from_Densitymatrix()
         call MeasSingleOcc_from_Densitymatrix()
         call MeasDoubleOcc_from_Densitymatrix()
         call Transform_DensityMatrix_EB_to_OCCB()
      endif

      if(allocated(expresdensitymatrix))then
        forall(ib = 1:NBands, is = 1:2, i=0:(size(Histo(1,1,:))-1), Histo(ib,is,i)/=0)
          expresdensitymatrix(ib,is,i,:,:)=expresdensitymatrix(ib,is,i,:,:)/Histo(ib,is,i) ! FIXME
        endforall
        write (0, *) "WARNING: unimplemented: the measured expansion order resolved density matrix MUST be transformed from eigenbasis to occupation number basis to get correct results!"
      endif
      Gtau(:,:,:)=-Gtau(:,:,:) * zScaling &
         *dble(NGtau)/(DTrace%beta)**2
      if (allocated(Gtau_mean_step)) Gtau_mean_step(:,:,:)=Gtau_mean_step(:,:,:)/(DTrace%beta)**2
      if (allocated(Gtau_mid_step)) Gtau_mid_step(:,:,:)=Gtau_mid_step(:,:,:)/(DTrace%beta)**2

      if (b_meas_susz_mat) then
         ntau_n0 = ntau_n0 * zScaling
         ntau_n0(:,:,:,:,NGtau)=ntau_n0(:,:,:,:,0)
      else
         ntau_n0(:,:,:,:,:)=ntau_n0(:,:,:,:,:) * zscaling / float(ngtau)
         ntau_n0(:,:,:,:,Ngtau)=ntau_n0(:,:,:,:,0)
      end if

      Gtau_full(:,:,:,:,:)=-Gtau_full(:,:,:,:,:) * zScaling &
         *dble(NGtau)/(DTrace%beta)**2

      Histo(:,:,:) = Histo(:,:,:)/dble(CntMeas(1))

      SignHistoSuperstates(:) = SignHistoSuperstates(:)/OuterSuperstateHisto(:)
      SignHistoStates(:) = SignHistoStates(:)/OuterStateHisto(:)

      OuterSuperstateHisto(:) = OuterSuperstateHisto(:)/dble(CntMeas(1))
      OuterStateHisto(:) = OuterStateHisto(:)/dble(CntMeas(1))

      TraceContribSuperstates(:) = TraceContribSuperstates(:)/dble(CntMeas(1))
      TraceContribStates(:) = TraceContribStates(:)/dble(CntMeas(1))

      do i=1,NLegMax
         GLeg(:,:,i)=GLeg(:,:,i)*sqrt(dble(2*(i-1))+1d0)/DTrace%beta
         GLeg_full(:,:,:,:,i)=GLeg_full(:,:,:,:,i)*sqrt(dble(2*(i-1))+1d0)/DTrace%beta
      enddo
      GLeg(:,:,:)=-GLeg(:,:,:) * zScaling 
      GLeg_full(:,:,:,:,:)=-GLeg_full(:,:,:,:,:) * zScaling 

      if(SimID.eq.0)then
         write(*,*) "globalmove_check ", globalmove_check 
         globalmove_check(:)=globalmove_check(:)/(sum(globalmove_check))
         write(*,*) "globalmove_check ", globalmove_check 
      endif

      ! Add the overall normalisation factors (Boehnke 2011, C24) and statistical
      ! weight 1/N to the 4-point measurement
      if( IAND(FourPnt,GF4_LEGENDRE) /= 0 ) then
        forall(l = 0:N4leg-1, lp = 0:N4leg-1, wm = -N4iwb:N4iwb)
          G4leg(:,:,:,:,l,lp,wm) = G4leg(:,:,:,:,l,lp,wm) * zScaling * &
            sqrt(2*l+1._KINDR)*sqrt(2*lp+1._KINDR)*m1pow(lp+1)/DTrace%beta
        end forall
      endif

      ! FIX: scaling of the imaginary-time GF4 as well
      if( IAND(FourPnt,GF4_IMAGTIME) /= 0 ) then
        G4tau(:,:,:,:,:,:,:) = G4tau(:,:,:,:,:,:,:) * zScaling * &
          (dble(N4tau)/(DTrace%beta))**3 /DTrace%beta
      endif

      ! Rescale the
      if( allocated(g4iw) ) then
        G4iw(:,:,:,:,:,:,:) = G4iw(:,:,:,:,:,:,:) * zScaling/DTrace%beta**2
      endif

      if( allocated(g4iw_pp) ) then
        G4iw_pp(:,:,:,:,:,:,:) = G4iw_pp(:,:,:,:,:,:,:) * zScaling/DTrace%beta**2
      endif

      if( allocated(g2iw) ) then
        G2iw(:,:,:,:) = G2iw(:,:,:,:) * zScaling/DTrace%beta
      endif

      ! timings
      if(allocated(timings_giw)) timings_giw(:) = timings_giw(:)/count_giw(:)
      if(allocated(timings_g4iw_add)) timings_g4iw_add(:) = timings_g4iw_add(:)/count_g4iw(:)
      if(allocated(timings_g4iw_ft)) timings_g4iw_ft(:) = timings_g4iw_ft(:)/count_g4iw(:)

#ifdef USE_NFFT
      if (b_Giw_lookup) then
         do ib = 1, NBands
            do is = 1, 2
               if (Giw_lookup_row(ib, is) > 1) call MeasGiwFlushLookup(ib, is)
            end do
         end do
      end if
#endif
      if( allocated(giw) ) then
        Giw(:,:,:) = Giw(:,:,:) * zScaling/DTrace%beta
      endif

      if( allocated(gsigmaiw) ) then
        gsigmaiw(:,:,:) = gsigmaiw(:,:,:) * zScaling/DTrace%beta
      endif
      
      if(iSector==2) wormNorm=4d0/DTrace%beta !(2^2 for spin)
      if(iSector==3) wormNorm=32d0/(DTrace%beta) !(2^5 for spin)
      if(iSector==4) wormNorm=16d0/(DTrace%beta**2) !(2^4 for spin)
      !FIXME: check beta factor for H
      if(iSector==5) wormNorm=128d0/(DTrace%beta**2) !(2^7 for spin)
      if(iSector==6) wormNorm=16d0/(DTrace%beta) !(2^4 for spin)
      if(iSector==7) wormNorm=16d0/(DTrace%beta) !(2^4 for spin)
      if(iSector==8) wormNorm=16d0/(DTrace%beta**2) !(2^4 for spin)
      if(iSector==9) wormNorm=16d0/(DTrace%beta**2) !(2^4 for spin)    

      if(CntMeas(iSector).ne.0) then
         wormScaling(iSector)=dble(CntSampling(iSector))/(dble(CntMeas(iSector))*dble(CntSampling(1)))&
         *wormNorm/mean_sign
      else 
         !capturing 0/0 events
         wormScaling(iSector)=0d0
      endif
      
      
      ! worm normalization
      if( allocated(Giw_worm) .or. allocated(Gtau_worm)) then

         if(allocated(Giw_worm)) then
            if(tau_fill .ne. 0) then
               allocate(matsubaras(4*NGiw),d(1))
               d(1)=size(matsubaras)
             
               call ft_nd(tau_worm(:tau_fill),val_worm(:tau_fill),matsubaras,d)
               giw_worm(:)=giw_worm(:)+matsubaras(2::2)          

               deallocate(matsubaras,d)
             endif
             
             Giw_worm(:) = -Giw_worm(:)*wormScaling(2)
         endif

         if(allocated(Gtau_worm)) then
            Gtau_worm(:) = Gtau_worm(:)*wormScaling(2)*(dble(NGtau)/Dtrace%beta)
         endif

      endif

      ! vertex normalization
      if( allocated(g4iw_worm) ) then
      
        !finalize the fourier transform for what is left in wormSum
         if(tau_fill .ne. 0) then
            allocate(matsubaras((4*N4iwf)*(4*N4iwf)*(4*N4iwb+3)),d(3))
            d(1)=4*N4iwf
            d(2)=4*N4iwf
            d(3)=4*N4iwb+3

            call ft_nd(tau_worm(:3*tau_fill),val_worm(:tau_fill),matsubaras,d)

            do in1=-N4iwf,N4iwf-1
               do in2=-N4iwf,N4iwf-1
                  do in3=-N4iwb,N4iwb

                   k=(4*N4iwf)*(4*N4iwb+3)*(2*(in1+N4iwf)+1)+(4*N4iwb+3)*(2*(in2+N4iwf)+1)+(2*(in3+N4iwb)+1)+1
                   g4iw_worm(in1,in2,in3)=g4iw_worm(in1,in2,in3) + matsubaras(k)

                  enddo
               enddo
            enddo

            deallocate(matsubaras,d)         
         endif

         g4iw_worm(:,:,:)=g4iw_worm(:,:,:)*wormScaling(4)
      endif
      
      ! IE-SIGMA normalization
      ! FIXME: we dont do this for the full estimator
      if( allocated(gsigmaiw_worm) ) then
         if(tau_fill .ne. 0) then
            allocate(matsubaras(4*NGiw),d(1))
            d(1)=size(matsubaras)
 
            call ft_nd(tau_worm(:tau_fill),val_worm(:tau_fill),matsubaras,d)
            gsigmaiw_worm(:)=gsigmaiw_worm(:)+matsubaras(2::2)         
 
            deallocate(matsubaras,d)
         endif

         Gsigmaiw_worm(:) = Gsigmaiw_worm(:)*wormScaling(3)
      endif
      
      ! IE-CHI normalization
      if( allocated(h4iw_worm) ) then

         !finalize the fourier transform for what is left in wormSum
         if(tau_fill .ne. 0) then
            allocate(matsubaras((4*N4iwf)*(4*N4iwf)*(4*N4iwb+3)),d(3))
            d(1)=4*N4iwf
            d(2)=4*N4iwf
            d(3)=4*N4iwb+3

            call ft_nd(tau_worm(:3*tau_fill),val_worm(:tau_fill),matsubaras,d)

            do in1=-N4iwf,N4iwf-1
               do in2=-N4iwf,N4iwf-1
                  do in3=-N4iwb,N4iwb

                      k=(4*N4iwf)*(4*N4iwb+3)*(2*(in1+N4iwf)+1)+(4*N4iwb+3)*(2*(in2+N4iwf)+1)+(2*(in3+N4iwb)+1)+1
                      h4iw_worm(in1,in2,in3)=h4iw_worm(in1,in2,in3) + matsubaras(k)

                  enddo
               enddo
            enddo

            deallocate(matsubaras,d)
         endif
        
         h4iw_worm(:,:,:)=h4iw_worm(:,:,:)*wormScaling(5)
      endif
          
      ! p2iw normalization
      if( allocated(p2iw_worm) ) then

        !finalize the fourier transform for what is left in wormSum
        if(tau_fill .ne. 0) then
           allocate(matsubaras(4*N2iwb+3),d(1))
           d(1)=size(matsubaras)
 
           call ft_nd(tau_worm(:tau_fill),val_worm(:tau_fill),matsubaras,d)
           p2iw_worm(:)=p2iw_worm(:)+matsubaras(2::2)

           deallocate(matsubaras,d)
        endif
   
        p2iw_worm(:)=p2iw_worm(:)*wormScaling(6)
      endif
      
      ! p2iwpp normalization
      if( allocated(p2iwpp_worm) ) then
      
        !finalize the fourier transform for what is left in wormSum
        if(tau_fill .ne. 0) then
           allocate(matsubaras(4*N2iwb+3),d(1))
           d(1)=size(matsubaras)

           call ft_nd(tau_worm(:tau_fill),val_worm(:tau_fill),matsubaras,d)
           p2iwpp_worm(:)=p2iwpp_worm(:)+matsubaras(2::2)

           deallocate(matsubaras,d)
        endif

        p2iwpp_worm(:)=p2iwpp_worm(:)*wormScaling(7)
      endif

      if( allocated(p2tau_worm)) then
        p2tau_worm(:)=p2tau_worm(:)*wormScaling(6)*(dble(NGtau)/Dtrace%beta)
      endif

      if( allocated(p2taupp_worm)) then
        p2taupp_worm(:)=p2taupp_worm(:)*wormScaling(7)*(dble(NGtau)/Dtrace%beta)
      endif
      
      ! p3iw normalization
      ! FIXME: we dont do this for the full estimator
      if( allocated(p3iw_worm) ) then
        
        !finalize the fourier transform for what is left in wormSum
        if(tau_fill .ne. 0) then
           allocate(matsubaras((4*N3iwf)*(4*N3iwb+3)),d(2))
           d(1)=4*N3iwf
           d(2)=4*N3iwb+3

           call ft_nd(tau_worm(:2*tau_fill),val_worm(:tau_fill),matsubaras,d)

           do i=-N3iwf,N3iwf-1
              do j=-N3iwb,N3iwb
                  k=(4*N3iwb+3)*(2*(i+N3iwf)+1)+(2*(j+N3iwb)+1)+1
                  p3iw_worm(i,j)=p3iw_worm(i,j)+matsubaras(k)
              enddo
           enddo

           deallocate(matsubaras,d)
         endif
      
         p3iw_worm(:,:)=p3iw_worm(:,:)*wormScaling(8)
      endif


      ! p3iwpp normalization
      ! FIXME: we dont do this for the full estimator
      if( allocated(p3iwpp_worm) ) then
        
        !finalize the fourier transform for what is left in wormSum
        if(tau_fill .ne. 0) then
           allocate(matsubaras((4*N3iwf)*(4*N3iwb+3)),d(2))
           d(1)=4*N3iwf
           d(2)=4*N3iwb+3

           call ft_nd(tau_worm(:2*tau_fill),val_worm(:tau_fill),matsubaras,d)

           do i=-N3iwf,N3iwf-1
              do j=-N3iwb,N3iwb
                  k=(4*N3iwb+3)*(2*(i+N3iwf)+1)+(2*(j+N3iwb)+1)+1
                  p3iwpp_worm(i,j)=p3iwpp_worm(i,j)+matsubaras(k)
              enddo
           enddo

           deallocate(matsubaras,d)
        endif

         p3iwpp_worm(:,:)=p3iwpp_worm(:,:)*wormScaling(9)
      endif

      write (0,"('Done post-processing, CTQMC ', A4, 'complete.')") simstr
   endif   
   
   call reset_notify()
   call clear_trap(SIGNAL_TERM)
   call clear_trap(SIGNAL_INT)
   call clear_trap(SIGNAL_USR1)
end subroutine ctqmc_measure

subroutine get_mc_config_scalars(N, outer_sst, outer_state)
   integer, intent(out)                                :: N, outer_sst, outer_state
   type(TTrace), pointer                               :: DTrace
   type(TStates_pointer)                               :: pDStates
   type(TTrace_pointer)                                :: pDTrace

   pDStates = transfer(ipDStates, pDStates)
   pDTrace = transfer(ipDTrace, pDTrace)
   DTrace => pDTrace%ptr

   N = DTrace%NOper
   outer_sst = DTrace%outer_sst
   outer_state = DTrace%outer_state
end subroutine get_mc_config_scalars

subroutine get_mc_config(N, taus, orbs, spins, cas, hashybs, outer_sst, outer_state)
   integer, intent(in)                                 :: N 
   real(KINDR), dimension(N), intent(out)              :: taus
   integer, dimension(N), intent(out)                  :: orbs, spins, cas, hashybs
   !f2py depend(N) taus, orbs, spins, cas, hashybs
   integer, intent(out)                                :: outer_sst, outer_state
   type(TTrace), pointer                               :: DTrace
   type(TTrace_pointer)                                :: pDTrace

   pDTrace = transfer(ipDTrace, pDTrace)
   DTrace => pDTrace%ptr

   call dump_mc_config_into_arrays(DTrace, taus, orbs, spins, cas, hashybs, outer_sst, outer_state)
end subroutine get_mc_config

subroutine set_mc_config(N, taus, orbs, spins, cas, hashybs, outer_sst, outer_state)
   !input
   integer                                :: N
   real(KINDR), dimension(N), intent(in)  :: taus
   integer, dimension(N), intent(in)      :: orbs, spins, cas, hashybs
   !f2py depend(N) taus, orbs, spins, cas, hashybs
   integer, intent(in)                    :: outer_sst, outer_state
   type(TStates), pointer                 :: DStates
   type(TTrace), pointer                  :: DTrace
   type(TStates_pointer)                  :: pDStates
   type(TTrace_pointer)                   :: pDTrace

   pDStates = transfer(ipDStates, pDStates)
   pDTrace = transfer(ipDTrace, pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   call set_mc_config_from_arrays(DTrace, DStates, taus, orbs, spins, cas,&
                                  hashybs, outer_sst, outer_state)
end subroutine set_mc_config

subroutine set_empty_mc_config(outer_sst, outer_state)
   !input
   integer, intent(in)                    :: outer_sst, outer_state
   type(TStates), pointer                 :: DStates
   type(TTrace), pointer                  :: DTrace
   type(TStates_pointer)                  :: pDStates
   type(TTrace_pointer)                   :: pDTrace
   real(KINDR), dimension(0)              :: taus
   integer, dimension(0)                  :: orbs, spins, cas, hashybs

   pDStates = transfer(ipDStates, pDStates)
   pDTrace = transfer(ipDTrace, pDTrace)
   DStates => pDStates%ptr
   DTrace => pDTrace%ptr

   call set_mc_config_from_arrays(DTrace, DStates, taus, orbs, spins, cas,&
                                  hashybs, outer_sst, outer_state)
end subroutine set_empty_mc_config

!>\brief initialise the random number generator before starting the qmc for the
!! first time
subroutine init_rnd(Nseed)
!input
  integer,intent(in)    :: Nseed
  call sgrnd(Nseed)
end subroutine init_rnd

!>\brief initialise the parameters with values defined in parastr
!! the format is the same as for the parameter file
subroutine init_paras(parastr)
!input
   character(*),intent(in) :: parastr
   call read_ParameterString(parastr)
end subroutine init_paras

!>\brief initialise the parameters with values defined in parastr
!! the format is the same as for the parameter file
subroutine dest_paras()
!input
   call dest_Parameters()
end subroutine dest_paras

!===============================================================================
end module MCTQMC
!===============================================================================

#ifdef CTQMC_Test

!===============================================================================
program Prog_CTQMC
!===============================================================================
use MParameters
use MCTQMC

integer,parameter          :: Neps=5000,iwMax=2000
!integer,parameter          :: Neps=1000,iwMax=10000
real(KINDR),parameter      :: piii=dacos(-1d0)
integer                    :: i,j,k,iB,iS,iB1,iB2,iS1,iS2,iw,ieps,NBin,indx
integer                    :: idt1,idt2,idt3,idt12p,idt12m,idt32p,idt32m
real(KINDR)                :: tau,tau1,tau2,tau3,mu,epsmax,deps,eps
real(KINDR),allocatable    :: FTau(:,:,:),G0tau(:),GGtau(:,:,:,:,:,:,:),muimp(:,:)
complex(KINDR)             :: G0iw(-iwMax:iwMax-1),Fiw(-iwMax:iwMax-1)
integer                    :: Nseed,Ntau
real(KINDR), allocatable :: u_matrix(:,:,:,:)

!beware if detrat is smaller than d-12 we might get numerical problems so
!we should check for that
call read_ParameterFile("Parameters.in")
NTau=get_Integer_Parameter("Ntau")
NSeed=get_Integer_Parameter("NSeed")
NBands=get_Integer_Parameter("Nd")
mu=get_Real_Parameter("mu")
beta=get_Real_Parameter("beta")
epsmax=get_Real_Parameter("half-bandwidth")
NBin=get_Integer_Parameter("NBin")
FourPnt=get_Integer_Parameter("FourPnt")

call init_rnd(Nseed)
allocate(Ftau(NBands,2,0:NTau-1))
allocate(muimp(NBands,2))
allocate(G0tau(0:NTau-1))

call u_allocate(NBands,u_matrix)
u_matrix=0d0

muimp(:,:)=0d0
Ftau=0d0
do iB=1,NBands
G0tau=0d0
G0iw=cmplx(0d0, 0d0, KINDR)
Fiw=cmplx(0d0, 0d0, KINDR)
deps=(2d0*epsmax/dble(Neps))
do iw=-iwmax,iwmax-1
   do ieps=1,Neps
      eps=-dble(epsmax)+deps*dble(ieps)
! semi-circular density
      Fiw(iw)=Fiw(iw) + deps*2d0/(piii*epsmax)*sqrt(1d0-(eps/epsmax)**2)/&
         cmplx(-eps,dble(2*iw+1)*piii/beta, KINDR)
! flat bath
!      Fiw(iw)=Fiw(iw)+deps/(2d0*epsmax)/cmplx(-eps,dble(2*iw+1)*piii/beta, KINDR)
   enddo
   G0iw(iw)=1d0/(cmplx(mu, piii*dble(2*iw+1)/beta, KINDR)-Fiw(iw))
!   write(1005,*)mu,beta,deps,epsmax,eps
!   write(1002,*)piii*dble(2*iw+1)/beta,real(G0iw(iw)),aimag(G0iw(iw))
!   write(1005,*)piii*dble(2*iw+1)/beta,real(Fiw(iw)),aimag(Fiw(iw))
enddo

do i=0,Ntau-1
   tau=beta/dble(Ntau-1)*dble(i)
   if(tau.eq.0d0)tau=0d0+1d-2
   if(tau.eq.beta)tau=beta-1d-2
   do iw=-iwmax,iwmax-1
      G0tau(i)=G0tau(i)+real(G0iw(iw)*exp(-cmplx(0d0,piii*dble(2*iw+1)/beta*tau, KINDR)))
      Ftau(iB,1,i)=Ftau(iB,1,i)+real(Fiw(iw)*exp(-cmplx(0d0,piii*dble(2*iw+1)/beta*tau, KINDR)))
   enddo
   G0tau(i)=G0tau(i)/beta
   Ftau(iB,:,i)=Ftau(iB,1,i)/beta
!   write(1090,*)tau,G0tau(i)
!   write(1003,*)tau,Ftau(iB,1,i)
enddo
!epsmax=epsmax*4d0
enddo
write(*,*)"*********************************************************************"
write(*,*)"                          CTQMC-Simulation"

!TODO: this makes problems now!!!
STOP "Please Fix U-Matrix"
call CTQMCSimulation(u_matrix,Ftau,muimp,NTau,NBands,0)

write(30,*)"Acceptance Meas:"
write(30,*)"Add: ",dble(AccAdd)/dble(TryAdd)
write(30,*)"total Add: ",dble(AccAdd)
write(30,*)"Rem: ",dble(AccRem)/dble(TryRem)
write(30,*)"total Rem: ",dble(AccRem)
write(30,*)"Glob: ",dble(AccGlob)/dble(TryGlob)
write(30,*)"total Glob: ",dble(AccGlob)
write(30,*)"AddQN: ",dble(AccQNAdd)/dble(TryAdd)
write(30,*)"RemQN: ",dble(AccQNRem)/dble(TryRem)
write(30,*)"GlobQN: ",dble(AccQNGlob)/dble(TryGlob)
write(30,*)"Mean Sign: ",dble(meas_sign)/dble(AccAdd+AccRem+AccGlob)
write(30,*)"<n n>"

do iB1=1,NBands
do iS1=1,2
   do iB2=1,NBands
   do iS2=1,2
      write(30,'(f10.7)',advance="no")occ(iB1,iS1,iB2,iS2)
   enddo
   enddo
   write(30,*)
enddo
enddo

do iB=1,NBands
do iS=1,2
do i=0,size(Histo(iB,iS,:))-1
   write(20,*)dble(Histo(iB,iS,i))
enddo
enddo
enddo

! Here we write the 1 particle GF using the tau-grid given by beta and NGtau.
do iB=1,NBands
do iS=1,2
do i=1,size(Gtau(iB,iS,:))
   tau=beta/dble(Ntau)*dble(i)-beta/(2d0*dble(Ntau))
!   write(100+iB*10+iS,*)tau,Gtau(iB,iS,i)
enddo
enddo
enddo

! Here we write the linearized 1 particle GF using the function get_GLin.
! Arbitrary tau´s are here allowed. The points off the tau-grid given by
! beta and NGtau are calculated by linearization.
! Additionally we mirrow the GF about beta/2 to meet the Toschi convention.
! This is done by writing G(beta-tau) instead of G(tau).
do iB=1,NBands
do iS=1,2
do i=1,2*size(Gtau(iB,iS,:))-1
   tau=beta/dble(2*NTau)*(dble(i))
!   write(900+iB*10+iS,*)tau,get_GLin(Gtau(iB,iS,:),NTau,beta,beta-tau)
!  write(9900+iB*10+iS,*)tau,get_GLin(Gtau(iB,iS,:,1),NTau,beta,tau)
enddo
enddo
enddo

if(FourPnt.ne.0)then
do iB1=1,NBands
do iB2=1,NBands
do iS1=1,2
do iS2=1,2
!  The maximum rank for Fortran arrays is 7, so here we use for (Spin1,Spin2)=(i,j) the
!  spin index convention s=2*int(s1/2)+s2: (1,1)=1, (1,2)=2, (2,1)=3, (2,2)=4 
   iS=2*int(iS1/2)+iS2
   indx=0
do i=1,NTau
   tau1=beta/dble(Ntau)*(dble(i)-0.5)
   idt1=int(tau1/beta*dble(Ntau))+1
do j=1,NTau
   tau2=beta/dble(Ntau)*(dble(j)-0.5)
   idt2=int(tau2/beta*dble(Ntau))+1
   write(40000+iB1*1000+iB2*100+iS1*10+iS2,'(A,I6)')'index=',indx
do k=1,NTau
   tau3=beta/dble(Ntau)*(dble(k)-0.5)
   idt3=int(tau3/beta*dble(Ntau))+1
   write(40000+iB1*1000+iB2*100+iS1*10+iS2,'(3(F8.3),2(E15.6))')&
      tau1,tau2,tau3,GGtau(iB1,iB2,iS,idt1,idt2,idt3,1),G4tau(iB1,iB2,iS,idt1,idt2,idt3,1)
enddo
   indx=indx+1
   write(40000+iB1*1000+iB2*100+iS1*10+iS2,*)
   write(40000+iB1*1000+iB2*100+iS1*10+iS2,*)
enddo
enddo
enddo
enddo
enddo
enddo
endif

! Here we write the 1 particle GF meassured in Legendre coefficients.
!do iB=1,NBands
!do iS=1,2
!do i=1,size(GLeg(iB,iS,:))
!!   write(300+iB*10+iS,*)i,GLeg(iB,iS,i)
!enddo
!enddo
!enddo

write(*,*)"**********************************************************************"
deallocate(Ftau)
deallocate(G0tau)
call dest_CTQMC()
call dest_Parameters()
end program Prog_CTQMC
#endif

!all the stuff from Interface.F90 still needed,
!shouldnt be here
!===============================================================================
subroutine invFT(Funtau,Funiw,w,tau,beta,NBands,Niw,Ntau,NNeq)
!===============================================================================
! input
   integer,intent(in)      :: Ntau,Niw,NBands,NNEq
   complex(KIND=KIND(0.D0)),intent(in)   :: Funiw(NNeq,NBands,2,-Niw:Niw-1)
   real(KIND=KIND(0.D0)),intent(in)       :: beta
   real(KIND=KIND(0.D0)),intent(in)       :: w(-Niw:Niw-1)
   real(KIND=KIND(0.D0)),intent(in)       :: tau(0:Ntau-1)
! output
   real(KIND=KIND(0.D0)),intent(out)      :: Funtau(NNeq,NBands,2,0:Ntau-1)
! local
   complex(KIND=KIND(0.D0))              :: Phase
   real(KIND=KIND(0.D0))                  :: taur
   integer                 :: itau,iw,ineq

   Funtau=0d0
   do ineq=1,NNeq
      do itau=0,Ntau-1
         taur=tau(itau)
         do iw=-Niw,Niw-1
            Phase=exp(-cmplx(0d0,w(iw)*taur, KIND = KIND(0.D0)))
            Funtau(ineq,:,:,itau)=Funtau(ineq,:,:,itau)+dble(Funiw(ineq,:,:,iw)*Phase)
         enddo
      enddo
   enddo
   Funtau=Funtau/beta
end subroutine invFT

!===============================================================================
subroutine FT(Funiw,Funtau,w,tau,beta,NBands,Ntau,Niw,NNeq)
!===============================================================================
! input
   integer,intent(in)                     :: Ntau,Niw,NBands,NNeq
   real(KIND=KIND(0.D0)),intent(in)       :: Funtau(NNeq,NBands,2,0:Ntau-1)
   real(KIND=KIND(0.D0)),intent(in)       :: beta
   real(KIND=KIND(0.D0)),intent(in)       :: w(-Niw:Niw-1)
   real(KIND=KIND(0.D0)),intent(in)       :: tau(0:Ntau-1)
! output
   complex(KIND=KIND(0.D0)),intent(out)   :: Funiw(NNeq,NBands,2,-Niw:Niw-1)
! local
   complex(KIND=KIND(0.D0))               :: Phase
   real(KIND=KIND(0.D0))                  :: taur
   integer                 :: itau,iw,ineq

   Funiw=cmplx(0d0, 0d0, KIND=KIND(0.D0))
   do ineq=1,NNEq
      do itau=0,Ntau-1
         taur=tau(itau)
         if(taur.eq.0d0)taur=0d0+3d-2
         if(taur.eq.beta)taur=beta-3d-2
         do iw=-Niw,Niw-1
            Phase=exp(cmplx(0d0,w(iw)*taur, KIND=KIND(0.D0)))
            Funiw(ineq,:,:,iw)=Funiw(ineq,:,:,iw)+Funtau(ineq,:,:,itau)*Phase
         enddo
      enddo
   enddo
   Funiw=Funiw*beta/dble(Ntau-1)
end subroutine FT

!===============================================================================
subroutine FlatHybr(w,mu,Fiw,G0iwinv,epsmax,Neps,NBands,Niw)
!===============================================================================
!input
   integer,intent(in)       :: NBands,Neps,Niw
   real(KIND=KIND(0.D0)),intent(in)        :: w(-Niw:Niw-1),mu
   real(KIND=KIND(0.D0)),intent(in)        :: epsmax(NBands)
!output  
   complex(KIND=KIND(0.D0)),intent(out)   :: Fiw(1,NBands,2,-Niw:Niw-1)
   complex(KIND=KIND(0.D0)),intent(out)   :: G0iwinv(1,NBands,2,-Niw:Niw-1)
!local
   integer                  :: iB,iS,iw,ieps
   real(KIND=KIND(0.D0))                   :: eps,deps

   Fiw=cmplx(0d0, 0d0, KIND=KIND(0.D0))
   do iB=1,NBands
      do iS=1,2 
         deps=(2d0*epsmax(iB)/dble(Neps-1))
         do iw=-Niw,Niw-1
            do ieps=0,Neps-1
               eps=-dble(epsmax(iB))+deps*dble(ieps)
               Fiw(1,iB,iS,-iw-1)=Fiw(1,iB,iS,-iw-1) + deps/(2d0*epsmax(iB))/cmplx(-eps, w(iw), KIND=KIND(0.D0))
            enddo
            G0iwinv(1,iB,iS,iw)=cmplx(mu, w(iw), KIND=KIND(0.D0))-Fiw(1,iB,iS,-iw-1)
         enddo
      enddo 
   enddo
end subroutine FlatHybr

!===============================================================================
subroutine SemiCircHybr(w,mu,Fiw,G0iwinv,epsmax,Neps,NBands,Niw)
!===============================================================================
!input
   integer,intent(in)      :: NBands,Neps,Niw
   real(kind=kind(0.D0)),intent(in)       :: w(-Niw:Niw-1),mu
   real(kind=kind(0.D0)),intent(in)       :: epsmax(NBands)
!output
   complex(kind=kind(0.D0)),intent(out)   :: Fiw(1,NBands,2,-Niw:Niw-1)
   complex(kind=kind(0.D0)),intent(out)   :: G0iwinv(1,NBands,2,-Niw:Niw-1)
!local
   integer                                :: iB,iS,iw,ieps
   real(kind=kind(0.D0))                  :: pi=dacos(-1d0),eps,deps

   Fiw = cmplx(0d0,0d0, kind(0.D0))
   do iB=1,NBands
      do iS=1,2
         deps=(2d0*epsmax(iB)/dble(Neps-1))
         do iw=-Niw,Niw-1
            do ieps=0,Neps-1
               eps=-dble(epsmax(iB))+deps*dble(ieps)
               Fiw(1,iB,iS,-iw-1)=Fiw(1,iB,iS,-iw-1)+deps*2d0/(pi*epsmax(iB))*sqrt(1d0-(eps/epsmax(iB))**2)/&
                             cmplx(-eps, w(iw), kind(0.D0))
            enddo
            G0iwinv(1,iB,iS,iw)=cmplx(mu, w(iw), KIND=KIND(0.D0))-Fiw(1,iB,iS,-iw-1)
         enddo
      enddo 
   enddo 
end subroutine SemiCircHybr

!===============================================================================
subroutine get_gbethe_para(glattiw,mu,sfw,w,d,Hmean,Ndp,Nneq)
!===============================================================================
!input
   integer,intent(in)                        :: Ndp,Nneq
   real(kind=kind(0.D0)),intent(in)          :: w,mu
   complex(kind=kind(0.D0)),intent(in)       :: sfw(Nneq,Ndp)
   real(kind=kind(0.D0)),intent(in)          :: d(Ndp)
   real(kind=kind(0.D0)),intent(in)          :: Hmean(NNEq,Ndp)
!output
   complex(kind=kind(0.D0)),intent(out)      :: glattiw(Nneq,Ndp)
!local
   integer                                   :: ib,ineq
   complex(kind=kind(0.D0))                  :: zeta

  do ib=1,Ndp
      do ineq=1,Nneq
         zeta=cmplx(mu, w, kind=kind(0.D0))-Hmean(ineq,ib)-sfw(ineq,ib)
         glattiw(ineq,ib)=2d0*zeta/d(ib)**2 * (1d0-sqrt(1d0-d(ib)**2/zeta**2))
      enddo 
   enddo
end subroutine get_gbethe_para

!===============================================================================
subroutine get_gsumhk_para(glattiw,mu,sfw,w,hk,neq2at,at2neq,dc,Ndp,Nneq,Nk,Nat,Ntb)
!===============================================================================
!input
   integer,intent(in)                       :: Ndp,Nneq,Nk,Nat,Ntb
   real(kind=kind(0.D0)),intent(in)         :: w,mu
   complex(kind=kind(0.D0)),intent(in)      :: sfw(Nneq,Ndp)
   complex(kind=kind(0.D0)),intent(in)      :: hk(Nk,Ntb,Ntb)
   integer,intent(in)                       :: neq2at(Nneq),at2neq(Nat)
   real(kind=kind(0.D0)),intent(in)         :: dc(Nneq,Ndp)
!output
   complex(kind=kind(0.D0)),intent(out)     :: glattiw(Nneq,Ndp)
!local
   integer                                  :: ib,ik,im,ineq,iat
   complex(kind=kind(0.D0))                 :: gfull(Ntb,Ntb)
   complex(kind=kind(0.D0))                 :: tmp(Ntb,Ntb)

   gfull=cmplx(0d0, 0d0, kind=kind(0.D0))
!!$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(SHARED)&
!!$OMP PRIVATE(ik,tmp,iat,ib,im,ineq)&
!!$OMP REDUCTION(+:gfull)
   do ik=1,Nk
     tmp(:,:)=-hk(ik,:,:)
     do iat=1,Nat
        ib=(iat-1)*Ndp
        do im=1,Ndp
           tmp(ib+im,ib+im)=tmp(ib+im,ib+im)+&
              cmplx(-dc(at2neq(iat)+1,im)+mu,w, kind=kind(0.D0))-sfw(at2neq(iat)+1,im)
        enddo
     enddo
     call linrg_cmplx(Ntb,tmp)
     gfull(:,:)=gfull(:,:)+tmp(:,:)
   enddo
!!$OMP END PARALLEL DO 
!fill in Gfw from Gfull (again assuming no band off-diagonal Sfw)
   do ineq=1,Nneq
     do im=1,Ndp
       glattiw(ineq,im)=gfull(im+Ndp*(neq2at(ineq)),im+Ndp*(neq2at(ineq)))/dble(Nk)
     enddo
   enddo
end subroutine get_gsumhk_para


!===============================================================================
subroutine get_Gk(Gk,Sfw,w,mu,Hk,at2neq,DC,NBands,NNEq,Nk,NAt,NTB)
!===============================================================================
!input
  integer,intent(in)         :: NBands,NNEq,Nk,NAt,NTB
  real(kind=kind(0.D0)),intent(in)          :: w,mu
  complex(kind=kind(0.D0)),intent(in)      :: Sfw(NNEq,NBands,2)
  complex(kind=kind(0.D0)),intent(in)      :: Hk(Nk,Nk,Nk,NTB,NTB)
  integer,intent(in)         :: at2neq(NAt)
  real(kind=kind(0.D0)),intent(in)          :: DC(NNeq,NBands,2)
!output
  complex(kind=kind(0.D0)),intent(out)     :: Gk(Nk,Nk,Nk,NTB,NTB)
!local
  integer                    :: iB,iM,ikx,iky,ikz,iat

  do ikx=1,Nk
  do iky=1,Nk
  do ikz=1,Nk
    Gk(ikx,iky,ikz,:,:)=-Hk(ikx,iky,ikz,:,:)
    do iat=1,NAt
      iB=(iat-1)*NBands
      do iM=1,NBands
        Gk(ikx,iky,ikz,iB+iM,iB+iM)=Gk(ikx,iky,ikz,iB+iM,iB+iM)+&
          cmplx(-DC(at2neq(iat)+1,iM,1)+mu,w, kind=kind(0.D0))-Sfw(at2neq(iat)+1,iM,1)
      enddo
    enddo
    call Linrg_cmplx(NTB,Gk(ikx,iky,ikz,:,:))
  enddo
  enddo 
  enddo 
end subroutine get_Gk



!===============================================================================
subroutine get_Polarization(Chi_q,Sfw,beta,w,mu,Hk,at2neq,DC,NBands,NNEq,Nk,NAt,NTB,Niw)
!===============================================================================
!input
  integer,intent(in)                       :: NBands,Niw,NNEq,Nk,NAt,NTB
  real(kind=kind(0.D0)),intent(in)         :: w(-Niw:Niw-1),mu,beta
  complex(kind=kind(0.D0)),intent(in)      :: Sfw(NNEq,NBands,2,-Niw:Niw-1)
  complex(kind=kind(0.D0)),intent(in)      :: Hk(Nk,Nk,Nk,NTB,NTB)
  integer,intent(in)                       :: at2neq(NAt)
  real(kind=kind(0.D0)),intent(in)         :: DC(NNeq,NBands,2)
!output
  real(kind=kind(0.D0)),intent(out)        :: Chi_q(NTB,NTB,0:Nk-1,0:Nk-1,0:Nk-1)
!local
  integer                        :: iw,ikx,iky,ikz,ikpqx,ikpqy,ikpqz,iqx,iqy,iqz,iM1,iM2
  complex(kind=kind(0.D0))                 :: Gk(0:Nk-1,0:Nk-1,0:Nk-1,NTB,NTB)

  Chi_q=0d0 
!$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(SHARED)&
!$OMP PRIVATE(iw,Gk,ikx,iky,ikz,iqx,iqy,iqz,ikpqx,ikpqy,ikpqz)&
!$OMP REDUCTION(+:Chi_q)
    do iw=0,Niw-1
      call get_Gk(Gk,Sfw(:,:,:,iw),w(iw),mu,Hk,at2neq,DC,NBands,NNEq,Nk,NAt,NTB)
      do iqx=0,Nk-1
      do iqy=0,Nk-1
      do iqz=0,Nk-1
        do ikx=0,Nk-1
        do iky=0,Nk-1
        do ikz=0,Nk-1
          ikpqx=ikx+iqx-sign(1,ikx+iqx-Nk/2)*Nk/2-Nk*((ikx+iqx-Nk/2)/Nk)
          ikpqy=iky+iqy-sign(1,iky+iqy-Nk/2)*Nk/2-Nk*((iky+iqy-Nk/2)/Nk)
          ikpqz=ikz+iqz-sign(1,ikz+iqz-Nk/2)*Nk/2-Nk*((ikz+iqz-Nk/2)/Nk)
          do iM1=1,NTB
          do iM2=1,NTB
            Chi_q(iM1,iM2,iqx,iqy,iqz)=Chi_q(iM1,iM2,iqx,iqy,iqz)+&
              (dble(Gk(ikx,iky,ikz,iM1,iM2))*dble(Gk(ikpqx,ikpqy,ikpqz,iM2,iM1))-&
              aimag(Gk(ikx,iky,ikz,iM1,iM2))*aimag(Gk(ikpqx,ikpqy,ikpqz,iM2,iM1)))
          enddo
          enddo
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
    enddo 
!$OMP END PARALLEL DO 
   Chi_q(:,:,:,:,:)=2d0/(beta*dble(Nk**3))*Chi_q(:,:,:,:,:)
   do iqx=0,Nk-1
   do iqy=0,Nk-1
   do iqz=0,Nk-1
      write(100,*)dble(iqx-Nk/2)/Nk,dble(iqy-Nk/2)/Nk,dble(iqz-Nk/2)/Nk
      do iM1=1,NTB
         write(100,'(3f17.10)')(Chi_q(iM1,iM2,iqx,iqy,iqz), iM2=1,NTB)
      enddo
   enddo
   enddo
   enddo

end subroutine get_Polarization

!===============================================================================
subroutine get_DOS(Glattiw,mu,wreal,Hk,neq2at,NBands,Niw,NNEq,Nk,NAt,NTB)
!===============================================================================
!input
   integer,intent(in)                        :: NBands,Niw,NNEq,Nk,NAt,NTB
   real(kind=kind(0.D0)),intent(in)          :: wreal(-Niw:Niw-1),mu
   complex(kind=kind(0.D0)),intent(in)       :: Hk(Nk,NTB,NTB)
   integer,intent(in)                        :: neq2at(NNeq)
!output
   real(kind=kind(0.D0)),intent(out)         :: Glattiw(NNEq,NBands,-Niw:Niw-1)
!local
   integer                                   :: iB,iw,ik,iM,ineq,iat
   complex(kind=kind(0.D0))                  :: GFull(NTB,NTB)
   complex(kind=kind(0.D0))                  :: tmp(NTB,NTB)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(iw,ik,tmp,GFull,iat,iB,iM,ineq)
     do iw=-Niw,Niw-1
       GFull=cmplx(0d0, 0d0, kind=kind(0.D0))
       do ik=1,Nk
         tmp(:,:)=-Hk(ik,:,:)
         do iat=1,NAt
            iB=(iat-1)*NBands
            do iM=1,NBands
               tmp(iB+iM,iB+iM)=tmp(iB+iM,iB+iM)+mu+wreal(iw)+cmplx(0.D0, 5D-2, kind(0.D0))!deltinho
            enddo
         enddo
         call Linrg_cmplx(NTB,tmp)
         GFull(:,:)=GFull(:,:)+tmp(:,:)
       enddo
!!fill in Gfw from Gfull (again assuming no band off-diagonal Sfw)
       do ineq=1,NNEq
         do iM=1,NBands
           Glattiw(ineq,iM,iw)=-1d0/dacos(-1d0)*aimag(Gfull(iM+NBands*(neq2at(ineq)),iM+NBands*(neq2at(ineq)))/dble(Nk))
         enddo
       enddo
     enddo 
!$OMP END PARALLEL DO
end subroutine get_DOS

!===============================================================================
subroutine get_fiwflat(fiw,w,mu,d,Neps,Niw)
!===============================================================================
!input
  integer,intent(in)       :: Neps,Niw
  real(kind=kind(0.D0)),intent(in)        :: w(-Niw:Niw-1),mu,d
!output  
  complex(kind=kind(0.D0)),intent(out)   :: fiw(-Niw:Niw-1)
!local
  integer                  :: iw,ieps
  real(kind=kind(0.D0))                   :: eps,deps

  fiw=cmplx(0d0, 0d0, kind=kind(0.D0))
  deps=(2d0*d/dble(Neps-1))
  do iw=-Niw,Niw-1
    do ieps=0,Neps-1
      eps=-dble(d)+deps*dble(ieps)
      fiw(iw)=fiw(iw)+deps/(2d0*d)/cmplx(mu-eps, w(iw), kind=kind(0.D0))
    enddo
  enddo
end subroutine get_fiwflat

!===============================================================================
subroutine get_fiwsemi(fiw,w,mu,d,Neps,Niw)
!===============================================================================
!input
  implicit none
  integer,intent(in)                     :: Neps, Niw
  real(kind=kind(0.D0)),intent(in)       :: w(-Niw:Niw-1)
  real(kind=kind(0.D0)), intent(in)      :: mu, d
!output
  complex(kind=kind(0.D0)),intent(out)   :: fiw(-Niw:Niw-1)
!local
  integer                                :: iw, ieps
  real(kind=kind(0.D0))                  :: pi=dacos(-1d0),eps,deps

  fiw=cmplx(0d0, 0d0, kind=kind(0.D0))
  deps=(2d0*d/dble(Neps-1))
  do iw=-Niw,Niw-1
    do ieps=0,Neps-1
      eps=-dble(d)+deps*dble(ieps)
      fiw(iw)=fiw(iw)+deps*2d0/(pi*d)*sqrt(1d0-(eps/d)**2)/cmplx(mu-eps, w(iw), kind=kind(0.D0))
    enddo 
  enddo 
end subroutine get_fiwsemi


!!===============================================================================
!subroutine get_gnano(gnano,mu,sfw,fiw,w,t,V,neq2at,at2neq,DC,Nneq,Nbands,Nat,Nlead)
!!===============================================================================
!!input
!  integer,intent(in)         :: Nneq,Nbands,Nat,Nlead
!  real*8,intent(in)          :: mu
!  complex*16,intent(in)      :: sfw(Nneq,Nbands,2)
!  complex*16,intent(in)      :: fiw(Nlead,2)
!  real*8,intent(in)          :: w
!  complex*16,intent(in)      :: t(Nat,Nbands,Nat,Nbands,2)
!  complex*16,intent(in)      :: V(Nat,Nbands,Nlead,2)
!  integer,intent(in)         :: neq2at(Nneq),at2neq(Nat)
!  real*8,intent(in)          :: DC(Nneq,Nbands,2)
!  !output
!  complex*16,intent(out)     :: gnano(Nneq,Nbands,2)
!  !local
!  integer                    :: Ntb
!  integer                    :: ib1,ib2,is1,is2,im1,im2,ineq,iat1,iat2,ilead
!  complex*16                 :: gfull(Nat*Nbands,Nat*Nbands)
!
!  Ntb=Nat*Nbands
!
!  gfull=cmplx(0d0,0d0, KIND=KIND(0.D0))
!  do iat1=1,NAt
!     ib1=(iat1-1)*NBands
!     do im1=1,NBands
!       do is1=1,2
!         do iat2=1,Nat
!           ib2=(iat2-1)*Nbands
!           do im2=1,Nbands
!             gfull(ib1+im1,ib2+im2)=-t(iat1,im1,iat2,im2,is1)
!             do ilead=1,Nlead
!               gfull(ib1+im1,ib2+im2)=gfull(ib1+im1,ib2+im2)-V(iat1,im1,ilead,is1)*V(iat2,im2,ilead,is1)*fiw(ilead,is1)
!             enddo
!           enddo
!         enddo
!         gfull(ib1+im1,ib1+im1)=gfull(ib1+im1,ib1+im1)+&
!           cmplx(mu-DC(at2neq(iat1)+1,im1,is1), w, KIND=KIND(0.D0))-sfw(at2neq(iat1)+1,im1,is1)
!       enddo
!     enddo
!   enddo
!   call Linrg_cmplx(Ntb,gfull)
!   do ineq=1,NNEq
!     do im1=1,NBands
!       do is1=1,2
!         gnano(ineq,im1,is1)=gfull(im1+NBands*(neq2at(ineq)),im1+NBands*(neq2at(ineq)))
!       enddo
!     enddo
!   enddo
!end subroutine get_gnano

!!===============================================================================
!subroutine get_w(w,beta,Niw)
!!===============================================================================
!!input
!   integer,intent(in)      :: Niw
!   real*8,intent(in)       :: beta
!!output
!   real*8,intent(out)      :: w(-Niw:Niw-1)
!!local
!   integer                 :: iw
!
!   do iw=-Niw,Niw-1
!      w(iw)=dacos(-1d0)*dble(2*iw+1)/beta
!   enddo
!end subroutine get_w

!!===============================================================================
!subroutine get_tau(tau,beta,NTau)
!!===============================================================================
!!input
!   integer,intent(in)      :: NTau
!   real*8,intent(in)       :: beta
!!output
!   real*8,intent(out)      :: tau(0:NTau-1)
!!local
!   integer                 :: itau
!
!   do itau=0,NTau-1
!      tau(itau)=beta/dble(NTau-1)*dble(itau)
!   enddo
!end subroutine get_tau

!!===============================================================================
!subroutine get_cttau(tau,beta,NTau)
!!===============================================================================
!! Here we write the tau-grid we use in the CTQMC simulation.
!! The first tau point is 0+dtau/2.
!!input
!   integer,intent(in)      :: NTau
!   real*8,intent(in)       :: beta
!!output
!   real*8,intent(out)      :: tau(1:NTau)
!!local
!   integer                 :: itau
!
!   do itau=1,NTau
!      tau(itau)=beta/dble(Ntau)*(dble(itau)-0.5)
!   enddo
!end subroutine get_cttau



!!===============================================================================
!subroutine get_Sfw(Sfw,Siw,Sfwhf,NNeq,Ndp,Nd,Niw)
!!===============================================================================
!!input
!   integer,intent(in)      :: Nd,Niw,NNeq,Ndp
!   complex*16,intent(in)   :: Siw(NNeq,Nd,2,-Niw:Niw-1)
!   complex*16,intent(in)   :: Sfwhf(NNeq,Ndp,Ndp,2,-Niw:Niw-1)
!!output
!   complex*16,intent(out)   :: Sfw(NNeq,Ndp,Ndp,2,-Niw:Niw-1)
!!local
!   integer                 :: iB
!
!   Sfw=Sfwhf
!   do iB=1,Nd
!      Sfw(:,iB,iB,:,:)=Sfw(:,iB,iB,:,:)+Siw(:,iB,:,:)
!   enddo
!end subroutine get_Sfw

!!===============================================================================
!subroutine get_Siw(Siw,Fiw,Giw,w,mu,Hmean,DC,NBands,Niw,NNEq,Ndp)
!!===============================================================================
!!input
!   integer,intent(in)      :: NBands,Niw,NNEq,Ndp
!   real*8,intent(in)       :: mu,w(-Niw:Niw-1)
!   real*8,intent(in)       :: Hmean(NNEq,Ndp,2)
!   real*8,intent(in)       :: DC(NNEq,Ndp,2)
!   complex*16,intent(in)   :: Fiw(NNEq,NBands,2,-Niw:Niw-1)
!   complex*16,intent(in)   :: Giw(NNEq,NBands,2,-Niw:Niw-1)
!!output
!   complex*16,intent(out)  :: Siw(NNEq,NBands,2,-Niw:Niw-1)
!!local
!   integer                 :: iw
!
!   do iw=-Niw,Niw-1
!      Siw(:,:,:,iw)=-Fiw(:,:,:,-iw-1)+cmplx(mu, w(iw), KIND=KIND(0.D0))-Hmean(:,:NBands,:)-DC(:,:NBands,:)-1d0/Giw(:,:NBands,:,iw)
!   enddo
!end subroutine get_Siw

!!===============================================================================
!subroutine get_fiw(fiw,sfw,glattiw,w,mu,hmean,dc,Nd,Niw,Nneq,Ndp)
!!===============================================================================
!!input
!   integer,intent(in)      :: Nd,Niw,Nneq,Ndp
!   real*8,intent(in)       :: mu,w(-Niw:Niw-1)
!   real*8,intent(in)       :: hmean(Nneq,Ndp,2)
!   real*8,intent(in)       :: dc(Nneq,Ndp,2)
!   complex*16,intent(in)   :: sfw(Nneq,Ndp,Ndp,2,-Niw:Niw-1)
!   complex*16,intent(in)   :: glattiw(Nneq,Ndp,2,-Niw:Niw-1)
!!output
!   complex*16,intent(out)  :: fiw(Nneq,Nd,2,-Niw:Niw-1)
!!local
!   integer                 :: iw,ib
!
!   do iw=-Niw,Niw-1
!      do ib=1,Nd
!        fiw(:,ib,:,-iw-1)=cmplx(mu, w(iw), KIND=KIND(0.D0))-hmean(:,ib,:)-dc(:,ib,:)-sfw(:,ib,ib,:,iw)-1d0/glattiw(:,ib,:,iw)
!      enddo
!   enddo
!end subroutine get_fiw

!===============================================================================
subroutine Linrg_cmplx(L,A)
!===============================================================================
!     MatrixInversion A = A ^ -1. A: LxL Matrix 
!     In/Output Parameters
      integer,intent(in)          :: L
      complex(KIND=KIND(0.D0))    :: A(L,L)
!     Local Parameters
      integer                     :: ipiv(L)
      complex(KIND=KIND(0.D0))    :: work(L*L*L*L)
      integer                     :: INFO
#ifdef LAPACK77_Interface      
      call ZGETRF(L,L,A,L,ipiv,INFO)
      call ZGETRI(L,A,L,ipiv,work,L*L*L*L,INFO)
#endif
#ifdef LAPACK95_Interface
      call GETRF(N,ipiv,info)
      call GETRI(N,ipiv,info)
#endif
      if(INFO.ne.0) then
         write(*,'("Error in Linrg_cmplx probably a singular matrix: ", I5)') INFO
         STOP
      end if
end subroutine Linrg_cmplx

!===============================================================================
subroutine inverse(B,A,L)
!===============================================================================
!     MatrixInversion A = A ^ -1. A: LxL Matrix 
!     In/Output Parameters
      integer,intent(in):: L
      complex(KIND=KIND(0.D0)), intent(in)  :: A(L,L)
      complex(KIND=KIND(0.D0)), intent(out) :: B(L,L)
!     Local Parameters
      integer                               :: ipiv(L)
      complex(KIND=KIND(0.D0))              :: work(2*L*L)
      integer                               :: INFO

#ifdef LAPACK77_Interface
      call ZLACPY('A', L, L, A, L, B, L)
      call ZGETRF(L,L,B,L,ipiv,INFO)
      call ZGETRI(L,B,L,ipiv,work,2*L*L,INFO)
#endif
#ifdef LAPACK95_Interface
      B = A
      call GETRF(N,ipiv,info)
      call GETRI(N,ipiv,info)
#endif
      if(INFO.ne.0) then
         write(*,'("Error in Linrg_cmplx probably a singular matrix: ", I5)') INFO
         STOP
      end if
end subroutine inverse


!===============================================================================
subroutine eigvals_cmplx(x,A,L)
!===============================================================================
!     Eigenvals for A: LxL complex matrix
!     MatrixInversion A = A ^ -1. A: LxL Matrix 
!     In/Output Parameters
      integer,intent(in)                    :: L
      complex(KIND=KIND(0.D0)), intent(in)  :: A(L,L)
      complex(KIND=KIND(0.D0)), intent(out) :: x(L)
!     Local Parameters
      complex(KIND=KIND(0.D0))              :: tmp(L,L)
      complex(KIND=KIND(0.D0))              :: vl(L,L), vr(L,L)
      complex(KIND=KIND(0.D0))              :: work(2*L*L)
      complex(KIND=KIND(0.D0))              :: rwork(2*L)
      integer                               :: INFO

      if(L == 1) then
         x = A(1,1)
         return
      endif

#ifdef LAPACK77_Interface
      call ZLACPY('A', L, L, A, L, tmp, L)
      call ZGEEV('N','N',L,tmp,L,x,vl,L,vr,L,work,2*L*L,rwork,INFO)
#endif
#ifdef LAPACK95_Interface
      tmp = A
      call LA_GEEV(tmp,x,INFO)
#endif
      if(INFO.ne.0) then
         write(*,'("Error in eigvals_cmplx: ", I5)') INFO
         STOP
      end if
end subroutine eigvals_cmplx

