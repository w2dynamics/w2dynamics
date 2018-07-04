!*****************************************************************************************
!*****************************************************************************************
!***
!***    The [basic idea] behind this code is, to find a spectral function (specFunc=A)
!***    which is obtained by Bayes' theorem maximizing the posteriori probability distribution
!***    for the spectral function, A(w), given a spectral model, m(w): 
!***    Pr(A|G)=Pr(G|A)*Pr(A), where Pr(G|A) is the likelihood function of the data 
!***    (greenFunc -> Pr(G|A)=DEXP(-chi²/2)) given the spectral function.
!***    Pr(A) is the prior knowledge about the spectral function concerning the model m
!***    in an entropic like fashion (DEXP(alpha*S)), where S=-integral(A-m-A*DLOG(A/m))dw 
!***    is the entropy (Ref[1,2,3]).
!***
!***    The [evaluation of specFunc, A], based on Bayesian statistics is done in an iterative
!***    fashion in subroutine <w2maxent> where first a guess for alpha and A is made and then Pr(A|G) 
!***    is maximized (here by an annealing Monte Carlo method) according to Bayes' theorem.
!***    The [parameter alpha] is determined via Bayes-statistics according to Ref[1,2] and
!***    this also is the key of the algorithm: alpha and A have to satisfy the equation
!***    { -2 alpha*S = Trace ( Lambda(A)_ij*(alpha*1_ij+Lambda(A)_ij)^-1 ) > }, 
!***    see Eq. 4.28 in Ref[1] and Eq. 22 in Ref[2], if they don't we alter alpha and 
!***    start again with maximizing Pr(A|G) leading (perhaps) to a new alpha and we repeat this
!***    procedure until convergence is reached.
!***
!***    MODULE MMaximumEntropy contains parameter and subroutines necessary to calculate the most probable spectral function
!***
!*** Reference [1] ... Physics Reports 269 (1996) 133-195
!***                   BAYESIAN INFERENCE AND THE ANALYTIC CONTINUATION OF IMAGINARY-TIME QUANTUM MONTE CARLO DATA
!***                   Mark JARRELL, J.E. GUBERNATIS
!*** Reference [2] ... Physical Review B, Volume 44, Number 12, 1991
!***                   Quantum MonteCarlo simulations and masimum entropy: Dynamics from imaginary-time data
!***                   Mark JARRELL, J.E. GUBERNATIS
!*** Reference [3] ... Physical Review B, Volume 57, Number 17, 1998
!***                   Stochastic method for analytic continuation of quantum Monte Carlo data
!***                   A.W. Sandvik
!*** Reference [4] ... Project thesis Feb-Aril 2015
!***                   Benedikt Hartl
!***
!*** The code is based on the former MaxEnt code written by Anders W. Sandvik but is
!*** reorderend and partly rewritten. Passages which where copy-pasted are explicitly 
!*** referenced to as such (method by Sandvik). 
!*** Two small bugs where fixed ( the smoothing process caused an array-reference out of
!*** it's definition area i\in[1,Nspec], the variation procedure in the MonteCarlo part
!*** did not consider the last index i=Nspec ).
!*** The parameter setting concerning annealing parameters in the Monte Carlo section,
!*** variation- and convergence criteria of the Bayesian inference method are now accessible
!*** by a list of parameter in the module.
!*** 
!*** The usage of imaginary-time Green's function (kernelMode=0)
!*** was extended to the usage of Matrsubara-frequency Green's function (kernelMode=1)
!*** and two-particle Green's function for susceptibility-calculations 
!*** (kernelMode=10 ... Bosonic kernel, kernelMode=11 ... Bosonic spectral function)
!***
!*****************************************************************************************
MODULE MMaximumEntropy
  IMPLICIT NONE

  INTEGER, PARAMETER         :: KINDR=KIND(0.D0)
  INTEGER, PARAMETER         :: KINDC=KIND(0.D0)

  !*** PARAMETER FOR BAYESIAN INVERENCE METHODE IN SUBROUTINE <bayesianInverence>
  REAL(KINDR)      :: bayesianConvergence         = 0.05  ! numerical accuracy for the Bayesian convergence criterion [1] Eq. (4.28): {-2alpha*S / Tr[...] - 1 < bayesianConvergence } ?
  REAL(KINDR)      :: bayesianUpdateAlpha         = 0.05  ! update factor for  Bayesian parameter alpha, see subroutine bayesianInverence
  !*** PARAMETER FOR BAYESIAN INVERENCE METHODE IN SUBROUTINE <bayesianInverence>

  !*** Metropolis and annealing parameters for determination of spectral function via Monte Carlo
  REAL(KINDR)      :: noiseFilter                 = 1D-1  ! noiseFilter, only variations above the noise {noiseFilter*max(A(omega))} are considered as important for try and acceptance rate in MonteCarlo run

  REAL(KINDR)      :: variationUpdate             = 1.5   ! magnitude for important sampling: varSpecFunc(omega) = variationStrength * noiseFilter *(spectralFunction(omega) + noiseFilter * max(spectralFunction(omega))):
                                                     ! variationStrength' = variationStrength * variationUpdate, methode by Sandivk
  REAL(KINDR)      :: criticaleAcceptRate         = 1D-1  ! critical parameter for performing update of variationStrength with respect to try/acceptance-rate in MonteCarlo Algorithm
  REAL(KINDR)      :: criticaleVariationInc       = 1D-2  ! critical parameter for how to update variational-strength-parameter
  REAL(KINDR)      :: criticaleVariationDec       = 1D-3  ! critical parameter for how to update variational-strength-parameter 

  REAL(KINDR)      :: variationStrengthWarmpup    = 1.    ! initial variational magnitude for MonteCarlo search of most probable spectrum in the warm-up-run (starting with initial randomly choosen spectral function)
  REAL(KINDR)      :: variationStrength           = 0.05  ! initial variational magnitude for MonteCarlo search of most probable spectrum in the main loop starting from spectral function of previous MonteCarlo-run
  REAL(KINDR)      :: variationStrengthSmooth     = 0.005 ! initial variational magnitude for MonteCarlo search of most probable spectrum in the smoothing procedure after convergence of spectral function is reached

  REAL(KINDR)      :: annealingUpdate             = 1.5   ! update factor for annealing parameter exp((-L+alpha*S)/Theta): Theta' = Theta * annealingUpdate, [3], methode by Sandivk
  REAL(KINDR)      :: annealingParamWarmup        = 10D0  ! initial annealing temperature for warm-up of a randomly choosen spectral function into a more probable one w.r.t. the data
  REAL(KINDR)      :: annealingParam              = 0.001 ! initial annealing temperature for main-loop of finding the most probable spectral function
  REAL(KINDR)      :: annealingParamSmooth        = 0.005 ! initial annealing temperature for smoothing procedure after convergence of the spectral function is reached
  !*** Metropolis and annealing parameters for determination of spectral function via Montecarlo


  REAL(KINDR)      :: kernelExponentCutoff        = 40D0  ! numerical cut-off when exponent of exponentially damped kernels (exp(-40)<<1) is considered to by asymptotic,
                                                     ! numerically limit to stabalize calculations exp(+40)>>1 (which can happen for negative frequencies) 
  REAL(KINDR)      :: accuracy                    = 1D-12 ! numerical accuracy used e.g. in the MonteCarlo-Subroutine montecarloSpectralFunction for numerical calculations of logarithm 
  
  INTEGER     :: integrationMethod           = 1 ! 1: trapezoid, else: assume equidistant (default)

  INTEGER      :: outputUnit                 = 6                     ! 6: terminal 
  CHARACTER(len=80) :: logMaxent                  = '_maxent.log'         ! postfix for maxent-log file 
  CHARACTER(len=80) :: logBayes                   = '_bayesInverence.log' ! postfix for bayes-log file 
  CHARACTER(len=80) :: logMontecarlo              = '_Montecarlo.log'     ! postfix for monte-carlo-log file

  LOGICAL     :: debug_maxent                = .TRUE.  ! shows maxent output on terminal if true
  LOGICAL     :: debug_mc                    = .FALSE. ! show monte-carlo in terminal
  LOGICAL     :: debug_bayes                 = .FALSE. ! show bayes-results in terminal     
  
  !*** in Benedikts original this is seperately allocatable stuff ***
  REAL(KINDR),DIMENSION(:),ALLOCATABLE           :: maxentCorrFunc
  REAL(KINDR),DIMENSION(:),ALLOCATABLE           :: maxentInversCov ! is the inverse covariance matrix of the measured Green's function.
  REAL(KINDR),DIMENSION(:,:),ALLOCATABLE         :: kernel          ! kernel, satisfying G(j)=Integral(kernel_ji(w,tau,beta,...)*SpecFunc_i(w) dw), Eq.(3.1)Ref[1]
  
  CONTAINS

!*****************************************************************************************
!*** w2maxent ****************************************************************************
!*****************************************************************************************
!***
!*** Calculates most probable spectral Function (specFunc) defined on a grid (specGrid) for a
!*** given Green's function (greenFunc) defined on an imaginary time or Matsubara-frequency grid 
!*** (greenGrid) for a given model (specMod) of the spectral function by means of Bayesian-statistics:
!***
!***
!*** PARAMETERS:
!***  specFunc is the spectral function A(w) (object of interest) which is in an integral relation to the data Eq.(3.1) in Ref[1].
!***  specModel supplies prior knowledge about the spectral Function. 
!***    A default choice for specModel would be a flat distribution (for kernel 0 and 1) but one can also choose 
!***    different models (e.g. results according to specific problems from perturbation theory). 
!***    This may be important for high frequency behavior of the spectral function, because 
!***    the likelihood function (based on Integral(Kernel*SpecFunc)) is insensetive in this region.
!***    Kernel 10 and 11 do not have a good, predefined default model, consider Ref[1,4]
!***  specGrid is the spectral grid where specFunc is defined (note w>=0 for kernel 10,11)
!***  Nspec defines the number of spectral grid points.
!***  greenFunc is the data-vector, measured values of the Green's function (e.q. by CT-QMC) and is a complex
!***            valued quantity.
!***            For kernel 0,10,11 only the real part of greenFunc is considered in calculation
!***            For kernel 1 (Matrsubara-kernel) the real- and imaginary part of the data are separately taken care of
!***  greenInverseCov: greenFunc exhibits errors due to the numerical evaluation (CT-QMC) which 
!***                   enter calculation in the inverse covariance matrix here assumed to be diagonal 
!***                   (greenInverseCov(:)=1/sigma(i)^2, no correlations between different grid-points).
!***                   Beware of very small errors, singular values (1/0), and very big errors in the preparation of greenInverseCov(:)
!***  greenGrid is the grid, where greenFunc is defined (either imaginary-time for kernels 0,10 and 11 or Matsubara-Frequencies for kernel 1
!***  NGreen defines the number of Green's function grid points.
!***  alphaPrior defines the initial guess for the Bayesian parameter alpha (should be rather big ~ 1000)
!***  kernelMode defines, which kind of Green's function we are dealing with (numbers correspond to Ref[4]):
!***             kernelMode=0       imaginary  Time grid
!***             kernelMode=1       Matsubara-frequency grid (Fermionic case implemented, pole in kernel at omega=0 not checked for Bosons)
!***             kernelMode=10      imaginary Time grid for two particle Green's function, Bosonic kernel
!***             kernelMode=11      imaginary Time grid for two particle Green's function, Bosonic spectral function
!***  kernelBeta defines inverse temperature from experiment (1/[kb*T]) 
!***  NMonteCarlo defines the number of Monte Carlo steps for numerical search of specFunc.
!***  NAnnealing is the number of Monte Carlo steps until the annealing parameter is updated.
!***  seed initialization point for random-number generator
!***  NSmooth defines the number of smoothing steps which are performed after Bayesian convergence is reached 
!***          different smoothing procedures are available through the parameter smoothingMethod=0,1,2
!***  smoothCnt defines number of spectral grid points, where the average is taken over
!***    smooth(i)=SUM(specFunc(i-smoothCnt : i+smoothCnt))/(i+smoothCnt-(i-smoothCnt) + 1) for local average method
!***  smoothingMethod defines the method, how the local average over the (noisy) spectral function is performed
!***                   0: locally average spectral function: smoothSpectrum(omega_i)=sum(spectralFunction(k:l))/(l-k+1), k=omega_i-smoothCnt, l=omega_i+smoothCnt, 
!***                   1: smoothing by locally integrate spectral function: smoothSpectrum(omega_i)=SUM(domega(k:l)*spectralFunction(k:l))/SUM(domgea(k:l), important for non-equidistant grid-points 
!***                   2: smoothing by locally integrate chi''(omega) propTo spectral function, important for non-equidistant grid points for kernel 10 and 11
!***  logPrefix defines the prefix of the DLOG-files for the Monte Carlo runs and the Bayesian convergence procedure.
!***
!*****************************************************************************************
!*****************************************************************************************
SUBROUTINE w2maxent ( specFunc, specModel, specGrid, NSpec, &
                      greenFunc, greenInverseCov, greenGrid, NGreen, &
                      alphaPrior, kernelMode, kernelBeta, &
                      NMonteCarlo, NAnnealing,seed, NSmooth, smoothCnt, &
                      smoothingMethod, logPrefix)
  USE MMersenneTwister                                      !random number generator
  IMPLICIT NONE

  !***************************************************************************************
  !*** Spectral Function Variables ***
  INTEGER,                      INTENT(IN)  :: NSpec           ! number of grid points for spectral function
  REAL(KINDR),     DIMENSION(Nspec), INTENT(INOUT) :: specFunc        ! spectral function, object of interest
  REAL(KINDR),     DIMENSION(Nspec)              :: specModel       ! (spectral) model for spectral function (enters in entropy)
  REAL(KINDR),     DIMENSION(Nspec)              :: specGrid        ! grid for spectral function
  REAL(KINDR),     DIMENSION(Nspec)              :: specUnit        ! grid weight for numerical integration
  !*** Spectral Function Variables ***

  !*** Green's Function Variables ***
  INTEGER,                       INTENT(IN) :: NGreen          ! number of grid points for Green's function
  COMPLEX(KINDC),  DIMENSION(Ngreen), INTENT(IN) :: greenFunc       ! Green's function from DEXPeriment / CT-QMC
  REAL(KINDR),     DIMENSION(Ngreen), INTENT(IN) :: greenInverseCov ! inverse covariance matrix (1/sigma(i)²), assumed to be diagonal
  REAL(KINDR),     DIMENSION(Ngreen), INTENT(IN) :: greenGrid       ! grid for Green's function

  INTEGER                                   :: maxentNCorr
  !*** Green's Function Variables ***
  
  !*** Kernel Variables ***
  INTEGER,                    INTENT(IN)    :: kernelMode      ! defines specific shape of kernel (0=imaginary time, )
  REAL(KINDR),                     INTENT(IN)    :: kernelBeta      ! beta from DEXPeriment (1/[kb*T]) which enters the kernel in imaginary time Green's function.
  !*** Kernel Variables ***

  !*** Maximum Entropy (Bayesian) Parameter ***
  REAL(KINDR),                       INTENT(IN)  :: alphaPrior          ! initial guess for alpha entering the priori pdf: DEXP(alpha*S), S being the entropy
  REAL(KINDR)                                    :: alpha               ! iteratively (hopefully convergent w.r.t. Eq. 4.28 in Ref[1]) alpha
  REAL(KINDR),     DIMENSION(Nspec,Nspec)        :: likelihoodCurvature ! curvature tensor of likelihood DEXPonent w.r.t. spectral Function, Eq. 4.9 in Ref[1]
  REAL(KINDR)                                    :: bayesPostExponent   ! = -2 alpha S, according to Eq. 4.28 in Ref[1]
  REAL(KINDR)                                    :: bayesTrace          ! = Trace ( Lambda(A)_ij*(alpha*1_ij+Lambda(A)_ij)^-1 ), Eq. 4.29, 4.30 in Ref[1]
                                                                !   equals number of good data ...o 
  REAL(KINDR)                                    :: chiSquared          ! least mean squares value obtained from Monte Carlo algorithm
  !*** Maximum Entropy (Bayesian) Parameter ***
  
  !*** Monte Carlo & Annelaling Parameter ***
  INTEGER,                      INTENT(IN)  :: NMonteCarlo     ! number of Monte Carlo steps for maximizing posteriori pdf for the spectral function Pr(A|G)
  INTEGER,                      INTENT(IN)  :: NAnnealing      ! number of Monte Carlo steps for one annealing (temperature) step
  INTEGER,                    INTENT(IN)  :: seed            ! initial seed for randomn number generator, 4 byte integer
  !*** Monte Carlo & Annelaling Parameter ***
  
  !*** Smoothening Parameter ***
  INTEGER,                      INTENT(IN)  :: NSmooth         ! number of smoothing steps after Bayesian convergence is reached (w.r.t. Eq. 4.28 in Ref[1])
  INTEGER,                      INTENT(IN)  :: smoothCnt       ! number of spectral grid points, where the average is taken 
  INTEGER,                      INTENT(IN)  :: smoothingMethod ! defines the method, how the local average over the (noisy) spectral function is performed
                                                               !   smooth(i)=SUM(specFunc(i-smoothCnt : i+smoothCnt))/(i+smoothCnt-(i-smoothCnt) + 1) for local average
  !*** Smoothening Parameter ***

  !*** Susceptibility values - only for Kernel 1 and 2 ***
  REAL(KINDR)                                    :: chiOmegaZero    ! real part of spectral susceptibility (via Kramers Kronig Relation)
  REAL(KINDR)                                    :: chiTauZero      ! susceptibility at origin of imaginary time (tau=0 in EQ.5.14 REf[1])
  !*** Susceptibility values - only for Kernel 1 and 2 ***

  !*** File IO Parameter ***
  CHARACTER(*),                 INTENT(IN)  :: logPrefix       ! prefix for DLOG-files for Monte Carlo & Bayesian Inverence
  CHARACTER(len=256)                             :: logMsg          ! DLOG message string
  !*** File IO Parameter ***
  
  !*** Local Helpers ***
  INTEGER                                   :: i_smooth,i_montecarlo,i_bayes ! running index for smoothing, Monte Carlo steps, Bayes Steps
  LOGICAL                                   :: converged       ! defines if Bayesian convergence is reached  (w.r.t. Eq. 4.28 in Ref[1])
  LOGICAL                                   :: working         ! defines if main loop is finished (convergence and smoothing)
  INTEGER                                   :: i               ! loop and index dummies
  INTEGER,DIMENSION(8)                      :: dateTime        ! timestamp-helper
  
  !*** Local Helpers ***
  !***************************************************************************************
    
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                                                  &
       TRIM('####################################################################################'),&
       'rewind',debug_maxent)
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                                                  &
       TRIM('### w2maxent #######################################################################'),&
       'rewind',debug_maxent)

  CALL writeLog (TRIM(logPrefix)//TRIM(logBayes),                                                   &
       TRIM('######################################################################'),&
       'rewind',debug_bayes)
  CALL writeLog (TRIM(logPrefix)//TRIM(logBayes),                                                   &
       TRIM('### w2maxent #########################################################'),&
       'rewind',debug_bayes)

  CALL writeLog (TRIM(logPrefix)//TRIM(logMontecarlo),                                              &
       TRIM('#################################################################################'),&
       'rewind',debug_mc)
  CALL writeLog (TRIM(logPrefix)//TRIM(logMontecarlo),                                              &
       TRIM('### w2maxent ####################################################################'),&
       'rewind',debug_mc)  
  
  !***************************************************************************************
  !*** print timestamp
  CALL DATE_AND_TIME(VALUES=dateTime)
  WRITE(logMsg,"(A16,I0.4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.3)")       &
               '### TIMESTAMP : ', dateTime(1),'.',dateTime(2),'.',dateTime(3),&
                              '-', dateTime(5),':',dateTime(6),':',dateTime(7),'.',dateTime(8)
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent)    ,TRIM(logMsg),'append',debug_maxent)
  CALL writeLog (TRIM(logPrefix)//TRIM(logBayes)     ,TRIM(logMsg),'append',debug_bayes)
  CALL writeLog (TRIM(logPrefix)//TRIM(logMontecarlo),TRIM(logMsg),'append',debug_mc)

  !***************************************************************************************
  !*** initializing specral function  
  CALL sgrnd(seed)
  
  WRITE(logMsg,*) '### Init spectral function : randomize (0,1]'
  DO i=1,NSpec
    specFunc(i)=grnd()                    !initially randomize spectral function
    IF(specFunc(i) .EQ. 0D0) CYCLE        !explicitely avoid 0D0 from MMersenneTwister
  END DO
  
  !*** define grid Units --> spectral weights for numerical integration
  CALL getGridUnit(specGrid,specUnit,Nspec)
  
  !*** normalizing spectral function 
  CALL normalizeToGrid(specFunc,specGrid,NSpec,kernelMode)
  IF(kernelMode .EQ. 0 .OR. kernelMode .EQ. 1) THEN
    WRITE(logMsg,"(A,F10.5)") "### Init spectral function : normalize to --> Sum[a(wi) dwi]=1:", &
                         SUM(specFunc(:)*specUnit(:))
    CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM(logMsg),'append',debug_maxent)
  END IF
  
  !*** normalizing spectral model, simplyfies calculation of entropy
  CALL normalizeToGrid(specModel,specGrid,NSpec,kernelMode)
  IF(kernelMode .EQ. 0 .OR. kernelMode .EQ. 1) THEN
    WRITE(logMsg,"(A,F10.5)") "### Init spectral model    : normalize to --> Sum[m(wi) dwi]=1:", &
                       SUM(specModel(:)*specUnit(:))
    CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM(logMsg),'append',debug_maxent)
  END IF
  
  !*** initialize kernel with mode kernelMode
  WRITE(logMsg,"(A17,I2)") '### kernelMode = ', kernelMode
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM(logMsg),'append',debug_maxent)
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM("### Init kernel k_ij=K_ij*dj"),'append',debug_maxent)
  CALL initKernel(specGrid,specUnit,NSpec,                         &
                  greenFunc,greenInverseCov,greenGrid,NGreen,      &
                  kernelMode,kernelBeta,maxentNCorr)   
  
  !*** initialize kernel with mode kernelMode
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                                                  &
                TRIM("### Init Likelihood curvature: (d^2 L / dA_i dA_j) = Sum_kl K_ki (C^-1)_kl K_lj"), &
                'append',debug_maxent)
  CALL initLikelihoodCurvature(likelihoodCurvature,specUnit,NSpec,maxentNCorr)
  !***************************************************************************************  
  
            
  !***************************************************************************************  
  !*** WARMUP spectral function
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                                                  &
                TRIM("### warmup spectral function ..."), &
                'append',debug_maxent)
  
  !*** logging - headliner - montecarlo.log
  WRITE(logMsg,"(A14,A15,A9,A8,A15,A20)") "#   mc(step)  ","chi2","acc","try","variatInit","annealInit"
  CALL writeLog (TRIM(logPrefix)//TRIM(logMontecarlo),TRIM(logMsg),'append',debug_mc)
  
  !*** warmup Monte Carlo run for random (but normalized) spectral function and alphaPrior
  alpha=alphaPrior  !initial guess for alpha
  i_montecarlo=1    ! index for logFile entry

  
  CALL montecarloSpectralFunction(specFunc, specModel, specUnit,                         &
                                  alpha, NSpec, maxentNCorr,                             & 
                                  NMonteCarlo, NAnnealing,                               &
                                  annealingParamWarmup,variationStrengthWarmpup,         & !initial annealing parameter
                                  chiSquared,logPrefix, i_montecarlo)
  !*** WARMUP spectral function
  !***************************************************************************************  
  
  
  !***************************************************************************************  
  !*** WARMUP bayesian parameter alpha
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                                                  &
                TRIM("### Start with search for most probable spectral function"), &
                'append',debug_maxent)

  !*** logging - headliner - bayes.log
  WRITE(logMsg,"(A5,A15,A15,A15,A20)") '#   i', 'alpha', 'Entropy=S', 'trace', '-2alpha*S/trace'
  CALL writeLog (TRIM(logPrefix)//TRIM(logBayes),TRIM(logMsg),'append',debug_bayes)

  WRITE(logMsg,"(A17,A12,A10,A16,A13,A10,A6)") '#      chiSquared',  &
                                        'Entropy=S', 'trace', '-2alpha*S/trace', 'alpha', 'action','Norm'
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM(logMsg),'append',debug_maxent)
  
  !*** make first check of Bayesian convergence with spectral function from warmup Monte Carlo run
  !*** and (most likely) update alpha
  i_bayes=1           !first DLOG-entry for warmup-run
  converged=.FALSE.   !init convergence - control logical
  CALL bayesianInverence (bayesTrace, bayesPostExponent, converged, alpha, &
                          likelihoodCurvature, specFunc, specModel, specUnit, &
                          NSpec,logPrefix,i_bayes,chiSquared,'warmup')
  ! converged is controled by bayesianInverence and itself controls the working branch  
  ! in the main loop (iteration, smoothing)
  !*** WARMUP bayesian parameter alpha
  !***************************************************************************************  


  !***************************************************************************************
  !*** MAIN LOOP * iteration until Bayesian convergence is reached and smoothining applied
  working   = .TRUE.  !controls if loop DEXPires or if there is still work to do (iteration, smoothing)
  converged = .FALSE. !make at least one more MonteCarlo run with lower annealing temperature
  i_smooth  = 1       !DLOG-index for Monte Carlo runs
  
  DO WHILE (working)  ! MAIN LOOP (endless loop, controlled by Bayesian convergence Eq. 4.28 in Ref[1] 
                      !            and by the convergence-level MODULE parameter bayesianConvergence)

    IF (converged) THEN  ! converged is controled by the SUBROUTINE bayesianInverence !
      !*** smoothing branch
      IF (i_smooth .LE. NSmooth) THEN                 ! converged, now proceed with smoothing

        !*** smoothing procedure        
        CALL smoothSpectrum(specFunc, specGrid, NSpec,smoothCnt,smoothingMethod,kernelMode,kernelBeta)        
        !*** smoothing procedure
        
        !*** apply Monte Carlo method to maximize posteriori probability with converged alpha
        !*** (smoothing can change a lot in the chi^2 value !)
        CALL montecarloSpectralFunction(specFunc, specModel, specUnit,                           & 
                                        alpha, NSpec, maxentNCorr, NMonteCarlo,                  &
                                        NAnnealing, annealingParamSmooth,variationStrengthSmooth,&
                                        chiSquared,logPrefix, i_montecarlo)
        
        !*** in principal one should also apply the Bayesian convergence again, but we run into 
        !*** the thread of an endless loop ...
        CALL bayesianInverence (bayesTrace, bayesPostExponent, converged, alpha, &    ! TODO, converging for smoothed spectrum: replace i_smooth .LE. NSmooth condition with logical,
                                likelihoodCurvature, specFunc, specModel, specUnit, & !       converged -> smoothConverged
                                NSpec,logPrefix,i_bayes,chiSquared,'smoothing')
        converged=.TRUE. ! therefore we artificially set converged = .TRUE. so alpha won't be updated by bayesianInverence
        
        i_smooth=i_smooth+1    ! go to next smoothing step
      ELSE                     ! smoothing finished, i_smooth > NSmooth

        CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                       &
                       TRIM("### Finished with search."),                      &
                       'append',debug_maxent)
        working=.FALSE.        ! main loop will exit

      END IF
      !*** smoothing branche      
    ELSE  
      !*** Bayesian convergence branch, converged .EQV. .FALSE.
      
      !*** apply Monte Carlo method to maximize posteriori probability with current alpha
      CALL montecarloSpectralFunction(specFunc, specModel, specUnit,                           &
                                      alpha, NSpec, maxentNCorr, NMonteCarlo,                  &
                                      NAnnealing, annealingParam,variationStrength,            &
                                      chiSquared,logPrefix, i_montecarlo)
                   
      !*** apply Bayesian convergence method to eighter update alpha or determine convergence
      !*** w.r.t. Eq. 4.28 in Ref[1] 
      CALL bayesianInverence (bayesTrace, bayesPostExponent, converged, alpha,    &
                              likelihoodCurvature, specFunc, specModel, specUnit, &
                              NSpec,logPrefix,i_bayes,chiSquared,'converge')
      !*** numerical convergence criterion defined in bayesianInverence(...)                             
      IF (converged) THEN ! fihished annealing
        i_smooth=1        ! start with smoothing
      END IF
      
      !*** Bayesian convergence branch
    END IF    
    
  END DO
  !*** MAIN LOOP * iteration until Bayesian convergence is reached and smoothining applied
  !***************************************************************************************
  

  !***************************************************************************************
  !*** logging user instructions 
  
  !*** convergence criterion
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM('###'),'append',debug_maxent)

  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent), &
                 TRIM('### 1.) Check convergence: -2alpha*S/trace --> 1'),'append',debug_maxent)
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent), &
                 TRIM('###     (also during smoothing steps)'),'append',debug_maxent)
  
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM('###'),'append',debug_maxent)
  !*** convergence criterion
  
  !*** chiSquared measure
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent), &
                 TRIM('### 2.) Check initial error estimate'),'append',debug_maxent)                                                               
                 
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM('###'),'append',debug_maxent)
  
  WRITE(logMsg,"(A,F10.3,A)") '###     normalized mean least square: chiSquared/N_data = ', &
        chiSquared/maxentNCorr, ', should be O(1)'
  
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM(logMsg),'append',debug_maxent)
  !*** chiSquared measure
  
  !*** error rescaling factor from Eq. 4.47 in Ref[1], Ng = bayesTrace -> Eq. 4.29, 4.30 in Ref[1]
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM('###'),'append',debug_maxent)
  
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                                                         &
                 TRIM('###     Error rescaling factor: quantifies, if error estimate of data is correct'), &
                 'append',debug_maxent)
                 
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent) ,                                           &
                 TRIM('###     gamma^2 = chiSquared/(N_data - N_gooddata), should be O(1)'), &
                 'append',debug_maxent)
                 
  WRITE(logMsg,"(A,F10.3)") '###     gamma= ', DSQRT(chiSquared/(maxentNCorr-bayesTrace))
  
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent) ,TRIM(logMsg),'append',debug_maxent)

  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM('###'),'append',debug_maxent)
  
  !*** error rescaling factor from Eq. 4.47 in Ref[1], Ng = bayesTrace -> Eq. 4.29, 4.30 in Ref[1]

  !*** calculate moments of susceptibility, Eq. 5.15 and 5.16 in Ref[1] for kernel 10 and 11
  IF (kernelMode .EQ. 0 .OR. kernelMode .EQ. 1) THEN
    CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                              &
                   TRIM('### 3.) Check if Norm of spectral function equals 1.0'), &
                   'append',debug_maxent)
  ELSE IF (kernelMode .EQ. 10 .OR. kernelMode .EQ. 11) THEN 
    CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                                                   &
                   TRIM('### 3.) Check model contraints (evaulated from spectral function):'), &
                   'append',debug_maxent)
    CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM('###'),'append',debug_maxent)
                   
    CALL chiMoments (specFunc, specGrid, specUnit, NSpec, kernelBeta, kernelMode, &
                     chiOmegaZero, chiTauZero)
                   
    WRITE(logMsg,"(A,F15.5)") '###      REALPART[chi(omega=0)] = ', chiOmegaZero
    CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM(logMsg),'append',debug_maxent)
    
    WRITE(logMsg,"(A,F15.5)") '###               chi(tau  =0)  = ', chiTauZero
    CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),TRIM(logMsg),'append',debug_maxent)
  END IF
  !*** calculate moments of susceptibility, Eq. 5.15 and 5.16 in Ref[1] for kernel 10 and 11
  
  CALL writeLog (TRIM(logPrefix)//TRIM(logMaxent),                                                  &
       TRIM('####################################################################################'),&
       'append',debug_maxent)

  CALL writeLog (TRIM(logPrefix)//TRIM(logBayes),                                                   &
       TRIM('######################################################################'),&
       'append',debug_bayes)

  CALL writeLog (TRIM(logPrefix)//TRIM(logMontecarlo),                                              &
       TRIM('#################################################################################'),&
       'append',debug_mc)

  
  !*** logging user instructions 
  !***************************************************************************************

  !*** normalizing spectral function 
  CALL normalizeToGrid(specFunc,specGrid,NSpec,kernelMode)
  !*** normalizing spectral function 

  DEALLOCATE(maxentCorrFunc,maxentInversCov,kernel)
  
END SUBROUTINE w2maxent
!*****************************************************************************************
!*****************************************************************************************



!*****************************************************************************************
!*** initKernel **************************************************************************
!*****************************************************************************************
!***
!*** Initializes kernel(Ngreen,NSpec) which defines the dependency between Green's function
!*** and spectral function:
!*** G(i)=Integral( Kernel(i,w)*specFunc(w) dw )
!*** depending on kernelMode, i may be imaginary time (kernelMode 0,10,11) or Matsubara frequencies (kernelMode 1)
!*** 
!*** In the discrete case (finite number of measured grid points), the integral becomes a sum
!*** G(i)=Sum( Kernel(i,j)*specFunc(j)*dj ), where dj = specUnit(j) is the spectral unit.
!***   unlike Eq. 3.9 (and below) in Ref[1] we absorb dj into the kernel (method by Sandvik): 
!***
!***   kernel(i,j)=Kernel(i,j)*dj
!***
!***   This simplifies the normalization procedures of the spectral function in the
!***   calculations (but the units have to be considered in SUBROUTINE initLikelihoodCurvature and bayesianInverence).
!***   There is no considerable enhancement in processing time otherwise.
!***
!*** PARAMETERS:
!*** specGrid defines the spectral grid.
!*** specUnit defines the spectral weight for numerical integration.
!*** Nspec defines the number of spectral grid points.
!*** greenFunc is the complex valued data-vector
!***           for kernel-modes 0,10,11 only real part of greenFunc is considered in maxentCorrFunc (dimension NGreen)
!***           for kernel-mode 1 real and imaginary part are considered seperately in maxentCorrFunc (extendet to dimension NGreen*2)
!*** greenInverseCov is the inverse covariance matrix of the data (considered to be diagonal)
!***           for kernel-modes 0,10,11 maxentInverseCorr(:)=greenInverseCov(:)
!***           for kernel-modes 1 maxentInverseCorr(1:NGreen)=greenInverseCov(:) and maxentInverseCorr(NGreen+1:NGreen*2)=greenInverseCov(:)
!*** greenGrid defines the grid where the Green's function is defined
!*** Ngreen defines the number of measurements
!*** kerneMode defines the actual functional form of the kernel
!*** kernelBeta defines the inverse temperature 1/(kb*T), hence the temperature used by evaluating the Green's function
!*** kernel ... defines the relation between Green's function and spectral function: G(i)=Integral( Kernel(i,w)*specFunc(w) dw )
!*** maxentCorrFunc is an allocated output parameter: the data-vector which is used in calculation (dimension of vector depends on kernelMode)
!*** maxentInverseCorr is an allocated output parameter: the inverse covariance-matrix (vector, diagonal) which is used in calculation (dimension of vector depends on kernelMode)
!*** maxentNCorr is the dimension of the data- and covariance-vector (diagonal) depending on kernelMode
!***
!*****************************************************************************************
SUBROUTINE initKernel(specGrid,specUnit,NSpec,                                     &
                      greenFunc,greenInverseCov,greenGrid,NGreen,                  &
                      kernelMode,kernelBeta,maxentNCorr)

!(kernel,specGrid,specUnit,greenGrid,NSpec,Ngreen, &
!                      kernelBeta,kernelMode)
  IMPLICIT NONE
  
  !***************************************************************************************
  !*** Spectral Function Variables ***
  INTEGER,                        INTENT(IN)  :: NSpec
  REAL(KINDR),DIMENSION(Nspec),        INTENT(IN)  :: specGrid
  REAL(KINDR),DIMENSION(Nspec),        INTENT(IN)  :: specUnit
  !*** Spectral Function Variables ***

  !*** Green's Function Variables ***
  INTEGER,                         INTENT(IN)  :: NGreen
  COMPLEX(KINDC),DIMENSION(Ngreen),     INTENT(IN) :: greenFunc
  REAL(KINDR),DIMENSION(Ngreen),        INTENT(IN)  :: greenGrid
  REAL(KINDR),DIMENSION(Ngreen),        INTENT(IN)  :: greenInverseCov
  !*** Green's Function Variables ***
  
  !*** Kernel Variables ***
  INTEGER,                           INTENT(OUT) :: maxentNCorr
  REAL(KINDR),                            INTENT(IN)  :: kernelBeta
  INTEGER,                           INTENT(IN)  :: kernelMode
  !*** Kernel Variables ***

  !*** Local Helpers ***
  INTEGER                                     :: i,j
  !*** Local Helpers ***
  !***************************************************************************************

  IF (kernelMode .EQ. 1) THEN     ! Fermionic Matsubara kernel --> extend dimension of correlation function
    maxentNCorr=2*NGreen          !                                to treat real & imag part seperately
  ELSE 
    maxentNCorr=NGreen            ! Other kernels, dimension is equal to dimension of initial correlation function
  END IF

  ALLOCATE(maxentCorrFunc (maxentNCorr)      )    ! allocate working array (vector) for correlation function
  ALLOCATE(maxentInversCov(maxentNCorr)      )    ! allocate working array (vector) for inverse covariance matrix
  ALLOCATE(kernel         (maxentNCorr,NSpec))    ! allocate working array (matrix) for the kernel

  !***************************************************************************************
  !*** init kernel
  kernel              (:,:)=0D0

  IF (kernelMode .EQ. 1) THEN                                              ! initialize kernel for kernelMode=1 - Fermionic Matsubara Kernel
    maxentCorrFunc        (       1:NGreen     )= REAL(greenFunc(:))       ! Real part of Green's function
    maxentCorrFunc        (NGreen+1:maxentNCorr)=AIMAG(greenFunc(:))       ! Imaginary part of Green's function

    maxentInversCov       (       1:NGreen     )=      greenInverseCov(:)  ! Covariance matrix is also extended where it is assumed that 
    maxentInversCov       (NGreen+1:maxentNCorr)=      greenInverseCov(:)  ! real and imaginary part of Green's function have same error !
  ELSE 
    maxentCorrFunc        (:)=REAL(greenFunc(:))                           ! initialize kernel for other kernelModes /= 1
    maxentInversCov       (:)=greenInverseCov(:)
  END IF
  
  !*** init kernel
  !***************************************************************************************

  !***************************************************************************************
  IF (kernelMode == 0) THEN
  !*** kernel 0 --> imaginary time green's function ***
    DO j=1,NSpec
      DO i=1,NGreen

        IF (specGrid(j)*greenGrid(i) .LE. -kernelExponentCutoff) THEN ! limit exp(-tau*omega) >> 1 -> Kernel ~ exp ((beta-tau)*omega)
          IF ((kernelBeta-greenGrid(i))*specGrid(j) .GT. -kernelExponentCutoff) THEN ! omega < 0, (beta-tau)*omega <= 0 
            kernel(i,j)=EXP((kernelBeta-greenGrid(i))*specGrid(j))
          ! ELSE asymptotic zero
          END IF
        ELSE IF (specGrid(j)*greenGrid(i) .LT. kernelExponentCutoff) THEN
          kernel(i,j)=EXP(-greenGrid(i)*specGrid(j))/(1D0+EXP(-kernelBeta*specGrid(j))) 
!        ELSE: 0D0
        END IF
      END DO
      
      ! method by Sandvik, spectral weight added to kernel, 
      kernel(:,j)=kernel(:,j)*specUnit(j)
      
    END DO
  !*** kernel 0 --> imaginary time green's function ***
  ELSE IF (kernelMode == 1) THEN 
  !*** kernel 1 --> Fermionic Matsubara Kernel (watch out for singular values of kernel in Bosonic case !) ***
    DO i=1, NGreen
      DO j=1,NSpec
      ! REAL PART OF KERNEL
        kernel(       i,j) = -1D0* specGrid(j)/(greenGrid(i)**2+specGrid(j)**2)*specUnit(j)
      ! IMAG PART OF KERNEL
        kernel(NGreen+i,j) = -1D0*greenGrid(i)/(greenGrid(i)**2+specGrid(j)**2)*specUnit(j)
      END DO
    END DO
  !*** kernel 1 --> Fermionic Matsubara Kernel (watch out for singular values of kernel in Bosonic case !) ***
  ELSE IF (kernelMode == 10) THEN
  !*** kernel 10 --> bosonic kernel for susceptibility ***
    DO j=1,NSpec
      WHERE(greenGrid(:)*specGrid(j).LT.kernelExponentCutoff)               ! neglect DEXP(-tau*omega)<<1, omega >= 0
        kernel(:,j)=EXP(-greenGrid(:)*specGrid(j))
      END WHERE 

      WHERE((kernelBeta-greenGrid(:))*specGrid(j).LT.kernelExponentCutoff)  ! neglect DEXP(-(beta-tau)*omega)<<1
        kernel(:,j)=kernel(:,j) + EXP(-(kernelBeta-greenGrid(:))*specGrid(j))
      END WHERE 

      IF (specGrid(j) .EQ. 0D0) THEN
        kernel(:,j)=2D0/kernelBeta !limit omega -> 0: K(tau,omega)-> 2/beta
      ELSE
        kernel(:,j)=specGrid(j)*kernel(:,j)/(1D0-EXP(-kernelBeta*specGrid(j)))
      END IF

      ! method by Sandvik, spectral weight added to kernel, 
      kernel(:,j)=kernel(:,j)*specUnit(j)

    END DO
    
  !*** kernel 10 --> bosonic kernel for susceptibility ***
  ELSE IF (kernelMode == 11) THEN
  !*** kernel = 11 --> bosonic spectral function -> 1/(1-exp(-beta*omega)) in spectral function ***
    DO j=1,NSpec
      WHERE(greenGrid(:)*specGrid(j).LT.kernelExponentCutoff)               ! neglect DEXP(-tau*omega)<<1
        kernel(:,j)=EXP(-greenGrid(:)*specGrid(j))
      END WHERE 

      WHERE((kernelBeta-greenGrid(:))*specGrid(j).LT.kernelExponentCutoff)  ! neglect DEXP(-(beta-tau)*omega)<<1
        kernel(:,j)=kernel(:,j) + EXP(-(kernelBeta-greenGrid(:))*specGrid(j))
      END WHERE 

      ! method by Sandvik, spectral weight added to kernel, 
      kernel(:,j)=kernel(:,j)*specUnit(j)
    END DO
  !*** kernel = 11 --> bosonic spectral function -> 1/(1-exp(-beta*omega)) in spectral function ***
  ELSE
    WRITE(outputUnit,*) 'ERROR: kernel mode ',kernelMode,' not (yet ?!) implemented ... maxent stop'
    DEALLOCATE(maxentCorrFunc,maxentInversCov,kernel)
    STOP
  END IF
  !***************************************************************************************
  
END SUBROUTINE initKernel
!*****************************************************************************************
!*****************************************************************************************



!*****************************************************************************************
!*** initLikelihoodCurvature *************************************************************
!*****************************************************************************************
!*** 
!*** Here the curvature-tensor of the loglikelihood function w.r.t. the spectral
!*** function, A(i), is calculated:
!***   d² L(i,j) / (dA(i)dA(j))
!***
!*** The quantity enters in the determination of alpha in the subroutine 
!*** BayesianInverence(..) via Eq. 4.8 in Ref[1]:
!***   Lambda_ij ~ DSQRT(A(i)) * d² L(i,j) / (dA(i)dA(j)) * DSQRT(A(j))
!***
!*** The curvature-tensor won't change in calculation (see Eq. 4.9 in Ref[1]):
!***   d² L(i,j) / (dA(i)dA(j)) = SUM_kl( kernel(k,i)*greenInverseCov(k,l)*kernel(l,j) )
!*** 
!*** PARAMETERS:
!***  likelihoodCurvature is the curvature tensor d² L(i,j) / (dA(i)dA(j)).
!***  specUnit is the spectral grid-weight for numerical integration
!***  kernel defines the dependency between spectral function and the Green's function 
!***         and should be of the form SUBROUTINE initKernel(..).
!***  maxentInverseCov is the inverse covariance matrix of the measured Green's function.
!***  NSpec is the number of spectral grid points.
!***  NGreen is the number of Green's function grid points.
!*****************************************************************************************
SUBROUTINE initLikelihoodCurvature(likelihoodCurvature,specUnit,NSpec,NGreen)
  IMPLICIT NONE

  !*** Spectral Function Variables ***
  INTEGER,                       INTENT(IN)  :: NSpec           ! is the number of spectral grid points.
  REAL(KINDR),DIMENSION(Nspec)                    :: specUnit        !  is the spectral grid-weight for numerical integration
  !*** Spectral Function Variables ***

  !*** Green's Function Variables ***
  INTEGER,                       INTENT(IN)  :: NGreen           ! is the number of Green's function grid points.
  !*** Green's Function Variables ***

  !*** Maximum Entropy (Bayesian) Parameter ***
  REAL(KINDR),DIMENSION(Nspec,Nspec), INTENT(OUT) :: likelihoodCurvature ! is the curvature tensor d² L(i,j) / (dA(i)dA(j)).
  !*** Maximum Entropy (Bayesian) Parameter ***

  !*** Local Helpers ***
  INTEGER                                     :: i,j
  !*** Local Helpers ***

  !method by Sandvik, spectral weight removed from log likelihood curvature
  DO i=1,NSpec
    DO j=1,NSpec
      likelihoodCurvature(i,j)=SUM(kernel(:,i)*maxentInversCov(:)*kernel(:,j))   ! see Eq. 4.9 in Ref[1]
      likelihoodCurvature(i,j)=likelihoodCurvature(i,j)/specUnit(i)/specUnit(j)      
    END DO                  
  END DO                                                                       

END SUBROUTINE initLikelihoodCurvature
!*****************************************************************************************
!*****************************************************************************************



!*****************************************************************************************
!*** bayesianInverence *******************************************************************
!*****************************************************************************************
!***
!*** Here the Bayesian convergence criterion is checked, Eq. 4.28 in Ref[1]:
!***   -2alpha*S = Trace(Lambda_ij*(alpha*1_ij + Lambda_ij)^-1)
!*** for a given alpha, a given spectral Function and a given likelihoodCurvature (hence 
!*** a given kernel and covariance matrix of the Green's function).
!*** S is the entropy and the matrix Lambda is defined in by Eq. 4.8 in Ref[1].
!*** 
!*** If the convergence (up to a numerical factor controlled by MODULE parameter bayesianConvergence)
!***   DABS(-2alpha*S / Trace(Lambda_ij*(alpha*1_ij + Lambda_ij)^-1) - 1D0) < bayesianConvergence
!*** is reached the parameter converged = .TRUE., 
!*** else converged = .FALSE. and the Bayesian parameter alpha will be updated iff converged=.FALSE. was the input-value, 
!***   in the main loop a new spectral function should be calculated with the new alpha.
!***
!*** PARAMETER:
!***   bayesTrace is an output parameter: Trace(Lambda_ij*(alpha*1_ij + Lambda_ij)^-1)
!***   bayesPostExponent is the convergence-criterion counter part of bayesTrace
!***     -2alpha*S
!***   converged determines if Eq. 4.28 in Ref[1] is fulfilled (up to a numerical factor) 
!***   alpha is the Bayesian Parameter which enters in the entropic prior, DEXP(alpha*S)
!***     and is updated if no convergence is reached and iff input value of converged=.FALSE. (to avoid endless loop in the smoothing procedure)
!***   likelihoodCurvature is the curvature tensor of the log-likelihood function (Eq. 4.9 in Ref[1])
!***   specFunc is the current spectral function which should fulfill the convergence-critetion
!***   specModel is the spectral model which enters in the entropy
!***   specUnit is the spectral grid weight for numerical integration of the entropy
!***   NSpec is the number of spectral grid points
!***   logName is the file-name, where the current covnergence properties will be stored
!***      (logIndex, alpha, entropy, bayesTrace, bayesPostExponent/bayesTrace)
!***   logIndex is written in the first column of the logFile and is afterward incremented by 1
!***   chiSquared is the chi-squared value obtained in the MonteCarlo run for a specific alpha and is meant for logging
!***   logState is a string meant for loggingwhich determines (helps the user to see) the state 
!***            of the maxent-procedure (warmup,main-loop=converge,smoothing)
!*****************************************************************************************
SUBROUTINE bayesianInverence (bayesTrace, bayesPostExponent, converged, alpha, &
                              likelihoodCurvature, specFunc, specModel,specUnit, NSpec, &
                              logName, logIndex,chiSquared,logState)
  USE MMersenneTwister
  IMPLICIT NONE

  !***************************************************************************************
  !*** Spectral Function Variables ***
  INTEGER,                      INTENT(IN)                 :: NSpec        ! is the number of spectral grid points
  REAL(KINDR),  DIMENSION(Nspec),    INTENT(IN)  :: specFunc     ! is the current spectral function which should fulfill the convergence-critetion
  REAL(KINDR),  DIMENSION(Nspec),    INTENT(IN)  :: specModel    ! is the spectral model which enters in the entropy
  REAL(KINDR),  DIMENSION(Nspec),    INTENT(IN)  :: specUnit     ! is the spacing between the grid points
  !*** Spectral Function Variables ***

  !*** Maximum Entropy (Bayesian) Parameter ***
  REAL(KINDR)                                    :: alpha        ! enters in the entropic prior DEXP(alpha*entropy) and is updated, if no convergence is reached
  REAL(KINDR),                       INTENT(IN)  :: chiSquared   ! chiSquared value for logging
  REAL(KINDR),DIMENSION(Nspec,Nspec),INTENT(IN)  :: likelihoodCurvature ! the curvaturetensor of the DLOG-likelihood function w.r.t. specFunc
  REAL(KINDR),DIMENSION(Nspec,Nspec)             :: lambdaMatrix ! enters the convergence criterion, defined in Eq. 4.8 in Ref[1]
  REAL(KINDR),DIMENSION(Nspec,Nspec)             :: gammaMatrix  ! is the curvature-tensor of the DLOG-posteriori pdf, Q=alpha*S-chi²/2 and also enters the
                                                            !    convergence criterion     -2alpha*S = Trace( lambdaMatrix_ij * (gammaMatrix_ij)^-1 )
                                                            !    gammaMatrix = DSQRT(Ai) * d²Q/(dAi*dAj) * DSQRT(Aj) = alpha*1_ij + Lambda_ij
                                                            !    see Eq. 4.4 and the definition of Gamma_ij above Eq. 4.8 in Ref[1]
  REAL(KINDR)                                    :: entropy      ! the entropy of the spectral function w.r.t. the model, also referenced as S 
  REAL(KINDR),                       INTENT(OUT) :: bayesTrace   ! Trace( lambdaMatrix_ij * (gammaMatrix_ij)^-1 )
  REAL(KINDR),                       INTENT(OUT) :: bayesPostExponent !-2alpha*Entropy
  LOGICAL                                   :: converged    ! true if Bayesian convergence is reached
  REAL(KINDR)                                    :: bayesRate    ! local helper: bayesTrace/bayesPostExponent
  !*** Maximum Entropy (Bayesian) Parameter ***

  !*** DLOG Variables ***
  CHARACTER(*),                 INTENT(IN)  :: logName         ! filename of the logFile
  CHARACTER(*),                 INTENT(IN)  :: logState        ! 
  INTEGER                                   :: logIndex        ! index / first column of a logFile entry, is incremented by 1
  CHARACTER(len=128)                            :: logMsg          ! logIndex, alpha, entropy, bayesTrace, bayesPostExponent/bayesTrace
  !*** DLOG Variables ***

  !*** Local Helpers ***
  INTEGER                                   :: i,j          ! loop helpers
  INTEGER                                   :: info         ! helper for LAPACK invoke
  INTEGER,  DIMENSION(NSpec)                :: ipiv         ! helper for LAPACK invoke
  REAL(KINDR),  DIMENSION(NSpec*64)              :: work         ! helper for LAPACK invoke
  !*** Local Helpers ***
  !***************************************************************************************

  !***************************************************************************************
  !*** part 1 --> calculate bayesTrace: Trace( lambdaMatrix_ij * (gammaMatrix_ij)^-1 )

  DO i=1,NSpec  
    ! lambdaMatrix_ij is calculated here : see Eq. 4.8 in Ref[1]
    DO j=1, NSpec
      lambdaMatrix(i,j)= SQRT(specFunc(i)*specUnit(i)) &
                        *likelihoodCurvature(i,j)           &
                        *SQRT(specFunc(j)*specUnit(j))
      gammaMatrix (i,j)=lambdaMatrix(i,j)
    END DO

    ! GAMMAMatrix_ij is calculated here : see Eq. 4.4 and above Eq. 4.8 in Ref[1]
    gammaMatrix(i,i)=gammaMatrix(i,i)+alpha
  END DO     
  
  CALL DGETRF(NSpec,NSpec,gammaMatrix,NSpec,ipiv,info)         ! prepare gammaMatrix for inversion
  CALL DGETRI(NSpec,gammaMatrix,NSpec,ipiv,work,NSpec*64,info) ! Invert gammaMatrix

  ! calculate trace 
  bayesTrace=0D0
  DO i=1,Nspec
    DO j=1,Nspec
      bayesTrace=bayesTrace+lambdaMatrix(i,j)*gammaMatrix(j,i)
    END DO
  END DO

  !  this trace is a measure for the number of good data --> Ng, see Eq. 4.29 in Ref[1]
  !  -> DIAGONALIZE lambda, eigenvalues lambda_i 
  !     Trace (lambda_ij * (alpha*1_ij + lambda_ij)^-1) -> SUM_i lambda_i/(alpha+lambda_i) = Ng
  !     for lambda_i >> alpha --> lambda_i/(alpha+lambda_i) = 1 ... good data point
  !     for lambda_i << alpha --> lambda_i/(alpha+lambda_i) ~ 0 assumed ... not good data
  
  !*** part 1 --> calculate number of good data Ng=SUM_i lambda_i/(alpha+lambda_i)
  !***************************************************************************************
  
  !***************************************************************************************
  !*** part 2 --> calculate Bayesian posteriori DEXPonent: DEXP^(-2D0*alpha*entopy))
  entropy =  -SUM(  specFunc(:)*DLOG(specFunc(:)/specModel(:))*specUnit(:) &
              ,1,specFunc(:).GT.accuracy.AND.specModel.GT.accuracy)

  ! include integral over [A(w) - m(w) ] if m(w) is not normalized ...
  entropy =  entropy + SUM( (specFunc(:)-specModel(:))*specUnit(:) )
            
  bayesPostExponent=-2D0*alpha*entropy
  !*** part 2 --> calculate Bayesian posteriori exponent: DEXP^(-2D0*alpha*entopy))
  !***************************************************************************************

  !***************************************************************************************
  !*** part 3 --> check convergence
  bayesRate=bayesTrace/bayesPostExponent

  !*** logging
  WRITE(logMsg,"(I0.5,F15.7,F15.4,F15.4,F20.10)") logIndex, alpha, entropy, bayesTrace, &
                                                       1D0/bayesRate
  CALL writeLog (TRIM(logName)//TRIM(logBayes),logMsg,'append',debug_bayes)    
  logIndex=logIndex+1  


  WRITE(logMsg,"(F17.3,F12.5,F10.5,F16.5,F13.5,A10,F6.3)") chiSquared,entropy,bayesTrace,&
                             1D0/bayesRate,alpha, TRIM(logState),SUM(specFunc(:)*specUnit(:))
  CALL writeLog (TRIM(logName)//TRIM(logMaxent),logMsg,'append',debug_maxent)    
  logIndex=logIndex+1  
  !*** logging
  
  IF (ABS(1D0/bayesRate-1D0) .GE. bayesianConvergence) THEN 
    IF (converged .EQV. .FALSE.) THEN
      !*** variation of alpha from Sandvik
      IF (bayesRate .LT. bayesianConvergence) THEN          ! very big differences: bayesTrace/bayesPostExponent << 1
        alpha=alpha*bayesianUpdateAlpha                     !   alpha should in general become much smaller ... e.g. alpha*0.05
      ELSE                                                  ! differences approach
        alpha=alpha*bayesRate*(1.+1d-3*(grnd()-0.5))        !
      END IF
      !*** variation of alpha from Sandvik
    END IF
    
    converged=.FALSE.
  ELSE 
    converged=.TRUE.                                     ! convergence reached for current alpha and specFunc
  END IF
  !*** part 3 --> check convergence
  !***************************************************************************************
  
  

END SUBROUTINE bayesianInverence
!*****************************************************************************************
!*****************************************************************************************



!*****************************************************************************************
!*** montecarloSpectralFunction **********************************************************
!*****************************************************************************************
!*** 
!*** Here the most probable spectral function w.r.t. to a fixed model and a fixed Bayesian
!*** parameter alpha is obtained by a simmulated annealing Monte Carlo algorithm.
!***
!*** The probability weighting of keeping a new configuration is based on the log-posteriori 
!*** probability density: Q = alpha*S - chi^2/2 (Eq. 3.22 in Ref[1]) where S is the entropy 
!*** and chi^2 is the sum of gaussian residuals between the "real" value of the Green's function, 
!*** G_i, (calculated with the kernel) and the measured Green's function, greenFunc.
!***   S     = - ingetral (specFunc(w) - specModel(w) - specFunc(w)*DLOG(specFunc(w)/specModel(w)) dw)
!***   chi^2 = SUM_ij [ (G_i - greenFunc(i)) * greenInverseCov_ij * (G_j - greenFunc(j)) ]
!***   where greenInverseCov is the inverse covariance matrix of the measurements (and is 
!***   assumed to be diagonal)
!*** 
!*** PARAMETER:
!***   specFunc is the spectral function, which is to be optimized
!***   specModel is the spectral model entering in the entropic prior
!***   specUnit is the spectral weight for numerical integration
!***   kernel defines the relationship between the Green's function and the spectral function
!***   greenFunc are the measured values for the Green's function, data
!***   greenInverseCov is the inverse covariance matrix of the measurements 
!***     assumed to be diagonal, no covariances between measurements
!***   alpha defines the fixed Bayesian parameter entering the entropic prior DEXP(alpha*S)
!***   NSpec is the number of spectral grid points
!***   NGreen is the number of measurements of the Green's function
!***   NMonteCarlo is number of Monte Carlo runs where new configurations of specFunc are configured
!***   NAnnealing is the number of MonteCarlo runs until a new annealing temperature is defined 
!***                     number of MonteCarlo runs until the system is cooled down
!***   annealingInit is the initial annealing parameter applied on the full posteriori pdf
!***     DEXP(Q/annealingInit) = DEXP[(alpha*S - chi^2/2)/annealingInit] 
!***     see Ref[3] for more details
!***   variationInit is a measure, how much the magnitude of the spectral function can 
!***     change in one MonteCarlo step
!***     big   variationInit -> big   variations are allowed
!***     small variationInit -> small variations are allowed
!***     this parameter is adapded afer every annealing step
!***   chiSquared is the sum of gaussian residuals between the "real" value of the Green's function, 
!***     G_i, (calculated with the kernel) and the measured Green's function, greenFunc.
!***   logName is the filename of the logFile
!***   logIndex is a running index meant as first column in the logFile
!****************************************************************************************
SUBROUTINE montecarloSpectralFunction (specFunc, specModel, specUnit,           &
                            alpha, NSpec, NGreen,                               &
                            NMonteCarlo,NAnnealing, annealingInit,variationInit,&
                            chiSquared,logName,logIndex)
  USE MMersenneTwister
  IMPLICIT NONE

  !**************************************************************************************
  !*** Spectral Function Variables ***
  INTEGER,                   INTENT(IN)  :: NSpec         ! number of spectral grid points
  REAL(KINDR),  DIMENSION(Nspec), INTENT(OUT) :: specFunc      ! the spectral function, which is to be optimized
  REAL(KINDR),  DIMENSION(Nspec), INTENT(IN)  :: specModel     ! the spectral model entering in the entropic prior
  REAL(KINDR),  DIMENSION(Nspec), INTENT(IN)  :: specUnit      ! the spacing between spectral grid points
  REAL(KINDR),  DIMENSION(Nspec)              :: specVariation ! allowed variation magnitude of the spectral function (for every grid point) 
                                                          !   for the Monte Carlo runs in a specific annealing configuration.
                                                          !   kind of importance sampling, if you want so ... from former version
  !*** Spectral Function Variables ***

  !*** Green's Function Variables ***
  INTEGER,                   INTENT(IN)  :: NGreen          ! number of measurements of the Green's function
  REAL(KINDR),  DIMENSION(NGreen)             :: greenCalc       ! the "real" values of the Green's function given the kernel and specFunc
  !*** Green's Function Variables ***

  !*** Kernel Variables ***

  !*** Bayesian Variables ***
  REAL(KINDR),                    INTENT(OUT) :: chiSquared    ! DLOG likelihood function, mean of squared deviations of greenFunc from greenCalc
  REAL(KINDR),                    INTENT(IN)  :: alpha         ! fixed Bayesian parameter entering the entropic prior
  !*** Bayesian Variables ***

  !*** Montecarlo Variables ***
  INTEGER,                   INTENT(IN)  :: NMonteCarlo   ! number of Monte Carlo runs for finding new configuration of specFunc
  INTEGER,                   INTENT(IN)  :: NAnnealing    ! number of MonteCarlo runs until a new annealing temperature is defined 
  REAL(KINDR)                                 :: annealingInit, &   ! initial value for the annealing parameter
                                            annealingBoltzmann ! actual annealing parameter: DEXP(Q/annealingBoltzmann)
  REAL(KINDR)                                 :: variationInit, &   ! initial parameter of magnitude of variation of the spectral function
                                            variationFact      ! actual variational parameter used in calculation
  REAL(KINDR)                                 :: try,    &     ! number of sweeps in a reagion of specFunc, where the variation is performed
                                                          !   above the noiseLevel
                                            accept, &     ! nubmer of accepted sweeps above the noiseLevel
                                            noiseLevel    ! defines the magnitude of specFunc, where a sweep is considered to be important
                                                          !   in combination with variationFact --> kind of importance sampling
  INTEGER                                :: sweepI,sweepJ ! indices where specFunc is variated in a single sweep
  REAL(KINDR)                                 :: sweepJValue,sweepIValue ! values of specFunc after variation at sweepI and sweepJ
  REAL(KINDR),  DIMENSION(NGreen)             :: sweepGreenCalc          ! changed "real" value of the Green's function because of sweep
  REAL(KINDR)                                 :: sweepChiSquared         ! changed DLOG-likelihood function because of sweep
  REAL(KINDR)                                 :: sweepEntropy            ! change in entropy because of sweep
  REAL(KINDR)                                 :: sweepWeight             ! probability ratio before/after sweep
  !*** Montecarlo Variables ***

  !*** DLOG Variables ***
  CHARACTER(*),              INTENT(IN)  :: logName       ! filename of the logFile
  INTEGER                                :: logIndex      ! index (first column) of DLOG entry
  CHARACTER(len=128)                         :: logMsg        ! message to be DLOGged
  !*** DLOG Variables ***
  
  !*** Local Helpers ***
  INTEGER                                :: i,j           ! loop helpers
  !*** Local Helpers ***
  !**************************************************************************************

  !**************************************************************************************
  !*** initializing
  annealingBoltzmann = annealingInit
  variationFact      = variationInit
  !**************************************************************************************


  !**************************************************************************************
  !*** initializing the annealing configuration
  try=0D0
  accept=0D0

  !noiseLevel is defined with respect to the global maximum of the current spectral function
  noiseLevel=noiseFilter*MAXVAL(specFunc)   
  
  !variation from Sandvik
  !the spectral variation is taken from the former version
  specVariation(:)=variationFact*(noiseFilter*specFunc(:)+noiseFilter*noiseLevel)
  !*** initializing the annealing configuration
  !**************************************************************************************

  !**************************************************************************************
  !*** MONTE CARLO LOOP *****************************************************************
  DO i=1,NMonteCarlo  

    !calculate "real" values for G(tau)
    DO j=1,NGreen
      greenCalc(j)=SUM(kernel(j,:)*specFunc(:)) !spectralGridUnit in kernel  
    ENDDO

    !calculate chi^2 value
    chiSquared=SUM(maxentInversCov(:)*(maxentCorrFunc(:)-greenCalc(:))**2)

    !******************************************************************
    !*** MONTE CARLO SWEEPS
    DO j=1,NSpec                        !find new configuration, if sweepI and sweepJ match
      !****************************************************************
      !*** RANDOM SWEEP POSITION AND CHANGE IN MAGNITUDE
      DO     ! loop: code from Sandvik
        sweepI=CEILING(grnd()*NSpec) ! changed by bhartl, now NSpec is also variated
        sweepJ=CEILING(grnd()*NSpec) ! changed by bhartl, now NSpec is also variated
        IF(sweepI .EQ. 0 .OR. sweepJ .EQ. 0 .OR. &  ! changed by bhartl: MMersenneTwister prepares random numebers in interval [0,1], so 0 could be possible
           sweepI.EQ.sweepJ) THEN                   ! sweepI and sweepJ must not be equal
           CYCLE                                    ! find new configuration ...
        END IF                   

        sweepIValue=specVariation(sweepI)*(grnd()-0.5)                   ! variate specFunc at position sweepI
        IF(specFunc(sweepI)+sweepIValue.LT.0D0) CYCLE                    ! specFunc must be semi positiv

        sweepJValue = -sweepIValue*specUnit(sweepI)/specUnit(sweepJ) &   ! [...]*specUnit(speepI)/specUnit(sweepJ) conserves spectral weight
               + specVariation(sweepJ)*0.5D0*noiseFilter*(grnd()-0.5)    ! NORMALIZATION violated. see !*** SWEEP CHANGED ENTROPY
        IF ((specFunc(sweepJ)+sweepJValue).GE.0D0) EXIT                  ! find total new configuration if negative spectral function
                                                                         ! otherwise endless loop may occur, changed by bhartl: .GT. -> .GE.
      END DO ! loop: code from Sandvik
      !*** RANDOM SWEEP POSITION AND CHANGE IN MAGNITUDE
      !****************************************************************

      !****************************************************************
      !*** SWEEP CHANGED "REAL" GREEN'S FUCNTION AND CHI-SQUARED VALUE
      IF (specFunc(sweepI).GT.noiseLevel) THEN !sweep in (concidered) important area of specFunc
        try=try+1D0
      END IF
      
      !calculate changed "real" Green's function given the spectral function after the sweep 
      sweepGreenCalc(:)=greenCalc(:)+kernel(:,sweepI)*sweepIValue+kernel(:,sweepJ)*sweepJValue 
      
      !calculate changed in the DL-likelihood function
      sweepChiSquared=SUM(maxentInversCov(:)*(maxentCorrFunc(:)-sweepGreenCalc(:))**2)
      !*** SWEEP CHANGED "REAL" GREEN'S FUCNTION AND CHI-SQUARED VALUE
      !****************************************************************
      
      !****************************************************************
      !*** SWEEP CHANGED ENTROPY: [Entropy(post sweep) - Entropy(pre sweep)] 
      !***                        at two positions sweep: sweepI, sweepJ

      sweepEntropy = sweepIValue*specUnit(sweepI) + & ! changed by bhartl, spectrum may not be
                     sweepJValue*specUnit(sweepJ)     ! normalized any more because of sweepJValue: see
                                                      ! *** RANDOM SWEEP POSITION AND CHANGE IN MAGNITUDE
!      sweepEntropy = 0D0

      ! method by sandvik
      ! old entropy at index sweepI                    
      IF (specFunc(sweepI).GT.accuracy)                                       THEN
        sweepEntropy = sweepEntropy                                              &
                     + ( specUnit(sweepI)*specFunc(sweepI)                       &
                         *LOG(specFunc(sweepI)/specModel(sweepI))               &
                        ) 
      END IF 
      
      ! sweeped entropy at index sweepI                    
      IF ((specFunc(sweepI)+sweepIValue).GT.accuracy)                         THEN
        sweepEntropy = sweepEntropy                                              &
                     - ( specUnit(sweepI)*(specFunc(sweepI)+sweepIValue)         &
                         *LOG((specFunc(sweepI)+sweepIValue)/specModel(sweepI)) &
                        ) 
      END IF
      
      ! old entropy at index sweepJ                    
      IF (specFunc(sweepJ).GT.accuracy)                                       THEN
        sweepEntropy = sweepEntropy                                              &
                     + ( specUnit(sweepJ)*specFunc(sweepJ)                       &
                         *LOG(specFunc(sweepJ)/specModel(sweepJ))               &
                        )
      END IF
     
      ! sweeped entropy at index sweepJ                    
      IF ((specFunc(sweepJ)+sweepJValue).GT.accuracy)                         THEN
        sweepEntropy = sweepEntropy                                              &
                     - ( specUnit(sweepJ)*(specFunc(sweepJ)+sweepJValue)         &
                         *LOG((specFunc(sweepJ)+sweepJValue)/specModel(sweepJ)) &
                        )
      END IF
      !*** SWEEP CHANGED ENTROPY
      !****************************************************************


      !****************************************************************
      !*** APPLY METROPOLIS TO SWEEP CHANGED PROBABILITY DEXPONENT
      ! depends only on change in chi^2 value and change in entropy
      ! probabilty weight:  
      ! exp ( -chi_sweeped^2 /2 + alpha*S_sweeped) / exp( -chi^2 / 2 + alpha*S)
      sweepWeight=((chiSquared-sweepChiSquared)/2D0+alpha*sweepEntropy)/annealingBoltzmann
      sweepWeight=EXP(MIN(sweepWeight,0D0))

      IF(grnd().LT.sweepWeight)THEN ! Metropolis !

        !accept new configuration
        IF (specFunc(sweepI).GT.noiseLevel) THEN
          !change in important are of specFunc
          accept=accept+1D0 
        END IF

        !keep new configuration
        specFunc(sweepI)=specFunc(sweepI)+sweepIValue
        specFunc(sweepJ)=specFunc(sweepJ)+sweepJValue
        greenCalc(:)=sweepGreenCalc(:)
        chiSquared=sweepChiSquared
      END IF
      !*** APPLY METROPOLIS TO SWEEP CHANGED PROBABILITY DEXPONENT
      !****************************************************************
    END DO
    !*** MONTE CARLO SWEEPS
    !******************************************************************
    
    !******************************************************************
    !*** CHANGE ANNEALING CONFIGURATION --> cool down

    IF(MOD(i,NAnnealing)  .EQ. 0  .AND. &  ! calculate new annealing parameters after NAnnealing sweeps
       i .NE. NMonteCarlo) THEN            ! --> to log final parameters in case of MOD(NMonteCarl,NAnnealing)=0
      
      !method by sandvik
      !depending on accept/try the magnitude of variation in specFunc is changed
      IF (accept/try.GT.criticaleAcceptRate.AND.variationFact.LT.criticaleVariationInc) THEN
        variationFact=variationFact*variationUpdate !variation from old version
      END IF
      
      IF (accept/try.LE.criticaleAcceptRate.AND.variationFact.GT.criticaleVariationDec) THEN
        variationFact=variationFact/variationUpdate !variation from old version
      END IF
      
      noiseLevel = noiseFilter*MAXVAL(specFunc)  !new configuration for new annealing procedure
                                                 !depending on current spectral density
      specVariation(:)=variationFact*(noiseFilter*specFunc(:)+noiseFilter*noiseLevel)

      annealingBoltzmann=annealingBoltzmann/annealingUpdate ! new annealing temperature
      
      accept=0D0 !reset accept and try counter for every new annelaing config
      try=0D0

      !*** CHANGE ANNEALING CONFIGURATION --> cool down
      !******************************************************************
    
    END IF
    !*** ANNEALING PROCEDURE
    !******************************************************************
    
  END DO
  !*** MONTE CARLO LOOP *****************************************************************
  !**************************************************************************************
  
  !****************************************************************
  !*** LOGGING after the total MonteCarlo procedure
  WRITE(logMsg,"(I0.12,A3,F15.1,I8,I8,F15.5,F20.5)") logIndex,"   ", &
                  chiSquared,INT(accept),INT(try),variationInit,annealingInit

  CALL writeLog (TRIM(logName)//TRIM(logMontecarlo),logMsg,'append',debug_mc)

  logIndex=logIndex+1

  !*** LOGGING
  !******************************************************************
  
END SUBROUTINE montecarloSpectralFunction
!*****************************************************************************************


!*****************************************************************************************
!*** smoothSpectrum **********************************************************************
!*** the spectral function is locally smoothed and normalized depending on the kernelMode
!*** and module-parameter smoothingMethod (and integration methode)
!***
!*** specFunc is the spectral function which is to be smoothed
!*** specGrid is the spectral grid where specFunc is defined
!*** NSpec is the number of spectal gridPoints
!*** smoothCnt is the number of gridPoints left and right from the grid-point (j), where the spectral
!***           function is to be smoothed over (j-smoothCnt  until j+smoothCnt)
!*** smoothingMethod defines the method, how the local average over the (noisy) spectral function is performed
!***           0: locally average spectral function: smoothSpectrum(omega_i)=sum(spectralFunction(k:l))/(l-k+1), k=omega_i-smoothCnt, l=omega_i+smoothCnt, 
!***           1: smoothing by locally integrate spectral function: smoothSpectrum(omega_i)=SUM(domega(k:l)*spectralFunction(k:l))/SUM(domgea(k:l), important for non-equidistant grid-points 
!***           2: smoothing by locally integrate chi''(omega) propTo spectral function, important for non-equidistand grid points for kernel 10 and 11
!*** kernelMode defines the kernel which is used
!*** kernelBeta defines the inverse temperature
!*****************************************************************************************
SUBROUTINE smoothSpectrum (specFunc, specGrid, NSpec,smoothCnt,smoothingMethod,kernelMode,kernelBeta)
IMPLICIT NONE

  INTEGER,                   INTENT(IN)  :: NSpec           ! number of grid points for spectral function
  REAL(KINDR),  DIMENSION(Nspec), INTENT(OUT) :: specFunc        ! spectral function, object of interest
  REAL(KINDR),  DIMENSION(Nspec), INTENT(IN)  :: specGrid        ! grid for spectral function
  
  INTEGER,                   INTENT(IN)  :: smoothingMethod ! defines the method, how the local average over the (noisy) spectral function is performed
  INTEGER,                   INTENT(IN)  :: smoothCnt       ! number of spectral grid points, where the average is taken 
                                                            ! smooth(i)=SUM(specFunc(i-smoothCnt : i+smoothCnt))/(i+smoothCnt-(i-smoothCnt))
  INTEGER,                   INTENT(IN)  :: kernelMode      ! 
  REAL(KINDR)                                 :: kernelBeta

  REAL(KINDR),  DIMENSION(Nspec)              :: specSmooth      ! dummy spectrum for smoothing
  REAL(KINDR)                                 :: specWeight      ! dummy spectrum for smoothing
  REAL(KINDR)                                 :: smoothRange     ! dummy spectrum for smoothing

  INTEGER                                :: i,j,k,l
  
  IF (smoothCnt .GT. 0) THEN

    specSmooth(:)=0D0                             ! default setting

    DO j=1,NSpec                                  ! take average of spectral function over grid points k:l
      k=MAX(j-ABS(smoothCnt),1)                   ! make sure, 1<=k     and k<l : DABS(smoothCnt)>0
      l=MIN(j+ABS(smoothCnt),NSpec)               ! make sure, l<=Nspec and k<l : DABS(smoothCnt)>0

      IF (smoothingMethod .EQ. 0) THEN                ! take average over j-|smoothWidth| : j+|smoothWidth|      
        specSmooth(j) = SUM(specFunc(k:l))/(l-k+1D0)  ! smoothCnt GT 0, l > k        
      ELSE IF (smoothingMethod .EQ. 1 .OR. smoothingMethod .EQ. 2) THEN

        smoothRange=0D0

        DO i=k, l                            ! locally integrate over j-|smoothWidth| : j+|smoothWidth|      
          IF (integrationMethod .EQ. 1) THEN ! trapezoid integration
            specWeight=specGrid(min(i+1,l))-specGrid(max(i-1,k))
          ELSE                               ! default, assume uniform grid, restores smoothingMedhod=0
            specWeight=specGrid(2)-specGrid(1)
          END IF
          
          smoothRange = smoothRange + specWeight

          IF (smoothingMethod .EQ. 1 .OR. kernelMode .EQ. 0 .OR. kernelMode .EQ. 1) THEN !local integral average
            IF (specGrid(i) .EQ. 0D0 .AND. (kernelMode .EQ. 10 .OR. kernelMode .EQ. 11)) THEN 
              specWeight=0D0 ! explicitly avoid spectral weight of pole (at omega=0) of the spectral function in the smoothing procedure
            END IF           ! for kernel-mode 10 and 11
            
            specSmooth(j) = specSmooth(j) + specFunc(i)*specWeight
          ELSE IF (kernelMode .EQ. 10 .AND. smoothingMethod .EQ. 2) THEN                 !local integral average of chi'', kernel 10
            specSmooth(j) = specSmooth(j) + specFunc(i)*specWeight*specGrid(i)           ! TODO: maybe one can interpolate for specGrid(i) --> 0
          ELSE IF (kernelMode .EQ. 11 .AND. smoothingMethod .EQ. 2) THEN                 !local integral average of chi'', kernel 11
            specSmooth(j) = specSmooth(j) + &                                       
                            specFunc(i)*specWeight*(1D0-EXP(-kernelBeta*specGrid(i)))   ! TODO: maybe one can interpolate for specGrid(i) --> 0
          ELSE                                                                           !not implemented, no smoothing
            smoothRange=1D0           ! no smoothing
            specSmooth(j)=specFunc(j) ! don't touch spectral function
          END IF        
        END DO
                
        specSmooth(j)=specSmooth(j)/smoothRange
                
      END IF
      
      IF (kernelMode .EQ. 10 .OR. kernelMode .EQ. 11) THEN  ! apply smoothing of spectral function A(omega)
        IF (smoothingMethod .EQ. 0) THEN
          specFunc(j)=specSmooth(j)
        ELSE IF (smoothingMethod .EQ. 1 .AND. specGrid(j) .NE. 0D0) THEN
          specFunc(j)=specSmooth(j)
        ELSE IF (smoothingMethod .EQ. 2 .AND. specGrid(j) .NE. 0D0) THEN
          IF (kernelMode .EQ. 10) THEN       ! Bosonic kernel, A(omega)=chi''(omega)/omega
            !avoid division by zero
            specFunc(j)=specSmooth(j)/specGrid(j)    
          ELSE IF  (kernelMode .EQ. 11) THEN ! Bosonic spectral function A(omega)=chi''(omega)/(1-exp(-beta*omega))
            !avoid division by zero
            specFunc(j)=specSmooth(j)/(1D0-EXP(-kernelBeta*specGrid(j)))
          END IF          
        END IF
      ELSE 
        specFunc(j)=specSmooth(j)
      END IF

    END DO
  END IF

  CALL normalizeToGrid(specFunc,specGrid,Nspec,kernelMode)  ! normalize spectral function again

END SUBROUTINE
!*****************************************************************************************




!*****************************************************************************************
!*** chiMoments **************************************************************************
!*****************************************************************************************
!*** calculates chiOmegaZero and chiTauZero given by Eq. 5.16 REF[1] and 5.17 REF[1]:
!*** Kernel1: chi(\omega=0) = 2*integral_0^{\infty} d\omega A(\omega)
!***          chi(\tau=0)   = 2*integral_0^{\infty} d\omega A(\omega)*\omega*coth(beta*omega/2)
!*** 
!*** Kernel2: chi(\omega=0) = 2*integral_0^{\infty} d\omega A(\omega)
!***          chi(\tau=0)   = 2*integral_0^{\infty} d\omega A(\omega)*\omega*coth(beta*omega/2)
!*****************************************************************************************
SUBROUTINE chiMoments (specFunct, grid, gridUnit, NGrid, beta, kernelMode, chiOmegaZero, chiTauZero)
  IMPLICIT NONE
  !*** Function Variables ***
  INTEGER,                   INTENT(IN)                 :: NGrid         !number of grid points
  REAL(KINDR),  DIMENSION(NGrid), INTENT(IN)  :: specFunct     !spectral function for specific kernel
  REAL(KINDR),  DIMENSION(NGrid), INTENT(IN)  :: grid          !grid where specFunct is defined, should range from [0,inf]
  REAL(KINDR),  DIMENSION(NGrid), INTENT(IN)  :: gridUnit      !grid integration weights
  REAL(KINDR),                    INTENT(IN)  :: beta          !Boltzmann factor
  INTEGER,                   INTENT(IN)  :: kernelMode!defines normalization mode
  !*** Function Variables ***

  
  !*** Return values ***
  REAL(KINDR),                    INTENT(OUT) :: chiOmegaZero  ! real part of spectral susceptibility (via Kramers Kronig Relation)
  REAL(KINDR),                    INTENT(OUT) :: chiTauZero    ! susceptibility at origin of imaginary time (tau=0 in EQ.5.14 REf[1])
  !*** Return values ***

  !*** Local Helpers ***
  INTEGER                                               :: i
  !*** Local Helpers ***

  chiOmegaZero=0D0
  chiTauZero =0D0

  IF (kernelMode .EQ. 10) THEN ! bosonic kernel
    !*** chi(\omega=0) = 2*integral_0^{\infty} d\omega A(\omega), EQ. 5.16 REF[1]
    chiOmegaZero=2D0*SUM(gridUnit(:)*specFunct(:)) 

    DO i=1, NGrid
      IF (grid(i) .NE. 0D0) THEN
        chiTauZero = chiTauZero + gridUnit(i)*specFunct(i)*grid(i)/TANH(beta*grid(i)/2D0) 
      ELSE 
        ! \omega / tanh(beta*\omega/2) = 2/beta for \omega -> 0
        chiTauZero = chiTauZero + gridUnit(i)*specFunct(i)*2D0/beta 
      END IF        
    END DO

  ELSE IF (kernelMode .EQ. 11) THEN ! bosonic spectral function
    !*** chi(\omega=0) = 2*integral_0^{\infty} d\omega A(\omega), EQ. 5.16 REF[1]
    chiTauZero=SUM(gridUnit(:)*specFunct(:)*(1D0+EXP(-beta*grid(:)))) 

    DO i=1, NGrid
      !*** chi(\omega=0) = 2*integral_0^{\infty} d\omega (1-e^{-\beta*\omega})/ \omega * A(\omega), EQ. 5.16 REF[1]
      IF (grid(i) .NE. 0D0) THEN
        chiOmegaZero = chiOmegaZero + 2D0*gridUnit(i)*specFunct(i)* &
                                          (1D0-EXP(-beta*grid(i)))/grid(i) 
      ELSE 
        ! (1-e^{-\beta*\omega})/ \omega = /beta for \omega -> 0
        chiOmegaZero = chiOmegaZero + 2D0*gridUnit(i)*specFunct(i)*beta 
      END IF        
    END DO
  ELSE
    ! not implemented, errorcode: 
    chiOmegaZero=-1D0
    chiTauZero =-1D0
  END IF 
  
END SUBROUTINE
!*****************************************************************************************




!*****************************************************************************************
!*** getGridUnit *************************************************************************
!*****************************************************************************************
!***
!*** evaluates spacing between grid-points (grid): integration weights or integration units (unit)
!*** integration method is controlled by MODEL-parameter integrationMethod
!*** grid ... grid where integration weights are calculated of
!*** unit ... integration weights (default: equal spacing assumed)
!*** NGrid ... number of grid-points = number of integration weights
!*****************************************************************************************
SUBROUTINE getGridUnit (grid,unit,Ngrid)
  IMPLICIT NONE
  INTEGER,                   INTENT(IN)  :: Ngrid
  REAL(KINDR),  DIMENSION(Ngrid), INTENT(IN)  :: grid
  REAL(KINDR),  DIMENSION(Ngrid), INTENT(OUT) :: unit
  INTEGER                                :: i
  
  IF (integrationMethod .EQ. 1) THEN
    unit(1) = (grid(2)-grid(1))/2D0     ! trapezoid integration weights for (non-)equidistant grid-spacing
    DO i=2,Ngrid-1                      ! grid assumed to contain more than 1 element ...
      unit(i)=(grid(i+1)-grid(i-1))/2D0
    END DO
    unit(Ngrid)=(grid(Ngrid)-grid(Ngrid-1))/2D0
  ELSE
    unit(:)=grid(2)-grid(1)            ! default: assume uniform
  END IF
END SUBROUTINE getGridUnit
!*****************************************************************************************


!*****************************************************************************************
!*** normalizeToGrid *********************************************************************
!*****************************************************************************************
!***
!*** normalize integral of funct over grid (if kernelMode equals 0 or 1):
!*** integral funct(w) dw ~ Sum funct(i)*dGrid(i) = 1, with dGrid(i) being a spectral unit
!***
!*** NGrid     number of grid points
!*** grid      grid where funct is defined
!*** funct     function to be normalized
!*** kernelMode defines if normalization is performed (only for kernelMode 0 and 1)
!*****************************************************************************************
SUBROUTINE normalizeToGrid(funct,grid,NGrid,kernelMode)
  IMPLICIT NONE

  !*** Function Variables ***
  INTEGER,                   INTENT(IN)                   :: NGrid     !number of grid points
  REAL(KINDR),  DIMENSION(NGrid), INTENT(INOUT) :: funct     !function to be normalized
  REAL(KINDR),  DIMENSION(NGrid), INTENT(IN)    :: grid      !grid where funct is defined
  INTEGER,                   INTENT(IN)                   :: kernelMode!defines normalization mode
  REAL(KINDR),  DIMENSION(NGrid)                :: gridUnit  !spacing between grid points
  !*** Function Variables ***

  IF (kernelMode .EQ. 0 .OR. kernelMode .EQ. 1) THEN
    CALL getGridUnit(grid,gridUnit,NGrid)               ! get spacing between grid points = integration weights
    funct(:)=funct(:)/SUM(funct(:)*gridUnit(:))         ! SUM(func(:)*gridUnit(:)) = 1 
  END IF
  
END SUBROUTINE normalizeToGrid 
!*****************************************************************************************
!*****************************************************************************************



!*****************************************************************************************
!*** writeLog ****************************************************************************
!*** 
!*** writes string <logstr> into file with name <filename> with POSITION=<pos>
!*** if <debugLog> is .true., <logstr> also written into UNIT=outputUnit
!*****************************************************************************************
SUBROUTINE writeLog (filename,logstr,pos,debugLog)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)    :: filename
  CHARACTER(*), INTENT(IN)    :: logstr
  CHARACTER(*), INTENT(IN)    :: pos
  LOGICAL, INTENT(IN)         :: debugLog
    
  !Patrik: writes to screen and file only when flag set to true
  IF (debugLog) THEN
    WRITE(outputUnit,"(A1,A)") " ", TRIM(logstr)
    OPEN(UNIT=101,FILE=TRIM(filename),POSITION=pos)
    WRITE(101,"(A1,A)") " ", TRIM(logstr)
    CLOSE(101)
  END IF
END SUBROUTINE
!*****************************************************************************************
!*****************************************************************************************


END MODULE
!*****************************************************************************************

