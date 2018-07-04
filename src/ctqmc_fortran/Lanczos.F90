!===============================================================================
module MLanczos
!===============================================================================
use MStates
use MSparseMatrix

contains

!===============================================================================
subroutine LanczosStep(DStates,iSSt,v1,v2,a,b,iter)
!===============================================================================
!this subroutine calculates the lanczos step:
!a(iter)=<v_iter|H|v_iter>
!b(iter)*v_iter+1=[H-a(iter)]*v_iter-b(iter-1)*v_iter-1
!input
   type(TStates)           :: DStates
   integer                 :: iSSt,iter
! lanczos vector at step i,beta_i * lanczos vector at step i+1
   real(KINDR)             :: v1(1:DStates%SubStates(iSSt)%NStates,PMAX),v2(:)
! diagonal and off diagonal lanczos paramaters
   real(KINDR)             :: a(PMAX),b(PMAX)

   if(iter.eq.1)then
! v1=v_0 (initial vector)
      v2(1:DStates%SubStates(iSSt)%NStates)=0d0
   else
! v1=v_iter-1; v2=b(iter-1)*v_iter
   v1(:,iter)=1d0/b(iter-1)*v2(:)
   v2(:)=-b(iter-1)*v1(:,iter-1)
! v1=v_iter; v2=-b(iter-1)*v_iter-1
   endif
   ! v2=H*v_iter-b(iter-1)*v_iter-1
  call SMatVecProdAdd(DStates%SubStates(iSSt)%Hamiltonian,v1(:,iter),v2)
! a(iter)=<v_iter|H|v_iter>
   a(iter)=dot_product(v1(1:DStates%SubStates(iSSt)%NStates,iter),v2(1:DStates%SubStates(iSSt)%NStates))
! v2=(H-a(iter))*v_iter-b(iter-1)*v_iter-1=b(iter)*v_iter+1
   v2(1:DStates%SubStates(iSSt)%NStates)=v2(1:DStates%SubStates(iSSt)%NStates)-a(iter)*v1(1:DStates%SubStates(iSSt)%NStates,iter)
! b(iter)=norm(v2) 
   b(iter)=sqrt(sum(v2(1:DStates%SubStates(iSSt)%NStates)**2))
   !write(*,*) "aus lanczosstep"
   !write(*,*) "a", a
   !write(*,*) "b", b
   
   
end subroutine LanczosStep

!===============================================================================
subroutine LanczosE0(DStates,iSSt,v1,v2,ee,ev,iter,nconvct,egs,tau)
!===============================================================================
!input
   type(TStates)            :: DStates
   integer                  :: iSSt
   real(KINDR)              :: egs,tau
! lanczos vector at step i,beta_i * lanczos vector at step i+1
   real(KINDR)              :: v1(1:DStates%SubStates(iSSt)%NStates,PMAX),v2(:)
! eigenvalues and eigenvectors of the tridiagonal matrix
   real(KINDR)              :: expev(PMAX)
   real(KINDR)              :: evinv(PMAX,PMAX),ct(PMAX)
   real(KINDR), intent(out) :: ee(PMAX),ev(PMAX,PMAX)
!local
   integer                  :: iter, nconvct, countct
! diagonal and off diagonal lanczos paramaters
   real(KINDR)              :: a(PMAX),b(PMAX)
   
   countct = nconvct
   do iter = 1, PMAX
      call LanczosStep(DStates,iSSt,v1,v2,a,b,iter)
      !write(*,*) "a", a
      !write(*,*) "b", b
      !write(*,*) "iter", iter
      call TriEigen(iter,a,b,ee,ev)
      !write(*,*) "ee", ee
      !write(*,*) "ev", ev
      call MatInvert(ev,evinv,iter)
      expev(1:iter)=exp(-tau*(ee(1:iter)-Egs))*evinv(1:iter,1)
      ct(iter)=sum(ev(iter,1:iter)*expev(1:iter))
      if (b(iter).lt.1d-50) then
         !write(*,*) "exit because of b !!!"
         exit
      endif

      if (countct .lt. nconvct) then
        if ( countct .eq. 0 ) exit
        countct = countct - 1
      endif

      if (ct(iter)**2.lt.EPSLANC.and.countct.eq.nconvct)then
         if ( nconvct .eq. 0 ) then
            !write(*,*) "ct(iter)**2", ct(iter)**2
            !write(*,*) "abbruch im dritten !!!"
            exit 
         endif
         countct = countct - 1
      endif
   enddo

end subroutine LanczosE0

!===============================================================================
subroutine LanczosV0(DStates,iSSt,v1,v2,Psi0,ee,ev,NIter)
!===============================================================================
!input
   type(TStates)           :: DStates
   integer                 :: iSSt,NIter
! lanczos vector at step i,lanczos vector at step i+1
   real(KINDR)             :: v1(1:DStates%SubStates(iSSt)%NStates,PMAX),v2(:)
! eigenvalues and eigenvectors of the tridiagonal matrix
   real(KINDR)             :: ee(PMAX),ev(PMAX,PMAX)
!output
! after this routine the ground state eigenvector
   real(KINDR)             :: Psi0(:)
!local
! diagonal and off diagonal lanczos paramaters
   real(KINDR)             :: a(PMAX),b(PMAX),resid
   integer                 :: iter

! calculate ground state vector Psi1
   Psi0=0d0
   do iter=1,min(NIter,PMAX)
      call LanczosStep(DStates,iSSt,v1,v2,a,b,iter)
      Psi0=Psi0+v1(1:DStates%SubStates(iSSt)%NStates,iter)*ev(iter,1)
   enddo
! check the accuracy of the ground state vector   
   Psi0=Psi0/sqrt(sum(Psi0**2))
! v1=H*Psi-E*Psi
   v1(:,1)=MatVecProd(DStates%SubStates(iSSt)%Hamiltonian,Psi0)-ee(1)*Psi0
   resid=sqrt(sum(v1(1:DStates%SubStates(iSSt)%NStates,1)**2))
   write(*,*)"resid= ",resid
end subroutine LanczosV0

!===============================================================================
subroutine LanczosTransform(DStates,iSSt,v1,ct,Psi,NIter)
!===============================================================================
! Transforms a vector in Krylov space back into the original space. See
! T.J.Park and al. J.Chem.Phys. 85(10) 5870
!input
   type(TStates)              :: DStates
   integer                    :: iSSt,NIter
! lanczos vector at step i
   real(KINDR)                :: v1(1:DStates%SubStates(iSSt)%NStates,PMAX)
! evolved Lanczos state
   real(KINDR)                :: ct(PMAX)
!output
! backtransformed evolved state
   real(KINDR)                :: Psi(:)
!local
! lanczos vector at step i+1
! diagonal and off diagonal lanczos paramaters
   integer                    :: iter
   
   Psi=0d0
   do iter=1,min(NIter,PMAX)
      Psi=Psi+v1(1:DStates%SubStates(iSSt)%NStates,iter)*ct(iter)
   enddo
end subroutine LanczosTransform

!===============================================================================
subroutine LanczosTimeEvolve(DStates,iSSt,Psi_t,Psi_tp,tau,Norm,Egs,nconvct,lsteps)
!===============================================================================
!input
   type(TStates)              :: DStates
   integer                    :: iSSt
   real(KINDR)                :: tau,Egs
!input Psi_t before time evolution
   real(KINDR)                :: Psi_t(:)
!output Psi_tp after time evolution
   real(KINDR)                :: Psi_tp(:)
!output
   real(KINDR)                :: Norm
! in/out if input 0 then output is the number of steps Lanczos needed for
! convergence. if the input is gt 0 steps is the number of steps Lanczos
! should use.
   integer                    :: lsteps, nconvct
!local
   real(KINDR)                :: ee(PMAX),ev(PMAX,PMAX),evinv(PMAX,PMAX),ct(PMAX)
   real(KINDR)                :: expev(PMAX),one(1)
   real(KINDR)                :: LanczosV(1:DStates%SubStates(iSSt)%NStates,PMAX)
   real(KINDR)                :: LanczosV3(1:DStates%SubStates(iSSt)%NStates)
   integer                    :: i

   Norm=sqrt(sum(Psi_t**2))
   if(Norm==0d0)return
   if(size(Psi_t).eq.1)then
      !write(*,*) "lanczos dimension one"
      one(1)=1d0
      call SMatVecProd(DStates%SubStates(iSSt)%Hamiltonian,one,ee)
      !write(*,*) "ee(1)", ee(1)
      !write(*,*) "egs", egs
      !write(*,*) "norm", norm
      !write(*,*) "tau", tau
      Psi_tp=Psi_t*exp(-tau*(ee(1)-Egs))/Norm
      lsteps = 1
      return
   endif 
   !write(*,*) "norm", norm
   LanczosV(1:DStates%SubStates(iSSt)%NStates,1)=Psi_t/Norm
   call LanczosE0(DStates,iSSt,LanczosV,LanczosV3,ee,ev,lsteps,nconvct,egs,tau)
   lsteps=min(PMAX,lsteps)

   call MatInvert(ev,evinv,LSteps)
   ! Time evolution in Krylov space
   !write(*,*) "tau", tau
   !write(*,*) "ee(1:lsteps)", ee(1:lsteps)
   !write(*,*) "egs", egs
   !write(*,*) "evinv(1:LSteps,1)", evinv(1:LSteps,1)
   expev(1:LSteps)=exp(-tau*(ee(1:LSteps)-Egs))*evinv(1:LSteps,1)
   do i=1,LSteps
      ct(i)=sum(ev(i,1:LSteps)*expev(1:LSteps))
   enddo
   call LanczosTransform(DStates,iSSt,LanczosV,ct,Psi_tp,lsteps)
end subroutine LanczosTimeEvolve


!===============================================================================
subroutine TriEigen(n,a,b,ee,ev)
!===============================================================================
! Calculates the eigenvalues (ee) and eigenvectors (ev) for a tridiagonal
! symmetric matrix (bands saved in a,b) with the dimension n.
!input
   integer, intent(in)      :: n
   real(KINDR), intent(in)  :: a(PMAX),b(PMAX)
   real(KINDR), intent(out) :: ee(PMAX),ev(PMAX,PMAX)
!local
   real(KINDR)              :: sd(PMAX),work(2*PMAX-2),d
   integer                  :: info
   
   if(n.eq.1)then
      ee(1)=a(1)
      ev(1,1)=1d0
   elseif(n.eq.2)then
      d=sqrt((a(1)-a(2))**2/4d0+b(1)**2)
      ee(1)=(a(1)+a(2))/2d0-d
      ee(2)=(a(1)+a(2))/2d0+d
      ev(1,1)=-1d0/sqrt(1d0+((ee(1)-a(1))/b(1))**2)
      ev(1,2)=-ev(1,1)*(ee(1)-a(1))/b(1)
      ev(2,1)=-ev(1,2)
      ev(2,2)=ev(1,1)
   else
!      if(n.eq.3)then
!      c(1)=-(a(1)+a(2)+a(3))
!      c(2)=-(b(1)**2+b(2)**2-a(1)*a(2)-a(1)*a(3)-a(2)*a(3))
!      c(3)=-(a(1)*a(2)*a(3)-a(1)*b(2)**2-a(3)*b(1)**2)
!      p=(3d0*c(2)-c(1)**2)/9d0
!      q=c(1)**3/27d0-c(1)*c(2)/6d0+c(3)/2d0
!      r=sign(sqrt(abs(p)),p)
!      phi=acos(q/r**3)
!      y(1)=-2d0*r*cos(phi/3d0)
!      y(2)=2d0*r*cos((pi-phi)/3d0)
!      y(3)=2d0*r*cos((pi+phi)/3d0)
!      ee(1)=y(2)-c(1)/3d0
!      ee(2)=y(3)-c(1)/3d0
!      ee(3)=y(1)-c(1)/3d0
!      ev(1,1)=sqrt(1d0/(1+((a(1)-ee(1))/b(1))**2+((a(1)-ee(1))/b(1))**2*(b(2)/(a(3)-ee(1)))**2))
!      ev(1,2)=sqrt(1d0/((b(1)/(a(1)-ee(1)))**2+1+(b(2)/(a(3)-ee(1)))**2))
!      ev(1,3)=-sqrt(1d0/(1+((a(3)-ee(1))/b(2))**2+((a(3)-ee(1))/b(2))**2*(b(1)/(a(1)-ee(1)))**2))
!!      ev(1,1)=b(1)*ev(1,2)/(a(1)-ee(1))
!!      ev(1,3)=b(2)*ev(1,2)/(a(3)-ee(1))
!!      ev(2,1)=b(1)*ev(1,2)/(a(1)-ee(2))
!!      ev(2,3)=b(2)*ev(1,2)/(a(3)-ee(2))
!      ev(3,1)=sqrt(1d0/(1+((a(1)-ee(3))/b(1))**2+((a(1)-ee(3))/b(1))**2*(b(2)/(a(3)-ee(3)))**2))
!      ev(3,2)=-sqrt(1d0/((b(1)/(a(1)-ee(3)))**2+1+(b(2)/(a(3)-ee(3)))**2))
!      ev(3,3)=-sqrt(1d0/(1+((a(3)-ee(3))/b(2))**2+((a(3)-ee(3))/b(2))**2*(b(1)/(a(1)-ee(3)))**2))
!!      ev(3,1)=b(1)*ev(1,2)/(a(1)-ee(3))
!!      ev(3,3)=b(2)*ev(1,2)/(a(3)-ee(3))
!      ev(2,1:3)=1d0
!      ev(2,1:3)=ev(2,1:3)-dot_product(ev(2,1:3),ev(1,1:3))*ev(1,1:3)-dot_product(ev(2,1:3),ev(3,1:3))*ev(3,1:3)
!      ev(2,1:3)=ev(2,1:3)/sqrt(dot_product(ev(2,1:3),ev(2,1:3)))
!      ev(2,1)=-ev(2,1)
!      ev(2,3)=-ev(2,3)
!!      ev(2,1)=sqrt(1d0/(1+((a(1)-ee(2))/b(1))**2+((a(1)-ee(2))/b(1))**2*(b(2)/(a(3)-ee(2)))**2))
!!      ev(2,2)=sqrt(1d0/((b(1)/(a(1)-ee(2)))**2+1+(b(2)/(a(3)-ee(2)))**2))
!!      ev(2,3)=sqrt(1d0/(1+((a(3)-ee(2))/b(2))**2+((a(3)-ee(2))/b(2))**2*(b(1)/(a(1)-ee(2)))**2))
!      write(*,*)"ana"
!      write(*,*)ev(1:3,1:3) 
!      endif
!   else
      ee(1:n)=a(1:n)
      sd(1:n)=b(1:n)
#ifdef LAPACK77_Interface      
      call dstev('V',n,ee,sd,ev,PMAX,work,info)
#endif
#ifdef LAPACK95_Interface
      call stev(ee,sd,ev,info)
#endif
!      if(n.eq.3)then 
!!         write(*,*)ee(1:3)
!         write(*,*)"num"
!         write(*,*)ev(1:3,1:3)
!      endif
   endif
end subroutine TriEigen

!===============================================================================
subroutine MatInvert(A,B,d)
!===============================================================================
! calculates the inverse of a matrix A with dimension d and stores it in A
! we only need the first column of the inverse matrix. this is used
! in the analytic calculation of the inverse matrix.
!input
   integer              :: d
   real(KINDR)          :: A(PMAX,PMAX)
   real(KINDR)          :: B(PMAX,PMAX)
!local
   integer              :: info,ipiv(d)
   real(KINDR)          :: work(d*64),det

   if(d.eq.1)then
      B(1,1)=1d0/A(1,1)
   elseif(d.eq.2)then
      det=(A(1,1)*A(2,2)-A(1,2)*A(2,1))
      B(1,1)=A(2,2)/det
      B(2,1)=-A(2,1)/det
   elseif(d.eq.3)then
      B(1,1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      B(2,1)=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      B(3,1)=A(2,1)*A(3,2)-A(2,2)*A(3,1)
      det=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      B(1:d,1)=B(1:d,1)/det
   else
      B(1:d,1:d)=A(1:d,1:d)
#ifdef LAPACK77_Interface
      call DGETRF(d,d,B,PMAX,ipiv,info)
      call DGETRI(d,B,PMAX,ipiv,work,d*64,info)
#endif
#ifdef LAPACK95_Interface
      call GETRF(B,ipiv,info)
      call GETRI(B,ipiv,info)
#endif
      if(info.ne.0)then
         write(stdout,*)"Error in MLanczos:MatInvert"
         stop
      endif
   endif
end subroutine MatInvert

!===============================================================================
end module MLanczos
!===============================================================================

#ifdef Lanczos_Test

!===============================================================================
program Prog_Lanczos
!===============================================================================
use MParameters
use MStates 
use MOperator
use MRandom
use MLanczos

type(TStates)                    :: DStates
real(KINDR),allocatable                :: HEValues(:), muimp(:,:)
type(TPsis)                            :: DPsis
type(TOperator)                        :: DH
type(TOperator)                        :: HEVectors

real(KINDR),pointer              :: Coeff(:),Psi1(:),Psi2(:),LanczosV(:),evinvsys(:,:)
real(KINDR)                      :: ee(PMAX),ev(PMAX,PMAX),evinv(PMAX,PMAX),ct(PMAX),c(PMAX)
real(KINDR)                      :: Norm
integer                          :: iSt,NIter,i,j
integer,parameter                :: SSt=4 !Substate used
real(KINDR)                      :: tau=2d0

call read_ParameterFile("Parameters.in")


call init_States(DStates)
call init_SubStates(DStates)
call init_Psi(DPsis,DStates)
allocate(muimp(DStates%NBands,2))
muimp(:,:) = get_Real_Parameter("mu") 
call init_Hamiltonian(DH,DStates,muimp)
call init_TOperator(HEVectors,DStates)
allocate(HEValues(0:DStates%NStates-1))
call diag_Operator(DH,DStates,HEVectors,HEValues)
call set_Psis(DPsis,DStates)
call set_Hamiltonian(DH,DStates)
call set_HEVal(HEValues,DStates)
call set_HEVec(HEVectors,DStates)
deallocate(HEValues)
call dest_Psis(DPsis)
call dest_TOperator(DH)
call dest_TOperator(HEVectors)


write(*,*)"Using the substate: ",SSt

allocate(LanczosV(0:DStates%SubStates(SSt)%NStates-1))
allocate(evinvsys(0:DStates%SubStates(SSt)%NStates-1,0:DStates%SubStates(SSt)%NStates-1))
allocate(Coeff(0:DStates%SubStates(SSt)%NStates-1))
allocate(Psi1(0:DStates%SubStates(SSt)%NStates-1))
allocate(Psi2(0:DStates%SubStates(SSt)%NStates-1))

do i=1,2
if(i.eq.1)then
! Small test case also written in Mathematica for comparison with the result one
! obtains with the Mathematica function Matexp[].
write(stdout,*)"================================================================"
write(stdout,*)"              Mathematica test case"
write(stdout,*)"================================================================"
if(get_Integer_Parameter("Nd").ne.2.and.SSt.ne.4)then
   write(stdout,*)"Not the right parameters for test case! Moving on."
   cycle
endif
Coeff(0)=0.50119979671648540d0
Coeff(1)=0.14763127393389930d0
Coeff(2)=0.79040199766821562d0
Coeff(3)=0.31979439146129807d0

elseif(i.eq.2)then
! Start with a random state.
write(stdout,*)"================================================================"
write(stdout,*)"           creating a random mixed state"
write(stdout,*)"================================================================"
call sgrnd(int(time()))
do iSt=0,size(Coeff)-1
   Coeff(iSt)=grnd()
enddo
endif

! We need the norm to a) normalize our vector before Lanczos and to b) multiply
! the resulting vector from Lanczos with the norm again. This is possible since
! every opertation commutes with this number.
Norm=sqrt(sum(Coeff**2))
write(stdout,*)"Mixed state: "
write(stdout,*)Coeff
write(stdout,*)"Norm: "
write(stdout,*)Norm

! Below the ground state energy and the ground state eigenvector is calculated with 
! Lanczos and the eigenvector.
write(stdout,*)"================================================================"
write(stdout,*)"           calculating the ground state of the mixed state"
write(stdout,*)"================================================================"

LanczosV=Coeff/Norm
call LanczosE0(DStates,SSt,LanczosV,ee,ev,NIter)
write(stdout,*)"Number of iterations ",NIter
write(stdout,*)"Lowest energy: ",minval(DStates%SubStates(SSt)%EVal)
write(stdout,*)"Lowest energy Lanczos: ", ee(1)

evinv=ev
call MatInvert(evinv(1:min(PMAX,NIter),1:min(PMAX,NIter)),min(PMAX,NIter))

LanczosV=Coeff/Norm
call LanczosV0(DStates,SSt,LanczosV,Psi1,ee,ev,NIter)
write(stdout,*)"Eigenvector to the lowest energy: "
write(stdout,*)DStates%SubStates(SSt)%EVec(:,minloc(DStates%SubStates(SSt)%EVal))
write(stdout,*)"Eigenvector to the lowest energy Lanczos:"
write(stdout,*)Psi1

! Below the exact evolution (obtained with eigenvectors and eigenenergies from
! the exact diagonalization stored in DStates) and the approximate Lanczos time
! evolution (see T.J. Park et al. J. Chem. Phys. 85(10) 5870) is calculated. 
write(stdout,*)"================================================================"
write(stdout,*)"           evolving the mixed state by tau = ",tau
write(stdout,*)"================================================================"

! Exact time evolution
evinvsys=DStates%SubStates(SSt)%EVec
call MatInvert(evinvsys,DStates%SubStates(SSt)%NStates)
forall(iSt=0:size(evinvsys(1,:))-1)evinvsys(iSt,:)=evinvsys(iSt,:)*exp(DStates%SubStates(SSt)%EVal(iSt)*(-tau))
Psi1=Matmul(Matmul(DStates%SubStates(SSt)%EVec,evinvsys),Coeff)
write(stdout,*)"Evolved Coeff: "
write(stdout,*)Psi1
evinvsys=DStates%SubStates(SSt)%EVec
call MatInvert(evinvsys,DStates%SubStates(SSt)%NStates)
forall(iSt=0:size(evinvsys(1,:))-1)evinvsys(iSt,:)=evinvsys(iSt,:)*exp(DStates%SubStates(SSt)%EVal(iSt)*(tau))
Psi2=Matmul(Coeff,Matmul(DStates%SubStates(SSt)%EVec,evinvsys))
write(stdout,*)"Evolved inverse coeff: "
write(stdout,*)Psi2
write(stdout,*)"Norm in normal space: ",dot_product(Psi1,Psi2)

! Time evolution in Krylov space
do j=1,min(NIter,PMAX)
   ct(j)=sum(ev(j,1:min(NIter,PMAX))*exp(-tau*ee(1:min(PMAX,NIter)))*evinv(1:min(PMAX,NIter),1))
   c(j)=sum(ev(j,1:min(NIter,PMAX))*exp(tau*ee(1:min(PMAX,NIter)))*evinv(1:min(PMAX,NIter),1))
enddo
write(stdout,*)"tau,norm,|ct_p|^2= ",tau,dot_product(ct(1:min(NIter,PMAX)),c(1:min(NIter,PMAX))),&
   ct(min(NIter,PMAX))*c(min(NIter,PMAX))

! Transformation from Krylov space
LanczosV=Coeff/Norm
call LanczosTransform(DStates,SSt,LanczosV,ct,Psi2,NIter)
write(stdout,*)"Evolved coeff Lanczos: "
write(stdout,*)Psi2*Norm
write(stdout,*)"Squared difference: ",tau, sum((Psi1-Psi2*Norm)**2)
LanczosV=Coeff/Norm
call LanczosTransform(DStates,SSt,LanczosV,c,Psi1,NIter)
write(stdout,*)"Evolved inverse coeff Lanczos: "
write(stdout,*)Psi1*Norm
write(stdout,*)"Norm in normal space Lanczos: ",dot_product(Psi1*Norm,Psi2*Norm)

! Below the eigenvector for the lowest energy is calculated for the evolved state.
Norm=sqrt(sum(Psi2**2))
LanczosV=Psi2/Norm
call LanczosE0(DStates,SSt,LanczosV,ee,ev,NIter)
LanczosV=Psi2/Norm
call LanczosV0(DStates,SSt,LanczosV,Psi2,ee,ev,NIter)
write(stdout,*)"Eigenvector to the lowest energy Lanczos for the evolved state:"
write(stdout,*)Psi2
enddo

deallocate(LanczosV)
deallocate(evinvsys)
deallocate(Coeff)
deallocate(Psi1)
deallocate(Psi2)

end program Prog_Lanczos

#endif
