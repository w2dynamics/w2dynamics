!===============================================================================
module MMatrixUpdate
!===============================================================================
use MParameters
   type TLogDet
      real(KINDR) :: log, sign
   end type TLogDet

   interface operator(*)
      module procedure multiplyLogDet
   end interface operator(*)

   interface operator(/)
      module procedure divideLogDet
   end interface operator(/)

contains

pure elemental function multiplyLogDet(det1, det2) result(det3)
   type(TLogDet), intent(in)  :: det1, det2
   type(TLogDet)              :: det3
   det3%log = det1%log + det2%log
   det3%sign = det1%sign * det2%sign
end function multiplyLogDet

pure elemental function divideLogDet(det1, det2) result(det3)
   type(TLogDet), intent(in)  :: det1, det2
   type(TLogDet)              :: det3
   det3%log = det1%log - det2%log
   det3%sign = det1%sign * det2%sign
end function divideLogDet

pure elemental real(KINDR) function detval(det)
   type(TLogDet), intent(in) :: det
   detval = det%sign * exp(det%log)
end function detval

!===============================================================================
! Calculate determinant det of N(k, k) and inverse (stored in N on exit).
! LU decomposition version
subroutine get_MatLogDetFull(k,N,det)
!===============================================================================
!input
   integer, intent(in)         :: k
   real(KINDR), intent(inout)  :: N(k,k)
!output
   type(TLogDet), intent(out)  :: det
!local
   integer                     :: ipiv(k),info,i
   real(KINDR)                 :: work(k*64)

   det%log = 0.D0
   det%sign = 1.D0
   if(k.eq.0)then
      return
   endif
   !write(*,*) " "
   !write(*,*) "_________________"
   !write(*,*) "in get_MatDetFull"
   !call Ausgabe(N,shape(N),"N",5,.false.)

#ifdef LAPACK77_Interface
   call DGETRF(k,k,N,k,ipiv,info)
#endif
#ifdef LAPACK95_Interface
   call GETRF(N,ipiv,info)
#endif

   if (info.lt.0) then
      write(stderr,*)"Error in MMatrixUpdate::get_MatLogDetFull DGETRF", info
      stop
   endif

   if (info.gt.0) then
      det%log = -huge(0.0_KINDR)
      det%sign = 0.0_KINDR
      return
   endif

   do i=1,k
      if (ipiv(i).ne.i) then
         det%sign = -det%sign * sign(1.0_KINDR, N(i, i))
      else
         det%sign = det%sign * sign(1.0_KINDR, N(i, i))
      end if
      det%log = det%log + log(abs(N(i, i)))
   enddo

#ifdef LAPACK77_Interface
   call DGETRI(k,N,k,ipiv,work,k*64,info)
#endif
#ifdef LAPACK95_Interface
   call GETRI(N,ipiv,info)
#endif
   if(info.ne.0)then
      write(stderr,*)"Error in MMatrixUpdate::get_MatLogDetFull DGETRI"
      stop
   endif
end subroutine get_MatLogDetFull

! !===============================================================================
! ! Calculate determinant det of N(k, k) and inverse (stored in N on exit).
! ! Singular value decomposition version (using QR iteration algorithm
! ! implemented in LAPACK subroutine DBDSQR)
! subroutine get_MatLogDetFull_DBDSQR(k,N,det)
! !===============================================================================
! !input
!    integer, intent(in)         :: k
!    real(KINDR),intent(inout)   :: N(k, k)
! !output
!    type(TLogDet), intent(out)  :: det
! !local
!    integer                     :: info, i, j, lwork
!    real(KINDR)                 :: D(k), E(k), tauQ(k), tauP(k), U(k, k), Vt(k, k)
!    real(KINDR), allocatable    :: work(:)

!    det = TLogDet(log=0.0_KINDR, sign=1.0_KINDR)
!    if(k.eq.0)then
!       return
!    endif

!    ! transform input matrix N to bidiagonal form
!    allocate(work(1))
!    call DGEBRD(k, k, N, k, D, E, tauQ, tauP, work, -1, info)
!    lwork = int(work(1), kind(lwork))
!    deallocate(work)
!    allocate(work(lwork))
!    call DGEBRD(k, k, N, k, D, E, tauQ, tauP, work, lwork, info)
!    deallocate(work)

!    ! multiply factors of the determinant of input matrix N originating
!    ! from elementary Householder matrices (-1 for each one present
!    ! one, i.e. with non-zero coefficient tau)
!    det = TLogDet(log = 0.0_KINDR, sign = 1.0_KINDR)
!    do i = 1, k - 1
!       if (tauQ(i) /= 0.0_KINDR) det%sign = -det%sign
!       if (tauP(i) /= 0.0_KINDR) det%sign = -det%sign
!    end do
!    if (tauQ(k) /= 0.0_KINDR) det%sign = -det%sign

!    ! multiply factors of the determinant of input matrix N in the diagonal
!    ! of the bidiagonal matrix
!    do i = 1, k
!       det%sign = det%sign * sign(1.0_KINDR, D(i))
!       det%log = det%log + log(abs(D(i)))
!    end do

!    do i = 1, k
!       do j = 1, k
!          U(i, j) = 0.0_KINDR
!          Vt(i, j) = 0.0_KINDR
!       end do
!       U(i, i) = 1.0_KINDR
!       Vt(i, i) = 1.0_KINDR
!    end do

!    ! generate matrix U for bidiagonal transformation such that N = U * B * Vt,
!    ! where N is the original input matrix, B the bidiagonal form, and U and Vt
!    ! orthogonal
!    U = N
!    allocate(work(1))
!    call DORGBR('Q', k, k, k, U, k, tauQ, work, -1, info)
!    lwork = int(work(1), kind(lwork))
!    deallocate(work)
!    allocate(work(lwork))
!    call DORGBR('Q', k, k, k, U, k, tauQ, work, lwork, info)
!    deallocate(work)

!    ! generate matrix Vt, see above
!    Vt = N
!    allocate(work(1))
!    call DORGBR('P', k, k, k, Vt, k, tauP, work, -1, info)
!    lwork = int(work(1), kind(lwork))
!    deallocate(work)
!    allocate(work(lwork))
!    call DORGBR('P', k, k, k, Vt, k, tauP, work, lwork, info)
!    deallocate(work)

!    ! generate singular value decomposition, N = U * S * Vt, where N is the
!    ! original input matrix, U and Vt orthogonal, and S diagonal with only
!    ! non-negative entries (stored in D)
!    allocate(work(4 * k))
!    call DBDSQR('U', k, k, k, 0, D, E, Vt, k, U, k, N, k, work, info)
!    deallocate(work)
!    if (info > 0) then
!       write (stderr, *) "Error in singular value decomposition in get_MatLogDetFull, DBDSQR returned info ", info
!       stop
!    end if

!    ! calculate pseudoinverse N^(-1) = V * S^(-1) * Ut
!    ! = transpose(Vt) * transpose(U * S^(-1)) using S = transpose(S)
!    ! divide columns of U by singular values to get U * S^(-1)
!    do i = 1, k
!       if (D(i) /= 0.0_KINDR) then
!          U(:, i) = U(:, i) / D(i)
!       else
!          U(:, i) = 0.0_KINDR
!       end if
!    end do

!    ! multiply transposed matrices
!    call DGEMM('C', 'C', k, k, k, 1.0_KINDR, Vt, k, U, k, 0.0_KINDR, N, k)
! end subroutine get_MatLogDetFull_DBDSQR

! !===============================================================================
! ! Calculate determinant det of N(k, k) and inverse (stored in N on exit).
! ! Singular value decomposition version (using divide-and-conquer
! ! algorithm implemented in LAPACK subroutine DBDSDC)
! subroutine get_MatLogDetFull_DBDSDC(k,N,det)
! !===============================================================================
! !input
!    integer, intent(in)         :: k
!    real(KINDR),intent(inout)   :: N(k, k)
! !output
!    type(TLogDet), intent(out)  :: det
! !local
!    integer                     :: info, i, j, lwork, iwork(8 * k)
!    real(KINDR)                 :: D(k), E(k), tauQ(k), tauP(k), U(k, k), Vt(k, k)
!    real(KINDR), allocatable    :: work(:)

!    det = TLogDet(log=0.0_KINDR, sign=1.0_KINDR)
!    if(k.eq.0)then
!       return
!    endif

!    ! transform input matrix N to bidiagonal form
!    allocate(work(1))
!    call DGEBRD(k, k, N, k, D, E, tauQ, tauP, work, -1, info)
!    lwork = int(work(1), kind(lwork))
!    deallocate(work)
!    allocate(work(lwork))
!    call DGEBRD(k, k, N, k, D, E, tauQ, tauP, work, lwork, info)
!    deallocate(work)

!    ! multiply factors of the determinant of input matrix N originating
!    ! from elementary Householder matrices (-1 for each one present
!    ! one, i.e. with non-zero coefficient tau)
!    det = TLogDet(log = 0.0_KINDR, sign = 1.0_KINDR)
!    do i = 1, k - 1
!       if (tauQ(i) /= 0.0_KINDR) det%sign = -det%sign
!       if (tauP(i) /= 0.0_KINDR) det%sign = -det%sign
!    end do
!    if (tauQ(k) /= 0.0_KINDR) det%sign = -det%sign

!    ! multiply factors of the determinant of input matrix N in the diagonal
!    ! of the bidiagonal matrix
!    do i = 1, k
!       det%sign = det%sign * sign(1.0_KINDR, D(i))
!       det%log = det%log + log(abs(D(i)))
!    end do

!    ! generate singular value decomposition, B = U * S * Vt, where B is the
!    ! bidiagonalized matrix, U and Vt orthogonal, and S diagonal with only
!    ! non-negative entries (stored in D)
!    allocate(work(3 * k**2 + 4 * k))
!    call DBDSDC('U', 'I', k, D, E, U, k, Vt, k, work, iwork, work, iwork, info)
!    deallocate(work)
!    if (info > 0) then
!       write (stderr, *) "Error in singular value decomposition in get_MatLogDetFull, DBDSDC returned info ", info
!       stop
!    end if

!    ! multiply orthogonal matrices from bidiagonalization and SVD to
!    ! get total orthogonal transformation matrices U and Vt such that
!    ! N = U * S * Vt, where N is the original input matrix and S diagonal
!    ! with only non-negative entries (stored in D)
!    allocate(work(1))
!    call DORMBR('Q', 'L', 'N', k, k, k, N, k, tauQ, U, k, work, -1, info)
!    lwork = int(work(1), kind(lwork))
!    deallocate(work)
!    allocate(work(lwork))
!    call DORMBR('Q', 'L', 'N', k, k, k, N, k, tauQ, U, k, work, lwork, info)
!    deallocate(work)

!    ! generate matrix Vt, see above
!    allocate(work(1))
!    call DORMBR('P', 'R', 'T', k, k, k, N, k, tauP, Vt, k, work, -1, info)
!    lwork = int(work(1), kind(lwork))
!    deallocate(work)
!    allocate(work(lwork))
!    call DORMBR('P', 'R', 'T', k, k, k, N, k, tauP, Vt, k, work, lwork, info)
!    deallocate(work)

!    ! calculate pseudoinverse N^(-1) = V * S^(-1) * Ut
!    ! = transpose(Vt) * transpose(U * S^(-1)) using S = transpose(S)
!    ! divide columns of U by singular values to get U * S^(-1)
!    do i = 1, k
!       if (D(i) /= 0.0_KINDR) then
!          U(:, i) = U(:, i) / D(i)
!       else
!          U(:, i) = 0.0_KINDR
!       end if
!    end do

!    ! multiply transposed matrices
!    call DGEMM('C', 'C', k, k, k, 1.0_KINDR, Vt, k, U, k, 0.0_KINDR, N, k)
! end subroutine get_MatLogDetFull_DBDSDC

!===============================================================================
! Calculate determinant det of N(k, k).
! LU decomposition version
subroutine get_LogDetFull(k,N,det)
!===============================================================================
!input
   integer, intent(in)        :: k
   real(KINDR), intent(in)    :: N(k,k)
   type(TLogDet), intent(out) :: det
!local
   integer                    :: ipiv(k),info,i
   real(KINDR)                :: TN(k,k)
   
   if(k.eq.1)then
      det%sign = sign(1.0_KINDR, N(1, 1))
      det%log = log(abs(N(1,1)))
      return
   endif
   TN=N
#ifdef LAPACK77_Interface
   call DGETRF(k,k,TN,k,ipiv,info)
#endif
#ifdef LAPACK95_Interface
   call GETRF(TN,ipiv,info)
#endif
   if (info.ne.0) then
      write(stderr,*)"Error in MMatrixUpdate::get_LogDetFull DGETRF"
      stop
   endif

   det%log = 0.D0
   det%sign = 1.D0
   do i=1,k
      if (ipiv(i).ne.i) then
         det%sign = -det%sign * sign(1.0_KINDR, TN(i, i))
      else
         det%sign = det%sign * sign(1.0_KINDR, TN(i, i))
      endif
      det%log = det%log + log(abs(TN(i, i)))
   enddo
end subroutine get_LogDetFull

! !===============================================================================
! ! Calculate determinant det of N(k, k).
! ! Bidiagonal decomposition version
! subroutine get_LogDetFull_BD(k,N,det)
! !===============================================================================
! !input
!    integer, intent(in)        :: k
!    real(KINDR), intent(in)    :: N(k, k)
!    type(TLogDet), intent(out) :: det
! !local
!    integer                    :: info, i, lwork
!    real(KINDR)                :: D(k), E(k), tauQ(k), tauP(k), TN(k, k)
!    real(KINDR), allocatable   :: work(:)

!    if (k == 1) then
!       det%sign = sign(1.0_KINDR, N(1, 1))
!       det%log = log(abs(N(1, 1)))
!       return
!    endif

!    TN = N
!    ! transform input matrix N to bidiagonal form
!    allocate(work(1))
!    call DGEBRD(k, k, TN, k, D, E, tauQ, tauP, work, -1, info)
!    lwork = int(work(1), kind(lwork))
!    deallocate(work)
!    allocate(work(lwork))
!    call DGEBRD(k, k, TN, k, D, E, tauQ, tauP, work, lwork, info)
!    deallocate(work)

!    ! multiply factors of the determinant of input matrix N originating
!    ! from elementary Householder matrices (-1 for each one present
!    ! one, i.e. with non-zero coefficient tau)
!    det = TLogDet(log = 0.0_KINDR, sign = 1.0_KINDR)
!    do i = 1, k - 1
!       if (tauQ(i) /= 0.0_KINDR) det%sign = -det%sign
!       if (tauP(i) /= 0.0_KINDR) det%sign = -det%sign
!    end do
!    if (tauQ(k) /= 0.0_KINDR) det%sign = -det%sign

!    ! multiply factors of the determinant of input matrix N in the diagonal
!    ! of the bidiagonal matrix
!    do i = 1, k
!       det%sign = det%sign * sign(1.0_KINDR, D(i))
!       det%log = det%log + log(abs(D(i)))
!    end do
! end subroutine get_LogDetFull_BD

! The following subroutines use inversion by partitioning see numerical recipes
! in C p. 77f second edition.
!===============================================================================
real(KINDR) function get_DetRatAdd(N,k,Q,R,S)
!===============================================================================
!input
   integer              :: k
   real(KINDR)          :: Q(k),R(k),S
   real(KINDR)          :: N(k,k)
! local
   integer              :: i,j

   get_DetRatAdd=S
   do i=1,k
      do j=1,k
         get_DetRatAdd=get_DetRatAdd-R(i)*N(i,j)*Q(j)
      enddo
   enddo
end function get_DetRatAdd

!===============================================================================
function get_MatAdd(N,k,Q,R,DetRat,pos)
!===============================================================================
   real(KINDR),pointer  :: get_MatAdd(:,:)
!input
   integer, intent(in)  :: k
   real(KINDR)          :: Q(k),R(k),DetRat
   real(KINDR)          :: N(k,k)
   integer, intent(in)  :: pos(2)
!local
   integer              :: i,j,inew,jnew

   allocate(get_MatAdd(k+1,k+1))
   get_MatAdd(pos(1),pos(2))=1d0/DetRat
   get_MatAdd(1:pos(1)-1,pos(2))=-matmul(N(1:pos(1)-1,:),Q)/DetRat
   get_MatAdd(pos(1)+1:k+1,pos(2))=-matmul(N(pos(1):k,:),Q)/DetRat
   get_MatAdd(pos(1),1:pos(2)-1)=-matmul(R,N(:,1:pos(2)-1))/DetRat
   get_MatAdd(pos(1),pos(2)+1:k+1)=-matmul(R,N(:,pos(2):k))/DetRat
   do i=1,k
      if(i.lt.pos(1))then 
         inew=i
      else
         inew=i+1
      endif
      do j=1,k
         if(j.lt.pos(2))then
            jnew=j
         else
            jnew=j+1
         endif
         get_MatAdd(inew,jnew)=N(i,j)+DetRat*get_MatAdd(inew,pos(2))*get_MatAdd(pos(1),jnew)
      enddo
   enddo
end function get_MatAdd

!===============================================================================
pure real(KINDR) function get_DetRatRem(N,k,pos)
!===============================================================================
!input
   integer, intent(in)     :: k,pos(2)
   real(KINDR), intent(in) :: N(k,k)

   get_DetRatRem=N(pos(1),pos(2))
end function get_DetRatRem

!===============================================================================
function get_MatRem(N,k,pos)
!===============================================================================
real(KINDR),pointer        :: get_MatRem(:,:)
!input
   integer                 :: k,pos(2)
   real(KINDR)             :: N(k,k),test(0:k+1,0:k+1)
!local
   integer                 :: i,j,iold,jold

   allocate(get_MatRem(k-1,k-1))
   
   do i=1,k-1
      if(i.lt.pos(1))then 
         iold=i
      else
         iold=i+1
      endif
      do j=1,k-1
         if(j.lt.pos(2))then
            jold=j
         else
            jold=j+1
         endif
         test(iold,pos(2))=N(iold,pos(2))
         get_MatRem(i,j)=N(iold,jold)-N(iold,pos(2))*N(pos(1),jold)/N(pos(1),pos(2))
      enddo
   enddo
end function get_MatRem


!===============================================================================
end module MMatrixUpdate
!===============================================================================

#ifdef MatrixUpdate_Test

!===============================================================================
program Prog_MatrixUpdate
!===============================================================================
use MParameters
use MMatrixUpdate

integer,parameter       :: MatSize=5
integer                 :: pos(MatSize,2),i,j,k,i1,j1
real(KINDR),pointer     :: Mat(:,:)=>null(),MatNew(:,:)=>null()
real(KINDR),pointer     :: InvMatExact(:,:),InvMatSherMorr(:,:),InvMatSherMorrOld(:,:)=>null()
real(KINDR)             :: DetExact,DetExactOld,DetRatSherMorr,DetOrigMatExact,DetOrigMatExactOld
real(KINDR)             :: FullMat(Matsize,Matsize)
pos=reshape((/1,1,1,2,2,1,1,3,2,4/),(/5,2/))
FullMat=reshape((/1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,11d0,12d0,13d0,14d0,&
   15d0,16d0,17d0,18d0,19d0,20d0,21d0,22d0,23d0,24d0,25d0/),(/5,5/))
!do i=1,MatSize
!do j=1,MatSize
!FullMat(i,j)=dble(rand())
!enddo
!enddo
do i=1,MatSize
write(*,*)FullMat(i,:)
enddo
DetExactOld=1d0
DetOrigMatExactOld=1d0
do k=1,MatSize
write(stdout,*)"================================================================"
write(stdout,*)"            adding an additional row and column at pos"
write(*,*)pos(k,:)
allocate(MatNew(k,k))
allocate(InvMatExact(k,k))
MatNew=0d0
MatNew(1:pos(k,1)-1,1:pos(k,2)-1)=Mat(1:pos(k,1)-1,1:pos(k,2)-1)
MatNew(pos(k,1)+1:k,1:pos(k,2)-1)=Mat(pos(k,1):k-1,1:pos(k,2)-1)
MatNew(1:pos(k,1)-1,pos(k,2)+1:k)=Mat(1:pos(k,1)-1,pos(k,2):k-1)
MatNew(pos(k,1)+1:k,pos(k,2)+1:k)=Mat(pos(k,1):k-1,pos(k,2):k-1)
MatNew(pos(k,1),pos(k,2))=FullMat(k,k)
MatNew(pos(k,1),1:pos(k,2)-1)=FullMat(k,1:pos(k,2)-1)
MatNew(pos(k,1),pos(k,2)+1:k)=FullMat(k,pos(k,2):k-1)
MatNew(1:pos(k,1)-1,pos(k,2))=FullMat(1:pos(k,1)-1,k)
MatNew(pos(k,1)+1:k,pos(k,2))=FullMat(pos(k,1):k-1,k)
write(stdout,*)"================================================================"
write(stdout,*)"            printing matrix"
write(stdout,*)"================================================================"
do i=1,k
write(*,*)MatNew(i,:)
enddo
write(stdout,*)"================================================================"
write(stdout,*)"        calculating the inverse and determinant exact"
InvMatExact=MatNew
call get_MatDetFull(k,InvMatExact,DetExact)
!write(*,*)FullMat(1:k,1:k)
!call get_DetFull(k,FullMat(1:k,1:k),DetOrigMatExact)
write(stdout,*)"================================================================"
write(stdout,*)"   calculating the inverse and determinant fast update"
DetRatSherMorr=get_DetRatAdd(InvMatSherMorr,k-1,FullMat(k,1:k-1),FullMat(1:k-1,k),FullMat(k,k))
InvMatSherMorr=>get_MatAdd(InvMatSherMorrOld,k-1,FullMat(k,1:k-1),FullMat(1:k-1,k),DetRatSherMorr,pos(k,:))
write(stdout,*)"================================================================"
write(stdout,*)"================================================================"
write(stdout,*)"     printing the determinant fast update and exact"
write(stdout,*)"================================================================"
write(stdout,*)"full: ",DetExact/DetExactOld
write(stdout,*)"fast update: ",DetRatSherMorr
!write(stdout,*)"full original: ",DetOrigMatExact/DetOrigMatExactOld
write(stdout,*)"================================================================"
write(stdout,*)"     printing inverse matrix fast update and exact"
write(stdout,*)"================================================================"
write(stdout,*)"Full: "
do i=1,k
write(stdout,*)InvMatExact(i,:)
enddo
write(stdout,*)"fast update: "
do i=1,k
write(stdout,*)InvMatSherMorr(i,:)
enddo
if(associated(Mat))deallocate(Mat)
allocate(Mat(k,k))
Mat=MatNew
deallocate(InvMatExact)
deallocate(MatNew)
DetExactOld=DetExact
DetOrigMatExactOld=DetOrigMatExact
if(associated(InvMatSherMorrOld))deallocate(InvMatSherMorrOld)
InvMatSherMorrOld=>InvMatSherMorr
enddo
!removing elements
do k=MatSize-1,1,-1
write(stdout,*)"================================================================"
write(stdout,*)"           removing a row and column at pos"
write(*,*)pos(k+1,:)
allocate(MatNew(k,k))
allocate(InvMatExact(k,k))
MatNew=0d0
MatNew(1:pos(k+1,1)-1,1:pos(k+1,2)-1)=Mat(1:pos(k+1,1)-1,1:pos(k+1,2)-1)
MatNew(pos(k+1,1):k,1:pos(k+1,2)-1)=Mat(pos(k+1,1)+1:k+1,1:pos(k+1,2)-1)
MatNew(1:pos(k+1,1)-1,pos(k+1,2):k)=Mat(1:pos(k+1,1)-1,pos(k+1,2)+1:k+1)
MatNew(pos(k+1,1):k,pos(k+1,2):k)=Mat(pos(k+1,1)+1:k+1,pos(k+1,2)+1:k+1)
write(stdout,*)"================================================================"
write(stdout,*)"            printing matrix"
write(stdout,*)"================================================================"
do i=1,k
write(*,*)MatNew(i,:)
enddo
write(stdout,*)"================================================================"
write(stdout,*)"        calculating the inverse and determinant exact"
InvMatExact=MatNew
call get_MatDetFull(k,InvMatExact,DetExact)
!call get_DetFull(k,FullMat(1:k,1:k),DetOrigMatExact)
write(stdout,*)"================================================================"
write(stdout,*)"   calculating the inverse and determinant fast update"
DetRatSherMorr=get_DetRatRem(InvMatSherMorrOld,k+1,pos(k+1,:))
InvMatSherMorr=>get_MatRem(InvMatSherMorrOld,k+1,pos(k+1,:))
write(stdout,*)"================================================================"
write(stdout,*)"     printing determinant fast update and exact"
write(stdout,*)"================================================================"
write(stdout,*)"full: ",DetExact/DetExactOld
write(stdout,*)"fast update: ",DetRatSherMorr
!write(stdout,*)"full original: ",DetOrigMatExact/DetOrigMatExactOld
write(stdout,*)"================================================================"
write(stdout,*)"     printing inverse matrix fast update and exact"
write(stdout,*)"================================================================"
write(stdout,*)"Full: "
do i=1,k
write(stdout,*)InvMatExact(i,:)
enddo
write(stdout,*)"fast update: "
do i=1,k
write(stdout,*)InvMatSherMorr(i,:)
enddo
if(associated(Mat))deallocate(Mat)
allocate(Mat(k,k))
Mat=MatNew
deallocate(InvMatExact)
deallocate(MatNew)
DetExactOld=DetExact
!DetOrigMatExactOld=DetOrigMatExact
if(associated(InvMatSherMorrOld))deallocate(InvMatSherMorrOld)
InvMatSherMorrOld=>InvMatSherMorr
enddo
end program Prog_MatrixUpdate

#endif
