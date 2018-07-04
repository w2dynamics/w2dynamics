!compiler call
!ifort Parameters.o Nfft_worm.F90 -DFT_Test -free -implicitnone -warn all -warn noerrors -warn nounused -limf -lm -lnfft3 -lfftw3

module ft_worm
!Fourier transform library to use with qmc estimators (naive and nfft)
!x stores the tau data
!f stores the values to be transformed
!fhat stores the transformed values 
!(2 time the inteded dimension, because we lose every other frequency)

use MParameters
use iso_c_binding

#ifdef USE_NFFT
use nfft_static
use nfft_interface
#endif

contains
#ifdef USE_NFFT

!===============================================================================
!n-dimensional nfft fourier transform
!Note: x,f,f_hat require contiguous memory layout when mapped to C functions
subroutine ft_nd(x,f,f_hat,d)
!===============================================================================
   !input
   complex(c_double_complex),target :: f(:),f_hat(:)
   real(c_double),target            :: x(:)
   integer(c_int)                   :: d(:)
   !temp
   integer                                :: i,power,k   

   !gfortran workaround
   complex(c_double_complex),pointer :: f_ptr,f_hat_ptr
   real(c_double),pointer            :: x_ptr
   

   !nfft internals
   integer(c_int),parameter               :: m=6
   integer(c_int)                         :: nfft_flags,fftw_flags
   integer(c_int),pointer                 :: nn(:),n(:)
   type(nfft_plan)                        :: plan 

   !re-initialize to zero
   f_hat=cmplx(0d0,0d0, kind=c_double_complex)
   
   !nfft internals
   allocate(nn(size(d)),n(size(d)))
   
   do i=1,size(d)
      nn(i)=d(i)
      power=int(log10(dble(nn(i)))/log10(2d0))
      n(i)=2*(2**(1+power))
      !constraint m < fftw_length/2 -> we choose a minimal fftw_length of 2m+1
      if(n(i)<(2*m+1)) n(i)=2*m+1
   enddo
   
   !we do not add the flags for the x,f,f_hat allocation since we do that in the fortran part
   !while the nfft routines suggest and fftw out of place we prefer in place to save on resources
   nfft_flags=IOR(PRE_PHI_HUT,IOR(PRE_PSI,FFTW_INIT))
   fftw_flags=FFTW_MEASURE

   !for large transforms we want to use FFTW_MEASURE instead of FFTW_ESTIMATE
   call nfft_init_guru(plan,size(d),c_loc(nn(1)),size(f),c_loc(n(1)),m,nfft_flags,fftw_flags)
    
   x_ptr=>x(lbound(x,1))
   f_ptr=>f(lbound(f,1))
   f_hat_ptr=>f_hat(lbound(f_hat,1))

   plan%x=c_loc(x_ptr)
   plan%f=c_loc(f_ptr)
   plan%f_hat=c_loc(f_hat_ptr)


   if(iand(plan%nfft_flags, PRE_ONE_PSI)/=0) then
      call nfft_precompute_one_psi(plan)
   endif       
   
   !do the adjoint transfrom
   call nfft_check(plan) 
   call nfft_adjoint(plan)
   call nfft_finalize(plan)

   deallocate(nn,n)
    
end subroutine
#endif

#ifndef USE_NFFT
!===============================================================================
!n-dimensional naive fourier transform
!temporary vector and omega vector may not be necessary
subroutine ft_nd(x,f,f_hat,d)
!===============================================================================
   !input
   complex(c_double_complex)      :: f(:),f_hat(:)
   real(c_double)                 :: x(:)
   integer(c_int)                 :: d(:)

   integer(c_int),allocatable     :: ind(:),r(:)
   integer                        :: i,j,k
   real(c_double),allocatable     :: omega(:),tmp(:)
   real(c_double)                 :: pi
   pi=dacos(-1d0)

   allocate(omega(size(d)*size(f_hat)))
   
   !nd mesh from 1d representation
   allocate(ind(size(d)),r(size(d)))
   do k=0,size(f_hat)-1
      ind(1)=k/product(d(2:))
      r(1)=mod(k,product(d(2:)))
      do j=2,size(d)
        ind(j)=r(j-1)/product(d(j+1:))
        r(j)=mod(r(j-1),product(d(j+1:)))
      enddo
      do i=1,size(d)
         omega(size(d)*k+i)=pi*(-(d(i)/2)+ind(i))
      enddo
   enddo
      
   !re-initialize to zero
   f_hat=cmplx(0d0, 0d0, kind=c_double_complex)
   
   allocate(tmp(size(f_hat)))
   do j=0,size(f)-1
      tmp=0d0
      do i=1,size(d)
         tmp=tmp+omega(i::size(d))*2d0*x(size(d)*j+i)
      enddo
      f_hat=f_hat+exp(cmplx(0d0,tmp,kind=c_double_complex))*f(j+1)
   enddo
   
   deallocate(omega,tmp,ind,r)
  
end subroutine
#endif

end module ft_worm

#ifdef FT_Test
program ft_test
use ft_worm
implicit none

real(c_double), allocatable             :: x(:)
complex(c_double_complex), allocatable  :: f(:),f_hat(:)
integer(c_int), allocatable             :: d(:)
integer                            :: i,j,base,inlen,dims
real(c_double),parameter           :: pi=dacos(-1d0)

dims=2

!output vectors
allocate(d(dims))
do j=1,dims
   d(j)=200
enddo
allocate(f_hat(product(d)))

!input vectos
inlen=20000
allocate(x(dims*inlen),f(inlen))

!mesh
call random_number(x)
x=x-0.5d0
f=cmplx(0d0,0d0)
do i=0,size(f)-1
   if(abs(x(2*i+1))<0.25 .and. abs(x(2*i+2))<0.25) then
     f(i+1)=1000d0
   endif
enddo

call ft_nd(x,f,f_hat,d)

f_hat=f_hat/dble(inlen)

if (dims .eq. 2) then
do i=0,d(1)-1
   do j=1,d(2)
      write(2,*) i+1,j,real(f_hat(d(2)*i+j)),aimag(f_hat(d(2)*i+j))
      !cut
      if(j.eq.d(2)/2) then
         write(3,*) i+1,real(f_hat(d(2)*i+j)),aimag(f_hat(d(2)*i+j))
      endif 
   enddo
   write(2,*) " "
enddo
endif

if (dims .eq. 1) then
do i=0,d(1)-1
      write(2,*) i+1,real(f_hat(i+1)),aimag(f_hat(i+1))
enddo
endif

deallocate(x,f,f_hat) 
end program
#endif
