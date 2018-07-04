!===============================================================================
!> This module generates Gaunt coefficients, 3jm symbols, Clebsch 
!! Gordon coefficients and factorials of
!! integer numbers and integer numbers+1/2. The factorials are stored in an 
!! array for a fast look up if it was calculated before.
!!
!! @author Nicolaus Parragh
!!
!! @todo The 3jm symbols could also be stored in a data structure. Currently
!! it does not really matter if things are fast in this module so there is no
!! real reason to implement this.
!!
!! @note There is a different algorithm to calculate the Gaunt coefficients 
!! developed by Xu (see 
!! MATHEMATICS OF COMPUTATION
!! Volume 65, Number 216
!! October 1996, Pages 1601–1612)
!! which is supposed to be more stable and faster than the algorithm currently
!! implemented which uses the 3jm symbols. But as already mentioned in the todo
!! list the speed really is no issue here and I also could not reproduce the
!! numerical instability Xu claims the approach implemented here has. The test
!! program contained in this file calculates some quantities also mentioned in
!! Xu paper.
module MAngularMomentum
!===============================================================================
use MParameters

!> This array stores the factorial up to 500! (the factorials of N.5! (N element
!! of natural numbers) can also be calculated) which have been calculated before.
   real(KINDR),save                 :: fact(0:1000)=0d0

!> See also the todo list for this module. It might be a good idea to store the
!! 3jm symbols which have alreay been calculated before. But this needs testing
!! and speed really is not the biggest issue in this module.
   type :: T3jm
!> The parameters specifying the 3jm symbol.
      real(KINDR)                   :: j1,j2,j3,m1,m2,m3 
!> The value for the parameters specified in the variables above.
      real(KINDR)                   :: value
!> The next 3jm symbol
      type(T3jm),pointer            :: next
   end type T3jm

   complex(KINDC),parameter,private :: re=cmplx(1.D0, 0.D0, kind=kindc), im=cmplx(0.D0, 1.D0, kind=kindc)
!> The Pauli spin matrices
   complex(KINDC),parameter,dimension(3,2,2) :: PauliSM=reshape(&
      (/0*re,0*re,re,re,im,0*re,re,-im,0*re,0*re,0*re,-re/),(/3,2,2/))
!> The spin creation/annihilation operators
   real(KINDR),parameter,dimension(2,2,2) :: Spm=reshape(&
      (/0,0,1,0,0,1,0,0/),(/2,2,2/))

!> This would be the pointer to the first 3jm symbol if it would get stored.
!! Currently rather useless...
!   type(T3jm),pointer,save,private  :: first_PT3jm=>null()

contains  
!> This functions returns the factorial of N and N.5 (N element 
!! of natural numbers). It stores factorials which have already been calculated
!! once in an array.
!! Formula for factorial taken from wikipedia.
!! @param N the number for which the factorial should be calculated. Either a
!! natural number or a natural number + 1/2.
!! @return the factorial of the number N.
!===============================================================================
real(KINDR) function get_fact(N)
!===============================================================================
!input
   real(KINDR)                      :: N
!local
   real(KINDR)                      :: x

   if(N<0d0)then
      write(stdout,*)"Error in MAngularMomentum::get_fact(N) N<0d0"
      stop
   endif
   if(fact(int(2d0*N))==0d0)then
      x=N
      get_fact=1d0
      do while(x.gt.1d0)
         get_fact=get_fact*x
         x=x-1d0
      enddo
      if(abs(x-0.5d0)<1d-10) get_fact=get_fact*sqrt(acos(-1d0))*0.5d0
      fact(int(2d0*N))=get_fact
   else
      get_fact=fact(int(2d0*N))
   endif
end function get_fact

!> This function calculates the 3jm symbol for the six defining numbers
!! j1,j2,j3,m1,m2,m3
!! The formula to calculate the 3jm symbol is taken http://dlmf.nist.gov/34.2.
!! @param j1 a non-negative integer, or a half-odd positive integer.
!! @param j2 a non-negative integer, or a half-odd positive integer.
!! @param j3 a non-negative integer, or a half-odd positive integer.
!! @param m1 |m1| <= j1
!! @param m2 |m2| <= j2
!! @param m3 |m3| <= j3
!! @return The value of the 3jm symbol for given parameters.
!===============================================================================
real(KINDR) function get_3jm(j1,j2,j3,m1,m2,m3)
!===============================================================================
!input
   real(KINDR),intent(in)           :: j1,j2,j3,m1,m2,m3
!local
   integer                          :: ik,lb,ub
   !type(T3jm),pointer               :: PT3jm
   real(KINDR)                      :: sum3jm,k

   get_3jm=0d0
   if(abs(m1+m2+m3).gt.1d-10)return
   if(abs(m1).gt.j1.or.abs(m2).gt.j2.or.abs(m3).gt.j3)return
   if(abs(j1-j2).gt.j3.or.j3.gt.j1+j2)return
   
!   PT3jm=>first_PT3jm
!   if(associated(first_PT3jm))then
!      do
!         if(eq(j1,PT3jm%j1).and.eq(j2,PT3jm%j2).and.eq(j3,PT3jm%j3)&
!            .and.eq(m1,PT3jm%m1).and.eq(m2,PT3jm%m2).and.eq(m3,PT3jm%m3))then
!            get_3jm=PT3jm%value
!            return
!         endif
!         if(.not.associated(PT3jm%next))exit
!         PT3jm=>PT3jm%next
!      enddo
!   endif
! j1-j2-m3 is always an integer number, even for the cases where half-odds are allowed
   get_3jm=(-1d0)**(int(j1-j2-m3))&
      *sqrt(get_fact(j1+j2-j3)*get_fact(j1-j2+j3)*get_fact(-j1+j2+j3)/get_fact(j1+j2+j3+1d0))&
      *sqrt(get_fact(j1+m1)*get_fact(j1-m1)*get_fact(j2+m2)*get_fact(j2-m2)*get_fact(j3+m3)*get_fact(j3-m3))
   lb=int(max(j2-j3-m1,j1-j3+m2,0d0))
   ub=int(min(j1+j2-j3,j1-m1,j2+m2))
!   write(*,*)'lb=',lb,j1-j2+m3,j1-m1,' ub=',ub,j2+j3-m1,j3-j1+j2,j3-m3
   sum3jm=0d0
   do ik=lb,ub
      k=dble(ik)
!      write(*,*)"factpars: ",k,j2+j3-m1-k,j1-m1+k,j3-j1+j2-k,j3-m3-k,k+j1-j2+m3
      sum3jm=sum3jm+(-1)**(ik)&
         /(get_fact(k)*get_fact(j1+j2-j3-k)*get_fact(j1-m1-k)*get_fact(j2+m2-k)&
         *get_fact(j3-j2+m1+k)*get_fact(j3-j1-m2+k))

   enddo
!   if(abs(sum3jm)<1d-10)sum3jm=1d0
   get_3jm=get_3jm*sum3jm
!   PT3jm=>first_PT3jm
!   allocate(first_PT3jm)
!   first_PT3jm%j1=j1
!   first_PT3jm%j2=j2
!   first_PT3jm%j3=j3
!   first_PT3jm%m1=m1
!   first_PT3jm%m2=m2
!   first_PT3jm%m3=m3
!   first_PT3jm%value=get_3jm
!   first_PT3jm%next=>PT3jm
end function get_3jm

!> This function calculates the Gaunt coefficient using the 3jm symbol.
!! The formula is taken from the Xu paper which is referenced in the 
!! documentation of the module.
!! @param j1
!! @param j2
!! @param j3
!! @param m1 In the papers by Xu this is an integer.
!! @param m2 In the papers by Xu this is an integer.
!! @return The Gaunt coefficient for the given parameters
!===============================================================================
real(KINDR) function get_gaunt(j1,j2,j3,m1,m2)
!===============================================================================
!input
   real(KINDR),intent(in)         :: j1,j2,j3,m1,m2

   get_gaunt=0d0
   if(((j1+m1).lt.0d0).or.((j2+m2).lt.0d0).or.((j3+m1+m2).lt.0d0))return
   if(((j1-m1).lt.0d0).or.((j2-m2).lt.0d0).or.((j3-m1-m2).lt.0d0))return
   
   get_gaunt=(-1d0)**(int(m1+m2))*(2d0*j3+1d0)*sqrt(get_fact(j1+m1)*get_fact(j2+m2)*get_fact(j3-m1-m2))&
      /sqrt(get_fact(j1-m1)*get_fact(j2-m2)*get_fact(j3+m1+m2))&
      *get_3jm(j1,j2,j3,0d0,0d0,0d0)*get_3jm(j1,j2,j3,m1,m2,-m1-m2)
end function get_gaunt

!> This function calculates the Clebsch Gordan coefficient for the given
!! parameters. The formula is taken from Wikipedia.
!! @param j1 a non-negative integer, or a half-odd positive integer.
!! @param j2 a non-negative integer, or a half-odd positive integer.
!! @param j3 a non-negative integer, or a half-odd positive integer.
!! @param m1 |m1| <= j1
!! @param m2 |m2| <= j3
!! @param m3 |m3| <= j3
!! @return The Clebsch Gordan coefficient for the given parameters.
!===============================================================================
real(KINDR) function get_clebschgordan(j1,j2,j3,m1,m2,m3)
!===============================================================================
!input
   real(KINDR),intent(in)         :: j1,j2,j3,m1,m2,m3
! The phase factor turns out to be truly real, by considering since |j1-j2| <= j3 <= j1+j2
   get_clebschgordan=(-1d0)**(int(j1-j2+m3))*sqrt(2d0*j3+1d0)*get_3jm(j1,j2,j3,m1,m2,-m3)
end function get_clebschgordan

!!This subroutine calculates the legendre polynomial of order l
!!at point x recursively
!===============================================================================
subroutine legendre_poly(lp1,x,Plp1x) 
!===============================================================================
   integer, intent(in)      :: lp1
   real(kindr), intent(in)  :: x
   real(kindr), intent(out) :: Plp1x
   integer                  :: n,l
   real(kindr)              :: Plx,Plm1x
    
   if(x.eq.-1d0.or.x.eq.1d0)then
      Plp1x=(x)**lp1
      return
   endif
   Plx=0d0
   Plm1x=0d0 
   do n=0,lp1
      if(n.eq.0)then
         Plp1x=1d0
      elseif(n.eq.1)then
         Plp1x=x
      else
         l=n-1
         Plp1x=(2d0*dble(l)+1d0)/dble(l+1)*x*Plx-dble(l)/dble(l+1)*Plm1x
      endif
      Plm1x=Plx
      Plx=Plp1x
   enddo
  
   return
end subroutine legendre_poly

!===============================================================================
end module MAngularMomentum
!===============================================================================

#ifdef AngularMomentum_Test
!===============================================================================
!>This program tests the AngularMomentum model
program AngularMomentum_Prog
!===============================================================================
use MAngularMomentum

!local
   integer                    :: j1,j2,j3,m1,m2,m3,i
   real(KINDR)                :: rj1,rj2,rj3,rm1,rm2,rm3

write(stdout,*)"================================================================"
write(stdout,*)"           printing the 0,1/2,1,...,10 factorials"
write(stdout,*)"================================================================"
do i=0,20
   write(stdout,*)dble(i)/2d0,get_fact(dble(i)/2d0)
enddo
write(stdout,*)"================================================================"
write(stdout,*)" printing wigner 3jm symbol j1=0,1 j2=0,1 j3=0,1"
write(stdout,*)"================================================================"
do j1=0,1
   do j2=0,1
      do j3=0,j1+j2
         do m1=-j1,j1
            do m2=-j2,j2
               do m3=-j3,j3
                  rj1=dble(j1)
                  rj2=dble(j2)
                  rj3=dble(j3)
                  rm1=dble(m1)
                  rm2=dble(m2)
                  rm3=dble(m3)
                  write(*,*)j1,j2,j3,m1,m2,m3,get_3jm(rj1,rj2,rj3,rm1,rm2,rm3)
               enddo
            enddo
         enddo
      enddo
   enddo
enddo

write(stdout,*)"================================================================"
write(stdout,*)"  printing clebsch gordan coefficients for j1=1/2 j2=1/2"
write(stdout,*)"================================================================"
write(stdout,*) get_clebschgordan(1d0/2d0,1d0/2d0,1d0,1d0/2d0,1d0/2d0,1d0)
write(stdout,*) get_clebschgordan(1d0/2d0,1d0/2d0,1d0,1d0/2d0,-1d0/2d0,0d0)
write(stdout,*) get_clebschgordan(1d0/2d0,1d0/2d0,0d0,1d0/2d0,-1d0/2d0,0d0)
write(stdout,*) get_clebschgordan(1d0/2d0,1d0/2d0,1d0,-1d0/2d0,-1d0/2d0,-1d0)
write(stdout,*)"================================================================"
write(stdout,*)"  printing clebsch gordan coefficients for j3=3"
write(stdout,*)"================================================================"
j3=3

do j1=0,j3
do j2=0,j3
do m1=-j1,j1
   do m2=-j2,j2
         do m3=-j3,j3
            if(m1+m2.ne.m3)cycle
            write(stdout,*)j1,j2,m1,m2,j3,m3,get_clebschgordan(dble(j1),dble(j2),dble(j3),dble(m1),dble(m2),dble(m3))
         enddo
      enddo
   enddo
enddo
enddo
write(stdout,*)"================================================================"
write(stdout,*)" printing gaunt coefficients to compare with xu paper"
write(stdout,*)"================================================================"
m1=1
m2=-1
do j1=1,20
   j2=j1
   j3=j1+j2
   rj1=dble(j1)
   rj2=dble(j2)
   rj3=dble(j3)
   rm1=dble(m1)
   rm2=dble(m2)
   rm3=dble(m3)
   write(*,*)j1,j2,j3,m1,m2,get_gaunt(rj1,rj2,rj3,rm1,rm2)
!i could not find the difference he claims but i am also using a different algorithm to
!calculate the 3jm symbol....
enddo
end program AngularMomentum_Prog

#endif
