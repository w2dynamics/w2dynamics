!____________________________________________________________________________
! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999

 module MMersenneTwister
! Default seed
    integer, parameter :: defaultsd = 4357
! Period parameters
    integer, parameter :: N = 624, N1 = N + 1

! the array for the state vector
    integer, save, dimension(0:N-1) :: mt
    integer, save                   :: mti = N1
    
! Overload procedures for saving and getting mt state
    interface mtsave
      module procedure mtsavef
      module procedure mtsaveu
    end interface
    interface mtget
      module procedure mtgetf
      module procedure mtgetu
    end interface

 contains

!Initialization subroutine
  subroutine sgrnd(seed)
    implicit none
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
    integer, intent(in) :: seed

    mt(0) = iand(seed,-1)
    do mti=1,N-1
      mt(mti) = iand(69069 * mt(mti-1),-1)
    enddo
!
    return
  end subroutine sgrnd

!Random number generator
  real(8) function grnd()
    implicit integer(a-z)

! Period parameters
    integer, parameter :: M = 397, MATA  = -1727483681
!                                    constant vector a
    integer, parameter :: LMASK =  2147483647
!                                    least significant r bits
    integer, parameter :: UMASK = -LMASK - 1
!                                    most significant w-r bits
! Tempering parameters
    integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544
    integer            :: kk 

    dimension mag01(0:1)
    data mag01/0, MATA/
    save mag01
!                        mag01(x) = x * MATA for x=0,1

    TSHFTU(y)=ishft(y,-11)
    TSHFTS(y)=ishft(y,7)
    TSHFTT(y)=ishft(y,15)
    TSHFTL(y)=ishft(y,-18)

    if(mti.ge.N) then
!                       generate N words at one time
      if(mti.eq.N+1) then
!                            if sgrnd() has not been called,
        call sgrnd( defaultsd )
!                              a default initial seed is used
      endif

      do kk=0,N-M-1
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      do kk=N-M,N-2
          y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
      enddo
      y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
      mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
      mti = 0
    endif

    y=mt(mti)
    mti = mti + 1 
    y=ieor(y,TSHFTU(y))
    y=ieor(y,iand(TSHFTS(y),TMASKB))
    y=ieor(y,iand(TSHFTT(y),TMASKC))
    y=ieor(y,TSHFTL(y))

    if(y .lt. 0) then
      grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
    else
      grnd=dble(y)/(2.0d0**32-1.0d0)
    endif

    return
  end function grnd

!State saving subroutines.
! Usage:  call mtsave( file_name, format_character )
!    or   call mtsave( unit_number, format_character )
! where   format_character = 'u' or 'U' will save in unformatted form, otherwise
!         state information will be written in formatted form.
  subroutine mtsavef( fname, forma )

!NOTE: This subroutine APPENDS to the end of the file "fname".

    character(*), intent(in) :: fname
    character, intent(in)    :: forma

    select case (forma)
      case('u','U')
       open(unit=10,file=trim(fname),status='UNKNOWN',form='UNFORMATTED', &
            position='APPEND')
       write(10)mti
       write(10)mt

      case default
       open(unit=10,file=trim(fname),status='UNKNOWN',form='FORMATTED', &
            position='APPEND')
       write(10,*)mti
       write(10,*)mt

    end select
    close(10)

    return
  end subroutine mtsavef

  subroutine mtsaveu( unum, forma )

    integer, intent(in)    :: unum
    character, intent(in)  :: forma

    select case (forma)
      case('u','U')
       write(unum)mti
       write(unum)mt

      case default
       write(unum,*)mti
       write(unum,*)mt

      end select

    return
  end subroutine mtsaveu

!State getting subroutines.
! Usage:  call mtget( file_name, format_character )
!    or   call mtget( unit_number, format_character )
! where   format_character = 'u' or 'U' will read in unformatted form, otherwise
!         state information will be read in formatted form.
  subroutine mtgetf( fname, forma )

    character(*), intent(in) :: fname
    character, intent(in)    :: forma

    select case (forma)
      case('u','U')
       open(unit=10,file=trim(fname),status='OLD',form='UNFORMATTED')
       read(10)mti
       read(10)mt

      case default
       open(unit=10,file=trim(fname),status='OLD',form='FORMATTED')
       read(10,*)mti
       read(10,*)mt

    end select
    close(10)

    return
  end subroutine mtgetf

  subroutine mtgetu( unum, forma )

    integer, intent(in)    :: unum
    character, intent(in)  :: forma

    select case (forma)
      case('u','U')
       read(unum)mti
       read(unum)mt

      case default
       read(unum,*)mti
       read(unum,*)mt

      end select

    return
  end subroutine mtgetu

 end module MMersenneTwister

#ifdef MersenneTwister_Test

 program Prog_MersenneTwister
! this main() outputs first 1000 generated numbers
!
    use MMersenneTwister
    implicit none
 
    integer, parameter      :: no=1000
    real(8), dimension(0:7) :: r
    integer j,k
!    real(8) grnd
!
!      call sgrnd(4357)
!                         any nonzero integer can be used as a seed
    do j=0,no-1
      r(mod(j,8))=grnd()
      if(mod(j,8).eq.7) then
        write(*,'(8(f9.6,'' ''))') (r(k),k=0,7)
      else if(j.eq.no-1) then
        write(*,'(8(f9.6,'' ''))') (r(k),k=0,mod(no-1,8))
      endif
    enddo

 end program Prog_MersenneTwister

#endif
