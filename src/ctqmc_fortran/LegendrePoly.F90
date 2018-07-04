!> Small module for serving a set of legendre coefficients
module LegendrePoly
    use MParameters, only: KINDR
    implicit none

    !> The convention for fermionic matsubara frequencies. MAT_OFFSET marks
    !! the first n for which the fermionic Matsubara frequency \f$ \omega_n \f$
    !! is positive, i.e., when you use
    !!
    !!   - \f$ \omega_n = (2n + 1)\pi\tau/\beta \f$, then set \c MAT_OFFSET=0
    !!   - \f$ \omega_n = (2n - 1)\pi\tau/\beta \f$, then set \c MAT_OFFSET=1
    !!
    !! Note that this alters the required dimension of the array passed to the
    !! \c matsubara_phase functions
    !!
    integer, parameter :: MAT_OFFSET = 0

    real(KINDR), parameter, public :: PI = acos(-1.0_KINDR), TWOPI = 2*PI
contains
    !> Computes the values for a set of legendre polynomials for a given x.
    subroutine legendrep(x, pl)
        real(KINDR), intent(in)  :: x
        real(KINDR), intent(out) :: pl(0:)
        integer :: i

        ! sanity checks
        if(abs(x) > 1._KINDR) stop '[legendrep] x must be between -1 and 1'
        if(size(pl) < 2) stop '[legendrep] pl must be at least of size 0:1'

        ! Check for special values -1,1 first
        if(x == 1._KINDR) then
            pl(:) = 1._KINDR
        elseif(x == -1._KINDR) then
            forall(i = 0:size(pl)) pl(i) = real(m1pow(i), KINDR)
        else
            ! Employ Bonnet´s recursion formula for the Legendre polynomials as
            ! forward recurrence:
            !   l P_l = (2l-1) x P_{l-1} - (l-1) P_{l-2};  P_0 = 1;  P_1 = x
            ! Should be stable according to Acton - "Num. Meth. that work" §16
            pl(0) = 1._KINDR
            pl(1) = x
            do i = 2,ubound(pl,1)
                pl(i) = ((2*i-1)*x*pl(i-1) - (i-1)*pl(i-2))/i
            enddo
        endif
    end subroutine

    !> Computes a set of complex phase factors exp(iwt) for bosonic Matsubara
    !! frequencies
    !!
    !! This function is about ten times faster than computation by hand, but
    !! suffers from somewhat poor stability.
    subroutine bmatsubara_phase(tauoverbeta, ph)
        real(KINDR), intent(in) :: tauoverbeta
        complex(KINDR), intent(out) :: ph(0:)

        integer :: i
        real(KINDR) :: twoc

        if(ubound(ph,1) < 2 .or. lbound(ph,1) > 0 .or. -lbound(ph,1) > ubound(ph,1)) &
            stop '[bmatsubara_phase] Array malformed, should be (-i<=j<1):(i>1)'

        ! This basically uses the complex form of the Chebyshev method (since
        ! this complex phase factors are related to the Chebyshev polynomials):
        !       exp[i(n+1)x] = 2 cos(x) exp(inx) - exp(i(n-1)x)
        ! as forward recursion. This is guaranteed to be stable.
        ph(0) = cmplx(1, kind=kindr)
        ph(1) = exp(cmplx(0,TWOPI*tauoverbeta, kind=kindr))   ! exp(i2pi*t/b)
        twoc = 2*real(ph(1))
        do i = 2, ubound(ph,1)
            ph(i) = twoc*ph(i-1) - ph(i-2)
        enddo
        !forall(i = lbound(ph,1):-1) ph(i) = conjg(ph(-i))
    end subroutine

    !> Computes a set of complex phase factors exp(iwt) for fermionic Matsubara
    !! frequencies
    subroutine fmatsubara_phase(tauoverbeta, ph)
        real(KINDR), intent(in) :: tauoverbeta
        complex(KINDR), intent(out) :: ph(MAT_OFFSET:)

        integer :: i
        real(KINDR) :: twoc

        if(ubound(ph,1) < MAT_OFFSET+2 .or. lbound(ph,1) > MAT_OFFSET &
                    .or. 2*MAT_OFFSET-1-lbound(ph,1) > ubound(ph,1)) &
            stop '[bmatsubara_phase] Array malformed, should be (-i<=j<1):(i>1)'

        ! This basically uses the complex form of the Chebyshev method (since
        ! this complex phase factors are related to the Chebyshev polynomials):
        !       exp[i(n+1)x] = 2 cos(x) exp(inx) - exp(i(n-1)x)
        ! as forward recursion. This is guaranteed to be stable and also
        ! extends to non-integer n.
        twoc = 2*cos(TWOPI*tauoverbeta)
        ph(MAT_OFFSET)   = exp(cmplx(0,PI*tauoverbeta, kind=kindr))
        ph(MAT_OFFSET+1) = twoc*ph(MAT_OFFSET) - conjg(ph(MAT_OFFSET))
        do i = MAT_OFFSET+2, ubound(ph,1)
            ph(i) = twoc*ph(i-1) - ph(i-2)
        enddo
        !forall(i = lbound(ph,1):MAT_OFFSET-1) ph(i) = conjg(ph(2*MAT_OFFSET-1-i))
    end subroutine

    !> Returns (-1)^n in a fast way, i.e., potentially much faster than (-1)**n
    pure integer function m1pow(n)
        integer, intent(in) :: n
        select case (modulo(n, 2))  ! also positive for n<0
            case (0);     m1pow = 1
            case default; m1pow = -1
        end select
    end function

end module LegendrePoly

#ifdef LegendrePoly_Test

!> Test programme provided for testing the Legendre coefficients
program Prog_LegendrePoly
    use LegendrePoly
    use MParameters, only: KINDR, PI

    ! Values from Mathematica
    real(KINDR), parameter :: lagm013(0:100) = (/ 1.000000,-0.013000,-0.499747,&
          0.019495, 0.374366,-0.024356,-0.311391, 0.028394, 0.271775,-0.031913,&
         -0.243810, 0.035063, 0.222619,-0.037931,-0.205766, 0.040574, 0.191884,&
         -0.043029,-0.180136, 0.045325, 0.169980,-0.047481,-0.161047, 0.049513,&
          0.153077,-0.051433,-0.145877, 0.053250, 0.139308,-0.054974,-0.133259,&
          0.056609, 0.127646,-0.058162,-0.122401, 0.059637, 0.117472,-0.061039,&
         -0.112815, 0.062369, 0.108393,-0.063632,-0.104177, 0.064829, 0.100143,&
         -0.065963,-0.096270, 0.067036, 0.092540,-0.068049,-0.088937, 0.069005,&
          0.085450,-0.069904,-0.082067, 0.070747, 0.078778,-0.071536,-0.075576,&
          0.072272, 0.072453,-0.072956,-0.069403, 0.073588, 0.066420,-0.074169,&
         -0.063500, 0.074701, 0.060638,-0.075183,-0.057831, 0.075618, 0.055076,&
         -0.076004,-0.052369, 0.076343, 0.049708,-0.076635,-0.047091, 0.076882,&
          0.044516,-0.077083,-0.041981, 0.077239, 0.039485,-0.077351,-0.037026,&
          0.077419, 0.034604,-0.077444,-0.032217, 0.077426, 0.029865,-0.077366,&
         -0.027546, 0.077264, 0.025261,-0.077121,-0.023008, 0.076937, 0.020788/)
    real(KINDR), parameter :: lag05429(0:100) = (/1.000000, 0.542900,-0.057889,&
         -0.414314,-0.350212,-0.010783, 0.281111, 0.292670, 0.051949,-0.206879,&
         -0.260152,-0.081561, 0.153604, 0.235655, 0.104104,-0.110676,-0.214015,&
         -0.121376, 0.073995, 0.193218, 0.134256,-0.041713,-0.172416,-0.143240,&
          0.012942, 0.151282, 0.148659, 0.012745,-0.129757,-0.150767,-0.035543,&
          0.107934, 0.149796, 0.055520,-0.085992,-0.145971,-0.072690, 0.064165,&
          0.139531, 0.087041,-0.042716,-0.130733,-0.098561, 0.021919, 0.119850,&
          0.107256,-0.002053,-0.107179,-0.113152,-0.016616, 0.093028, 0.116310,&
          0.033836,-0.077723,-0.116819,-0.049379, 0.061596, 0.114807, 0.063049,&
         -0.044983,-0.110434,-0.074681, 0.028218, 0.103891, 0.084147,-0.011629,&
         -0.095403,-0.091360,-0.004469, 0.085219, 0.096275, 0.019780,-0.073609,&
         -0.098887,-0.034031, 0.060863, 0.099234, 0.046976,-0.047282,-0.097396,&
         -0.058400, 0.033174, 0.093488, 0.068124,-0.018847,-0.087666,-0.076007,&
          0.004605, 0.080114, 0.081946, 0.009259,-0.071048,-0.085883,-0.022466,&
          0.060705, 0.087796, 0.034760,-0.049343,-0.087709,-0.045909, 0.037233/)
    real(KINDR), parameter :: LPREC = 1E-5, MPREC = 1E-20
    real(KINDR) :: legvalues(0:100) ! Output array

    real(KINDR), parameter :: FVAL = 0.99999999999
    integer, parameter :: MAXMAT = 1000000
    real :: startexpl, endexpl, endrec
    complex(KINDR), dimension(0:MAXMAT) :: fph, ftest
    complex(KINDR), dimension(0:MAXMAT) :: bph, btst
    integer :: i,b,e

    call legendrep(-0.013_KINDR, legvalues)
    write (0,"('Precision @x=-0.013: ',ES11.3,' to ',ES11.3)") &
            minval(abs(legvalues - lagm013)), maxval(abs(legvalues - lagm013))

    call legendrep(0.5429_KINDR, legvalues)
    write (0,"('Precision @x=-0.013: ',ES11.3,' to ',ES11.3)") &
            minval(abs(legvalues - lag05429)), maxval(abs(legvalues - lag05429))

    write (0,*) 'Fermionic'
    call cpu_time(startexpl)
    forall(i = 0:MAXMAT) ftest(i) = exp(cmplx(0_KINDR,(2*i+1)*PI*FVAL))
    call cpu_time(endexpl)
    call fmatsubara_phase(FVAL, fph)
    call cpu_time(endrec)

    do i = 0, int(log10(MAXMAT*1.0))
        b = 10**i
        e = min(10*b, MAXMAT)
        write (0,"('Precision in [',I7,':',I7,']: ',ES11.3,' to ',ES11.3)") &
                b,e,minval(abs(fph(b:e) - ftest(b:e))), maxval(abs(fph(b:e) - ftest(b:e)))
    enddo
    write (0,"('Explicit:',F9.3,'s; Non-explicit:',F9.3,'s')") &
            endexpl-startexpl, endrec-endexpl


    write (0,*) 'Bosonic'
    call cpu_time(startexpl)
    forall(i = 0:MAXMAT) btst(i) = exp(cmplx(0_KINDR,(2*i)*PI*FVAL))
    call cpu_time(endexpl)
    call bmatsubara_phase(FVAL, bph)
    call cpu_time(endrec)

    do i = 0, int(log10(MAXMAT*1.0))
        b = 10**i
        e = min(10*b, MAXMAT)
        write (0,"('Precision in [',I7,':',I7,']: ',ES11.3,' to ',ES11.3)") &
                b,e,minval(abs(bph(b:e) - btst(b:e))), maxval(abs(bph(b:e) - btst(b:e)))
    enddo
    write (0,"('Explicit:',F9.3,'s; Non-explicit:',F9.3,'s')") &
            endexpl-startexpl, endrec-endexpl

    write(0,*) 'Complete.'
end program


#endif
