
module nfft_fallback
    !< Fallback interface in case no NFFT is there
    use nfft_static
    use iso_c_binding
    use LegendrePoly, only: PI
    interface
        function my_malloc(bytes) bind(C,name="malloc") result(ptr)
            use iso_c_binding
            integer(c_size_t), value :: bytes
            type(c_ptr) :: ptr
        end function

        subroutine my_free(ptr) bind(C,name="free")
            use iso_c_binding
            type(c_ptr), value :: ptr
        end subroutine
    end interface
contains
    subroutine warn_fallback(func)
        character(len=*), intent(in) :: func
        logical, save :: warn = .true.
        if (warn) then
            warn = .false.
            write(0,"('[',A,'] WARNING: Using slow fallback for NFFT')") func
        endif
    end subroutine

    subroutine nfft_init_1d(ths, n1, m)
        type(nfft_plan), intent(inout) :: ths
        integer(c_size_t) :: n1, m
        real(c_double) :: cdouble
        complex(c_double_complex) ::cdoublecomplex

        call warn_fallback('nfft_init_1d')
        ths%f_hat = my_malloc(n1 * c_sizeof(cdoublecomplex))
        ths%f = my_malloc(m * c_sizeof(cdoublecomplex))
        ths%x = my_malloc(m * c_sizeof(cdouble))
        ths%nn_total = n1
        ths%mm_total = m
    end subroutine

    subroutine nfft_init_2d(ths, n1, n2, m)
        type(nfft_plan), intent(inout) :: ths
        integer(c_int) :: n1, n2, m

        call warn_fallback('nfft_init_2d')
        stop '[nfft_init_2d] Not implemented in fallback interface (1D only'
    end subroutine

    subroutine nfft_finalize(ths)
        type(nfft_plan), intent(inout) :: ths

        call warn_fallback('nfft_finalize')
        call my_free(ths%f_hat)
        call my_free(ths%f)
        call my_free(ths%x)
    end subroutine

    subroutine nfft_precompute_one_psi(ths)
        type(nfft_plan), intent(inout) :: ths

        call warn_fallback('nfft_precompute_one_psi')
    end subroutine

    subroutine nfft_adjoint(ths)
        !< Computes the adjoint fourier transform:
        !! \f[
        !!     \hat f(n) = \sum_m \exp(2\pi i n t_m) f(t_m)
        !! \f]
        !! where t_m must be between -1/2 and 1/2.
        type(nfft_plan), intent(inout) :: ths

        real(c_double), pointer :: x(:)
        complex(c_double_complex), pointer :: f(:), f_hat(:)
        real(c_double) :: freq(ths%nn_total)
        integer :: m, n

        call warn_fallback('nfft_adjoint')
        call c_f_pointer(ths%f_hat, f_hat, (/ths%nn_total/))
        call c_f_pointer(ths%f, f, (/ths%mm_total/))
        call c_f_pointer(ths%x, x, (/ths%mm_total/))
        forall (n=1:ths%nn_total) &
            freq(n) = PI*(2*n - ths%nn_total - 2)

        f_hat = 0
        do n = 1,ths%nn_total
          do m = 1,ths%mm_total
            f_hat(n) = f_hat(n) + &
                    exp(cmplx(0.d0, freq(n)*x(m), kind=c_double_complex))*f(m)
          enddo
        enddo
    end subroutine
end module nfft_fallback

module ft_common
    use MParameters, only: KINDR
    integer :: nfreq = 0, n2freq = 0
    real(KINDR) :: beta = 0.0
contains
    subroutine ft_common_init(new_nfreq, new_n2freq, new_beta)
        integer, intent(in) :: new_nfreq, new_n2freq
        real(KINDR), intent(in) :: new_beta

        if(new_nfreq < 0) stop '[ft_common_init] nfreq must be positive'
        if(mod(new_nfreq,2) /= 0) stop '[ft_common_init] nfreq must be even'
        if(new_n2freq < 0) stop '[ft_common_init] nfreq must be positive'
        if(mod(new_n2freq,2) /= 0) stop '[ft_common_init] nfreq must be even'
        if(new_beta <= 0.) stop '[ft_common_init] beta must be positive'

        nfreq = new_nfreq
        n2freq = new_n2freq
        beta = new_beta
    end subroutine
end module

#ifdef USE_NFFT

module ft_fast
    use MParameters
    use LegendrePoly, only: PI
    use nfft_interface
    use nfft_static
    use iso_c_binding
    use ft_common
contains
    subroutine ft_init()
    end subroutine

    subroutine ft_cleanup()
    end subroutine

    function get_giw(order, m, tau_c, tau_a, m_prefactors) result(giw)
        !< Computes the Green's function in Matsubara frequencies, related to
        !! the Fourier transform of the Hybridisation matrix using NFFT.
        !! \f[
        !!     \beta G(iw) = \beta \int_0^\beta dt \exp(-iwt) G(t)
        !!                 = \sum_{ij} \exp(iw(t_i - t^+_j)) M_ji
        !! \f]
        integer, intent(in) :: order
        real(KINDR), intent(in) :: m(order,order), tau_c(order), tau_a(order)
        real(KINDR), intent(in), optional :: m_prefactors(order)
        complex(KINDC) :: giw(nfreq, 2)

        type(nfft_plan) :: plan
        real(c_double), pointer :: dtau(:,:)
        complex(c_double_complex), pointer :: func(:,:), func_hat(:)

        integer :: i, j
        real(KINDR) :: o2beta, shift
        complex(KINDC) :: exp_c(order), exp_a(order)

        o2beta = 1/(2*beta)

        ! Generate a new NFFT plan
        call nfft_init_1d(plan, 2*nfreq, order*order)

        ! Initialise with scaled tau differences (t^+_j - t_i)/2beta, because
        ! NFFT requires the sampling domain to be (-1/2, 1/2).
        call c_f_pointer(plan%x, dtau, (/order,order/))
        forall(i=1:order, j=1:order) &
            dtau(i,j) = (tau_a(i) - tau_c(j))*o2beta

        ! Do the pre-computation
        if(iand(plan%nfft_flags, PRE_ONE_PSI) /= 0) &
            call nfft_precompute_one_psi(plan)

        ! Absorb the fermionic shift in the frequencies so that they become
        ! bosonic and remember that the formula uses M_ji instead of M_ij. Use
        ! that:
        !
        !  G(iw_n) = sum_ij exp[iw_n(t_i - t+_j)] M_ji
        !          = sum_ij exp[i*pi/beta(2n+k)(t_i - t+_j)] M_ji
        !          = sum_ij exp[2pi*i*(2n)(t_i-t+_j)/(2beta)] *
        !                     exp[-pi*i*k/beta*t+_j] M_ji exp[pi*i*k/beta*t_i]
        !          = NFFT_2n[t_i - t+_j; (exp(t+) * M * exp(t))^T]
        !
        ! where k is the difference between our first frequency (-N/2+1/2) and
        ! the first NFFT frequency (-N/2).
        shift = PI/beta
        exp_c = exp(cmplx(0., -tau_c*shift, kind=KINDC))
        exp_a = exp(cmplx(0., tau_a*shift, kind=KINDC))
        call c_f_pointer(plan%f, func, (/order,order/))
        forall(i=1:order, j=1:order) &
            func(i,j) = exp_c(j)*m(j,i)*exp_a(i)   ! transposition

        ! Do the NFFT
        call nfft_adjoint(plan)
        call c_f_pointer(plan%f_hat, func_hat, (/2*nfreq/))
        giw(:,1) = func_hat(1::2)  ! only every other frequency is meaningful

        ! Now do the same thing for the improved estimators, which are just
        ! a pre-factor to the tau'_j coordinate of the tau Green's function
        if(present(m_prefactors)) then
            forall(j=1:order) &
                func(:,j) = func(:,j) * m_prefactors(j)
            call nfft_adjoint(plan)
            call c_f_pointer(plan%f_hat, func_hat, (/2*nfreq/))
            giw(:,2) = func_hat(1::2)  ! only every other frequency is meaningful
        endif

        call nfft_finalize(plan)
    end function

    function get_giw2(order, m, taus, m_prefactors) result(giw)
        !< Computes the Green's function in Matsubara frequencies, related to
        !! the Fourier transform of the Hybridisation matrix using NFFT.
        !! \f[
        !!     \beta G(iw) = \beta \int_0^\beta dt \exp(-iwt) G(t)
        !!                 = \sum_{ij} \exp(iw(t_i - t^+_j)) M_ji
        !! \f]
        integer, intent(in) :: order
        complex(KINDC), intent(in) :: m(order)
        real(KINDR), intent(in) :: taus(order)
        real(KINDR), intent(in), optional :: m_prefactors(order)
        complex(KINDC) :: giw(nfreq, 2)

        type(nfft_plan) :: plan
        real(c_double), pointer :: dtau(:)
        complex(c_double_complex), pointer :: func(:), func_hat(:)

        integer :: i, j
        real(KINDR) :: o2beta, shift
        complex(KINDC) :: exp_c(order), exp_a(order)

        o2beta = 1/(2*beta)

        ! Generate a new NFFT plan
        call nfft_init_1d(plan, 2*nfreq, order)

        ! Initialise with scaled tau differences (t^+_j - t_i)/2beta, because
        ! NFFT requires the sampling domain to be (-1/2, 1/2).
        call c_f_pointer(plan%x, dtau, (/order/))
        !forall(i=1:order, j=1:order) &
            !dtau(i,j) = (tau_a(i) - tau_c(j))*o2beta
        dtau=taus !*o2beta

        ! Do the pre-computation
        if(iand(plan%nfft_flags, PRE_ONE_PSI) /= 0) &
            call nfft_precompute_one_psi(plan)

        ! Absorb the fermionic shift in the frequencies so that they become
        ! bosonic and remember that the formula uses M_ji instead of M_ij. Use
        ! that:
        !
        !  G(iw_n) = sum_ij exp[iw_n(t_i - t+_j)] M_ji
        !          = sum_ij exp[i*pi/beta(2n+k)(t_i - t+_j)] M_ji
        !          = sum_ij exp[2pi*i*(2n)(t_i-t+_j)/(2beta)] *
        !                     exp[-pi*i*k/beta*t+_j] M_ji exp[pi*i*k/beta*t_i]
        !          = NFFT_2n[t_i - t+_j; (exp(t+) * M * exp(t))^T]
        !
        ! where k is the difference between our first frequency (-N/2+1/2) and
        ! the first NFFT frequency (-N/2).
        !shift = PI/beta
        !exp_c = exp(cmplx(0., -tau_c*shift, kind=KINDC))
        !exp_a = exp(cmplx(0., tau_a*shift, kind=KINDC))
        call c_f_pointer(plan%f, func, (/order/))
        !forall(i=1:order, j=1:order) &
            !func(i,j) = exp_c(j)*m(j,i)*exp_a(i)   ! transposition
        func=m

        ! Do the NFFT
        call nfft_adjoint(plan)
        call c_f_pointer(plan%f_hat, func_hat, (/2*nfreq/))
        giw(:,1) = func_hat(1::2)  ! only every other frequency is meaningful

        !! Now do the same thing for the improved estimators, which are just
        !! a pre-factor to the tau'_j coordinate of the tau Green's function
        !if(present(m_prefactors)) then
            !forall(j=1:order) &
                !func(:,j) = func(:,j) * m_prefactors(j)
            !call nfft_adjoint(plan)
            !call c_f_pointer(plan%f_hat, func_hat, (/2*nfreq/))
            !giw(:,2) = func_hat(1::2)  ! only every other frequency is meaningful
        !endif

        call nfft_finalize(plan)
    end function get_giw2

    function get_g2iw(order, m, tau_c, tau_a) result(g2iw)
        !< Computes the complete two-time Green's function in Matsubara
        !! frequencies (without using time differences), related to the two-
        !! dimensional Fourier transform of the Hybridisation matrix using NFFT:
        !! \f[
        !!     \beta^2 G(iw,iw')
        !!         = \beta^2 \int dt \exp(iwt) \int dt' \exp(-iw't') G(t,t')
        !!         = \sum_{ij} \exp(-iw't'_j) M_ji \exp(+iwt_i)
        !! \f]
        integer, intent(in) :: order
        real(KINDR), intent(in) :: m(order,order), tau_c(order), tau_a(order)
        complex(KINDC) :: g2iw(n2freq,n2freq)

        type(nfft_plan) :: plan
        real(c_double), pointer :: x(:,:,:)
        complex(c_double_complex), pointer :: func(:,:), func_hat(:,:)

        integer :: i, j
        real(KINDR) :: shift
        complex(KINDC) :: exp_c(order), exp_a(order)

        ! Generate a new NFFT plan
        call nfft_init_2d(plan, n2freq, n2freq, order*order)

        ! Initialise with tau coordinates corresponding to the two dimensions
        ! tau/omega (1st in Fortran, 2nd in C) and tau'/omega' (1st in C, 2nd in
        ! Fortran), which are stored as final in C, first in Fortran dimension
        ! in the x vector. Initialise with shifted taus -> tau/beta - 1/2,
        ! because NFFT requires the sampling domain to be (-1/2, 1/2).
        call c_f_pointer(plan%x, x, (/2, order, order/))
        forall(i=1:order, j=1:order)
            x(2,i,j) = tau_a(i)/beta - 0.5
            x(1,i,j) = -tau_c(j)/beta + 0.5
        end forall

        ! Absorb the fermionic shift in the frequencies so that they become
        ! bosonic and remember that the formula uses M_ji instead of M_ij. Use
        ! The two-frequency part is trickier, because of the constant tau shift,
        ! but can be related to the 2D NFFT with the same transformation of M_ji
        !
        !  beta^2 G(iw_m,iw_n) =
        !    = sum_ij exp[-iw_n*t'_j] M_ji exp[iw_m*t_i]
        !    = sum_ij exp[-i*pi/beta*[(2n+k)*t'_j-(2n+k)*t_i] M_ji]
        !    = sum_ij exp[i*pi*2m*t_i/beta] * exp[-i*pi*2n*t'_j/beta] *
        !               exp[-i*pi/beta*k*t'_j]*M_ji*exp[i*pi/beta*k*t_i]
        !    = exp[i*pi*m] * exp[-i*pi*n] *
        !      sum_ij exp[2*i*pi*[m*(t_i/beta-1/2) + n*(-t'_j/beta+1/2)]] *
        !               exp[-i*pi/beta*k*t'_j]*M_ji*exp[i*pi/beta*k*t_i]
        !    = (-1)^(m+n) *
        !      NFFT_{m,n}[(t_i/beta-1/2) 1_ij, (-t'_j/beta+1/2) 1_ij;
        !               exp[-i*pi/beta*k*t'_j]*M_ji*exp[i*pi/beta*k*t_i]_ij ]
        !
        ! where k is the difference between our first frequency (-N/2+1/2) and
        ! the first NFFT frequency (-N/2), 1_ij is the diagonal matrix.
        shift = PI/beta
        exp_c = exp(cmplx(0., -tau_c*shift, kind=KINDC))
        exp_a = exp(cmplx(0., tau_a*shift, kind=KINDC))
        call c_f_pointer(plan%f, func, (/order,order/))
        forall(i=1:order, j=1:order) &
            func(i,j) = exp_c(j)*m(j,i)*exp_a(i)   ! transposition

        ! Do the pre-computation
        if(iand(plan%nfft_flags, PRE_ONE_PSI) /= 0) &
            call nfft_precompute_one_psi(plan)

        ! Do the NFFT
        call nfft_adjoint(plan)
        call c_f_pointer(plan%f_hat, func_hat, (/n2freq,n2freq/))

        ! this is okay because x(2, ...) corresponds to annihilator times, which
        ! the correspond to C dimension 2 in the result, which corresponds to
        ! Fortran dimension 1 because of other pivoting
        g2iw = func_hat

        ! Include the overall sign induced by the shift, the beginning of i,j
        ! being wrong do not matter MOD 2, so the condition is OK.
        forall(i=1:n2freq, j=1:n2freq, mod(i+j,2)==1) &
            g2iw(i,j) = -g2iw(i,j)

        call nfft_finalize(plan)
    end function
end module ft_fast

#endif

module ft_naive
    use MParameters
    use LegendrePoly, only: PI
    use ft_common
contains
    subroutine ft_init()
    end subroutine

    subroutine ft_cleanup()
    end subroutine

    function get_giw(order, m, tau_c, tau_a, m_prefactors) result(giw)
        !< Computes the Green's function in Matsubara frequencies, related to
        !! the Fourier transform of the Hybridisation matrix in a naive way.
        !! \f[
        !!     \beta G(iw) = \beta \int_0^\beta dt \exp(iwt) G(t)
        !!                 = \sum_{ij} \exp(-iwt^+_j) M_ji \exp(+iwt_i)
        !! \f]
        integer, intent(in) :: order
        real(KINDR), intent(in) :: m(order,order), tau_c(order), tau_a(order)
        real(KINDR), intent(in), optional :: m_prefactors(order)
        complex(KINDC) :: giw(nfreq,2)

        real(KINDR) :: freq(nfreq)
        complex(KINDC) :: exp_c(order,nfreq), exp_a(order,nfreq)
        integer :: ifreq, ic, ia

        ! Initialize the Matsubara frequencies
        forall (ifreq=1:nfreq)
            freq(ifreq) = PI/beta*(2*ifreq - nfreq - 1)
        end forall
        forall (ic=1:order,ifreq=1:nfreq)
            exp_c(ic,ifreq) = exp(cmplx(0, -freq(ifreq)*tau_c(ic), KINDC))
            exp_a(ic,ifreq) = exp(cmplx(0, freq(ifreq)*tau_a(ic), KINDC))
        end forall
        ! Do the Fourier transform in a very, very slow way
        giw = 0
        do ifreq = 1, nfreq
          do ic = 1, order
            do ia = 1, order
              giw(ifreq,1) = giw(ifreq,1) + m(ic,ia)*exp_c(ic,ifreq)*exp_a(ia,ifreq)
              ! There must be ic here because the M matrix is transposed before
              ! Fourier transforming
              if(present(m_prefactors)) &
                giw(ifreq,2) = giw(ifreq,2) + m(ic,ia)*m_prefactors(ic)*&
                                    exp_c(ic,ifreq)*exp_a(ia,ifreq)
            enddo
          enddo
        enddo
    end function

    function get_g2iw(order, m, tau_c, tau_a) result(g2iw)
        !< Naively computes the complete two-time Green's function in Matsubara
        !! frequencies (without using time differences), related to the two-
        !! dimensional Fourier transform of the Hybridisation matrix:
        !! \f[
        !!     \beta^2 G(iw,iw')
        !!         = \beta^2 \int dt \exp(iwt) \int dt' \exp(-iw't') G(t,t')
        !!         = \sum_{ij} \exp(-iw't'_j) M_ji \exp(+iwt_i)
        !! \f]
        integer, intent(in) :: order
        real(KINDR), intent(in) :: m(order,order), tau_c(order), tau_a(order)
        complex(KINDC) :: g2iw(n2freq,n2freq)

        real(KINDR) :: freq(n2freq)
        complex(KINDC) :: exp_c(order,n2freq), exp_a(order,n2freq)
        integer :: ifreq, ifc, ifa, ic, ia

        ! Initialize the Matsubara frequencies
        forall (ifreq=1:n2freq)
            freq(ifreq) = PI/beta*(2*ifreq - n2freq - 1)
        end forall
        forall (ic=1:order,ifreq=1:n2freq)
            exp_c(ic,ifreq) = exp(cmplx(0, -freq(ifreq)*tau_c(ic), KINDC))
            exp_a(ic,ifreq) = exp(cmplx(0, freq(ifreq)*tau_a(ic), KINDC))
        end forall
        ! Do the Fourier transform in a very, very slow way. Here we see the
        ! slowness of the two-particle measurement in its full beauty
        g2iw = 0
        do ifc = 1, n2freq
          do ifa = 1, n2freq
            do ic = 1, order
              do ia = 1, order
                g2iw(ifa,ifc) = g2iw(ifa,ifc) + m(ic,ia)*exp_c(ic,ifc)*exp_a(ia,ifa)
              enddo
            enddo
          enddo
        enddo
    end function
end module

#ifdef nfft_Test

program test_nfft
    use nfft_interface
    use nfft_static
    use iso_c_binding
    type(nfft_plan) :: plan
    real(c_double), pointer :: x_ptr(:) => null()
    real(c_double) :: x(100) = (/ &
       -0.49817432, -0.47229947, -0.45946755, -0.44596662, -0.43078845, &
       -0.42326877, -0.41375774, -0.40957637, -0.40340169, -0.4026997 , &
       -0.39271297, -0.37480459, -0.3716082 , -0.36735026, -0.33777499, &
       -0.33383619, -0.29905347, -0.29758886, -0.2952534 , -0.29468725, &
       -0.2875901 , -0.28694288, -0.28096557, -0.27994664, -0.2789838 , &
       -0.26170112, -0.25960972, -0.25513435, -0.25261837, -0.24719763, &
       -0.24177315, -0.21219892, -0.20757987, -0.20067964, -0.20028451, &
       -0.19086766, -0.18504922, -0.18363041, -0.17339269, -0.17334526, &
       -0.15282455, -0.14688955, -0.14660677, -0.13488203, -0.11790803, &
       -0.10591728, -0.10507694, -0.08497412, -0.06373546, -0.03971589, &
       -0.03346923, -0.03313057, -0.02912421, -0.02018274,  0.0036763 , &
        0.03795309,  0.03800827,  0.04349427,  0.0715811 ,  0.08224311, &
        0.08870178,  0.10653369,  0.11300763,  0.11647832,  0.1244734 , &
        0.13097748,  0.15243964,  0.15360636,  0.15416568,  0.16296895, &
        0.17451897,  0.17730087,  0.19297311,  0.19439945,  0.20358798, &
        0.2145724 ,  0.25531845,  0.26285454,  0.28053482,  0.28462854, &
        0.30282576,  0.32007291,  0.3350159 ,  0.34550186,  0.34691832, &
        0.34841742,  0.36486721,  0.3757441 ,  0.37827176,  0.38742559, &
        0.39219971,  0.39488018,  0.40701395,  0.41654002,  0.42279463, &
        0.4668105 ,  0.47637237,  0.48577949,  0.49159105,  0.49232989 /)
    complex(c_double_complex), pointer :: f_ptr(:) => null()
    complex(c_double_complex) :: f(100) = (/ &
        2.84021551e-21,   4.49876634e-21,   2.33634104e-20,   1.96317458e-19, &
        3.88526044e-19,   2.96182070e-17,   4.99780848e-17,   6.44484936e-17, &
        1.05655759e-16,   2.43992770e-16,   3.09471720e-16,   9.96346960e-16, &
        1.03610111e-15,   1.13829987e-14,   1.08795577e-13,   1.68817110e-13, &
        6.25800158e-13,   1.29039322e-12,   1.74540925e-10,   6.45864923e-10, &
        2.48397357e-09,   8.45828695e-09,   1.72803501e-08,   1.67556950e-07, &
        2.93764052e-07,   3.87223354e-07,   6.32013691e-06,   1.20723144e-05, &
        1.57333265e-05,   2.47489217e-05,   2.86780874e-05,   4.31605639e-05, &
        5.49609076e-05,   6.18924423e-05,   1.18504825e-04,   2.10070369e-04, &
        1.87842171e-03,   4.35459002e-03,   9.75486392e-03,   1.51629535e-02, &
        2.44870218e-02,   2.61766208e-02,   2.62022655e-02,   3.25350524e-02, &
        4.71101812e-02,   1.19909968e-01,   1.22264390e-01,   4.40452111e-01, &
        5.00423019e-01,   5.46248171e-01,   5.79070042e-01,   9.99923224e-01, &
        9.92884034e-01,   9.79853839e-01,   9.09223048e-01,   8.92467934e-01, &
        7.61826901e-01,   7.43396665e-01,   7.19297549e-01,   6.65485444e-01, &
        6.03961049e-01,   4.01636239e-01,   3.19600471e-01,   3.18476501e-01, &
        2.37740462e-01,   1.84649287e-01,   1.05632379e-01,   8.11325512e-02, &
        6.49614192e-02,   5.50493924e-02,   3.18150309e-02,   2.30058701e-02, &
        1.33885091e-02,   1.14622953e-02,   3.02138686e-03,   1.45501422e-03, &
        2.70140276e-04,   1.23082501e-05,   6.59682881e-07,   1.15856149e-07, &
        3.05931082e-08,   8.44768072e-09,   2.55808131e-09,   8.16632253e-10, &
        3.52538382e-10,   1.45347541e-10,   2.45182919e-11,   1.22903620e-11, &
        9.40901797e-13,   7.24026177e-13,   1.17777874e-13,   2.11234352e-15, &
        1.89086433e-15,   7.55084748e-17,   7.32853302e-19,   5.85061810e-19, &
        5.16993174e-20,   2.38881158e-20,   1.07077265e-21,   6.30053965e-22 /)
    complex(c_double_complex), pointer :: f_hat_ptr(:) => null()
    complex(c_double_complex) :: f_hat(100) = (/ &
        cmplx( 0.78498172,-0.57681648),  cmplx( 1.18476572,-0.91500728), &
        cmplx( 1.84273712,-1.248964)  ,  cmplx( 2.36715497,-1.21266138), &
        cmplx( 2.28304278,-0.67684259),  cmplx( 1.37721544,+0.16378875), &
        cmplx(-0.11336763,+0.94756317),  cmplx(-1.63384619,+1.40259045), &
        cmplx(-2.65413426,+1.49407277),  cmplx(-2.95091014,+1.37337496), &
        cmplx(-2.64035372,+1.21028388),  cmplx(-1.99568213,+1.06038359), &
        cmplx(-1.23793941,+0.85737798),  cmplx(-0.45960541,+0.50556483), &
        cmplx( 0.31053843,-0.02063512),  cmplx( 1.01307982,-0.64908956), &
        cmplx( 1.54199668,-1.26911584),  cmplx( 1.82090669,-1.82792738), &
        cmplx( 1.89810395,-2.3827698) ,  cmplx( 1.96757954,-3.04783385), &
        cmplx( 2.27296108,-3.85418903),  cmplx( 2.93780901,-4.63338937), &
        cmplx( 3.83598450,-5.0499646) ,  cmplx( 4.61024429,-4.79969554), &
        cmplx( 4.85115347,-3.83890151),  cmplx( 4.32847229,-2.45958509), &
        cmplx( 3.12216156,-1.136306)  ,  cmplx( 1.56864402,-0.25096371), &
        cmplx( 0.06593725,+0.10520242),  cmplx(-1.12734492,+0.12320003), &
        cmplx(-1.97342715,+0.0979178) ,  cmplx(-2.58848058,+0.23229355), &
        cmplx(-3.10253081,+0.55431078),  cmplx(-3.55244289,+0.96205592), &
        cmplx(-3.86602506,+1.32310274),  cmplx(-3.92619212,+1.55361586), &
        cmplx(-3.66298155,+1.64955554),  cmplx(-3.11265672,+1.67485804), &
        cmplx(-2.40402193,+1.72123299),  cmplx(-1.68584677,+1.8615945) , &
        cmplx(-1.06308504,+2.11883688),  cmplx(-0.59602956,+2.45067541), &
        cmplx(-0.32753077,+2.73482107),  cmplx(-0.24254855,+2.76387988), &
        cmplx(-0.13931011,+2.30300445),  cmplx( 0.44548009,+1.24615807), &
        cmplx( 2.06037355,-0.20201536),  cmplx( 4.90458362,-1.49590866), &
        cmplx( 8.43702671,-1.99662164),  cmplx(11.43581322,-1.39580958), &
        cmplx(12.61823645,+0.)        ,  cmplx(11.43581322,+1.39580958), &
        cmplx( 8.43702671,+1.99662164),  cmplx( 4.90458362,+1.49590866), &
        cmplx( 2.06037355,+0.20201536),  cmplx( 0.44548009,-1.24615807), &
        cmplx(-0.13931011,-2.30300445),  cmplx(-0.24254855,-2.76387988), &
        cmplx(-0.32753077,-2.73482107),  cmplx(-0.59602956,-2.45067541), &
        cmplx(-1.06308504,-2.11883688),  cmplx(-1.68584677,-1.8615945) , &
        cmplx(-2.40402193,-1.72123299),  cmplx(-3.11265672,-1.67485804), &
        cmplx(-3.66298155,-1.64955554),  cmplx(-3.92619212,-1.55361586), &
        cmplx(-3.86602506,-1.32310274),  cmplx(-3.55244289,-0.96205592), &
        cmplx(-3.10253081,-0.55431078),  cmplx(-2.58848058,-0.23229355), &
        cmplx(-1.97342715,-0.0979178) ,  cmplx(-1.12734492,-0.12320003), &
        cmplx( 0.06593725,-0.10520242),  cmplx( 1.56864402,+0.25096371), &
        cmplx( 3.12216156,+1.136306)  ,  cmplx( 4.32847229,+2.45958509), &
        cmplx( 4.85115347,+3.83890151),  cmplx( 4.61024429,+4.79969554), &
        cmplx( 3.83598450,+5.0499646) ,  cmplx( 2.93780901,+4.63338937), &
        cmplx( 2.27296108,+3.85418903),  cmplx( 1.96757954,+3.04783385), &
        cmplx( 1.89810395,+2.3827698) ,  cmplx( 1.82090669,+1.82792738), &
        cmplx( 1.54199668,+1.26911584),  cmplx( 1.01307982,+0.64908956), &
        cmplx( 0.31053843,+0.02063512),  cmplx(-0.45960541,-0.50556483), &
        cmplx(-1.23793941,-0.85737798),  cmplx(-1.99568213,-1.06038359), &
        cmplx(-2.64035372,-1.21028388),  cmplx(-2.95091014,-1.37337496), &
        cmplx(-2.65413426,-1.49407277),  cmplx(-1.63384619,-1.40259045), &
        cmplx(-0.11336763,-0.94756317),  cmplx( 1.37721544,-0.16378875), &
        cmplx( 2.28304278,+0.67684259),  cmplx( 2.36715497,+1.21266138), &
        cmplx( 1.84273712,+1.248964)  ,  cmplx( 1.18476572,+0.91500728) /)

    call nfft_init_1d(plan, size(f_hat), size(f))
    call c_f_pointer(plan%x, x_ptr, (/ size(x) /))
    call c_f_pointer(plan%f, f_ptr, (/ size(f) /))
    x_ptr = x
    f_ptr = f
    if(iand(plan%nfft_flags, PRE_ONE_PSI) /= 0) then
        write (*,*) 'I should precompute (this is supposed to happen)'
        call nfft_precompute_one_psi(plan)
    endif
    call nfft_adjoint(plan)
    call c_f_pointer(plan%f_hat, f_hat_ptr, (/ size(f_hat) /))
    write (*,*) x_ptr
    write (*,*) f_hat_ptr
    write (*,*) f_hat
end program

#endif
