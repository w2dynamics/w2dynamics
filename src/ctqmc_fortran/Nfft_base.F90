module nfft_static
    !< Provides static types for NFFT
    use iso_c_binding

    type, bind(C) :: nfft_plan
        ! c_intptr_t should ideally be c_ptrdiff_t as NFFT_INT is
        ! currently defined in nfft3.h as ptrdiff_t, but this is not
        ! supported by Intel yet
        integer(c_intptr_t) :: nn_total, mm_total
        type(c_ptr)    :: f_hat, f
        !FIXME: these dont show up in the newest documentation, but in
        !the
        !nfft3.h -> we should really make this version independent
        type(c_funptr) :: mv_trafo, mv_adjoint
        integer(c_intptr_t) :: d
        type(c_ptr)    :: nn,sigma,n
        integer(c_intptr_t) :: n_total,m
        type(c_ptr)    :: b
        integer(c_intptr_t) :: kk
        integer(c_int) :: nfft_flags, fftw_flags
        type(c_ptr)    :: x
        real(c_double) :: measure_time_t(3)
        type(c_ptr)    :: my_fftw_plan1, my_fftw_plan2
        type(c_ptr)    :: c_phi_inv, psi, psi_index_g, psi_index_f
        type(c_ptr)    :: g, g_hat, g1, g2, spline_coeffs, index_x
    end type

    !#define PRE_PHI_HUT   (1U<< 0) 
    !#define FG_PSI   (1U<< 1) 
    !#define PRE_LIN_PSI   (1U<< 2)
    !#define PRE_FG_PSI   (1U<< 3)
    !#define PRE_PSI   (1U<< 4) 
    !#define PRE_FULL_PSI   (1U<< 5) 
    !#define MALLOC_X   (1U<< 6) 
    !#define MALLOC_F_HAT   (1U<< 7) 
    !#define MALLOC_F   (1U<< 8) 
    !#define FFT_OUT_OF_PLACE   (1U<< 9) 
    !#define FFTW_INIT   (1U<< 10) 
    !#define FFTW_DESTROY_INPUT (1U << 0)
    !#define FFTW_MEASURE (0U)
    !#define FFTW_ESTIMATE (1U << 6)
    !#define PRE_ONE_PSI (PRE_LIN_PSI| PRE_FG_PSI| PRE_PSI|
    !PRE_FULL_PSI) 

    integer(c_int), parameter :: PRE_PHI_HUT = ISHFT(1,0)
    integer(c_int), parameter :: FG_PSI = ISHFT(1,1)
    integer(c_int), parameter :: PRE_LIN_PSI = ISHFT(1,2)
    integer(c_int), parameter :: PRE_FG_PSI = ISHFT(1,3)
    integer(c_int), parameter :: PRE_PSI = ISHFT(1,4)
    integer(c_int), parameter :: PRE_FULL_PSI = ISHFT(1,5)
    integer(c_int), parameter :: MALLOC_X = ISHFT(1,6)
    integer(c_int), parameter :: MALLOC_F_HAT = ISHFT(1,7)
    integer(c_int), parameter :: MALLOC_F = ISHFT(1,8)

    integer(c_int), parameter :: FFT_OUT_OF_PLACE = ISHFT(1,9)
    integer(c_int), parameter :: FFTW_INIT = ISHFT(1,10)
    integer(c_int), parameter :: FFTW_DESTROY_INPUT = ISHFT(1,0)
    integer(c_int), parameter :: FFTW_ESTIMATE = ISHFT(1,6)
    integer(c_int), parameter :: FFTW_MEASURE = 0
    integer(c_int), parameter :: PRE_ONE_PSI = 60
    !integer(c_int), parameter :: PRE_ONE_PSI =
    !OR(PRE_LIN_PSI,OR(PRE_FG_PSI,OR(PRE_PSI,PRE_FULL_PSI))

end module nfft_static

module nfft_interface
    !< Interfaces to the C NFFT routines
    use iso_c_binding
    implicit none
    interface
        subroutine nfft_init_1d(ths, n1, m) bind(C)
            use nfft_static
            type(nfft_plan), intent(inout) :: ths
            integer(c_int), value :: n1, m
        end subroutine 
        
        subroutine nfft_init_2d(ths, n1, n2, m) bind(C)
            use nfft_static
            type(nfft_plan), intent(inout) :: ths
            integer(c_int), value :: n1, n2, m
        end subroutine 
        
        subroutine nfft_init_3d(ths, n1, n2, n3, m) bind(C)
            use nfft_static
            type(nfft_plan), intent(inout) :: ths
            integer(c_int), value :: n1, n2, n3, m
        end subroutine
    
        subroutine nfft_precompute_one_psi(ths) bind(C)
            use nfft_static
            type(nfft_plan), intent(inout) :: ths
        end subroutine
    
        subroutine nfft_adjoint(ths) bind(C)
            use nfft_static
            type(nfft_plan), intent(inout) :: ths
        end subroutine
        
        subroutine nfft_check(ths) bind(C)
            use nfft_static
            type(nfft_plan), intent(inout) :: ths
        end subroutine
        
        subroutine nfft_finalize(ths) bind(C)
            use nfft_static
            type(nfft_plan), intent(inout) :: ths
        end subroutine
    
        subroutine nfft_init_guru(ths, d, nn, m_total, n, m, nfft_flags,fftw_flags) bind(C)
            use nfft_static
            type(nfft_plan), intent(inout) :: ths
            integer(c_int),value           :: d
            type(c_ptr),value              :: nn
            integer(c_int),value           :: m_total 
            type(c_ptr),value              :: n
            integer(c_int),value           :: m 
            integer(c_int),value           :: nfft_flags, fftw_flags
        end subroutine        
   end interface
end module nfft_interface

