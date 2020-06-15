!> Module providing accumulator types
module accumulators
    use iso_c_binding, only: c_int64_t, c_double
    implicit none

    !> Base accumulator which tracks nothing
    !!
    !!  - Time cost:    O(N)
    !!  - Memory cost:  O(1)
    !!  - Provides:     count
    !!
    type :: accumulator
        real(c_double), allocatable :: buffer(:)
        integer(c_int64_t) :: n
    contains
        procedure :: add => accumulator_add
        procedure :: reset => accumulator_reset
        procedure :: delete => accumulator_delete

        procedure :: num_comp => accumulator_num_comp
        procedure :: num_blocks => accumulator_na

        procedure :: has_mean => accumulator_no
        procedure :: has_blocks => accumulator_no
        procedure :: get_mean => accumulator_handle1d
        procedure :: get_blocks => accumulator_handle2d
    end type

    !> Accumulator which tracks mean
    !!
    !!  - Time cost:    O(k * N)
    !!  - Memory cost:  O(k)
    !!  - Provides:     count, mean
    !!
    type, extends(accumulator) :: meanacc
    contains
        procedure :: add => meanacc_add
        procedure :: reset => meanacc_reset
        procedure :: delete => meanacc_delete

        procedure :: has_mean => meanacc_has_mean
        procedure :: get_mean => meanacc_mean
    end type meanacc

    ! Grouped to improve locality of reference
    type :: varatom
        real(c_double) :: spill
        real(c_double) :: mean
        real(c_double) :: sqdiff
    end type varatom

    !> Accumulator which tracks mean and a stack of variances
    !!
    !!  - Time cost:    O(k * N)
    !!  - Memory cost:  O(k * log N)
    !!  - Provides:     mean, variance, autocorrelation time, spectrum
    !!
    type, extends(accumulator) :: autocorracc
        type(varatom), allocatable :: dat(:,:)
    contains
        procedure :: add => autocorracc_add
        procedure :: reset => autocorracc_reset
        procedure :: delete => autocorracc_delete

        procedure :: has_mean => autocorracc_has_mean
        procedure :: has_blocks => autocorracc_has_blocks
        procedure :: num_blocks => autocorracc_num_blocks
        procedure :: get_mean => autocorracc_mean
        procedure :: get_blocks => autocorracc_blocks
    end type autocorracc

#if defined(__GNUC__) && __GNUC__ > 9
    ! gfortran 10 and up does not allow the specific code used just
    ! below any more (BOZ literal constant cannot appear in an array
    ! constructor), this is essentially the same hack with the integer
    ! literal taken from ieee_quiet_nan on x86_64
    real(c_double), parameter :: nan = transfer(-2251799813685248_c_int64_t, 1.0_c_double)
#else
    ! Define NAN from bitpattern
    ! HACK: we should be using ieee_arithmetic for this, but need to work
    !       around the fact that most clusters last updated their toolchain
    !       when the dinosaurs became extinct.
    real(c_double), parameter :: nan = &
                        transfer((/ Z'00000000', Z'7FF80000' /), 1.0_c_double)
#endif

contains

    pure function ilog2(val) result(res)
        integer(c_int64_t), intent(in) :: val
        integer(c_int64_t) :: res, tmp

        res = -1
        if (val < 1) return

        tmp = val
        do while (tmp /= 0)
            res = res + 1
            tmp = tmp / 2
        enddo
    end function

    ! ######## METHODS FOR BASE ACCUMULATOR ########

    function accumulator_create() result(this)
        type(accumulator) :: this

        this%n = 0
    end function accumulator_create

    subroutine accumulator_reset(this)
        class(accumulator), intent(inout) :: this

        this%n = 0
    end subroutine accumulator_reset

    subroutine accumulator_delete(this)
        class(accumulator), intent(inout) :: this
    end subroutine accumulator_delete

    subroutine accumulator_add(this)
        class(accumulator), intent(inout) :: this

        this%n = this%n + 1
    end subroutine accumulator_add

    function accumulator_no(this) result(answer)
        class(accumulator), intent(in) :: this
        logical :: answer

        answer = .false.
    end function accumulator_no

    function accumulator_num_comp(this) result(answer)
        class(accumulator), intent(in) :: this
        integer(c_int64_t) :: answer

        if (allocated(this%buffer)) then
            answer = size(this%buffer)
        else
            answer = -1
        endif
    end function accumulator_num_comp

    function accumulator_na(this) result(answer)
        class(accumulator), intent(in) :: this
        integer(c_int64_t) :: answer

        answer = -1
    end function accumulator_na

    subroutine accumulator_handle1d(this, out_)
        class(accumulator), intent(in) :: this
        real(c_double), intent(out) :: out_(:)

        out_ = nan
    end subroutine accumulator_handle1d

    subroutine accumulator_handle2d(this, out_)
        class(accumulator), intent(in) :: this
        real(c_double), intent(out) :: out_(:,:)

        out_ = nan
    end subroutine accumulator_handle2d

    ! ######## METHODS FOR MEAN ACCUMULATOR ########

    function meanacc_create(ncomp) result(this)
        type(meanacc) :: this
        integer(c_int64_t), intent(in) :: ncomp

        allocate(this%buffer(ncomp))
        call meanacc_reset(this)
    end function meanacc_create

    subroutine meanacc_reset(this)
        class(meanacc), intent(inout) :: this

        this%n = 0
        this%buffer(:) = 0
    end subroutine meanacc_reset

    function meanacc_has_mean(this) result(answer)
        class(meanacc), intent(in) :: this
        logical :: answer

        answer = .true.
    end function meanacc_has_mean

    subroutine meanacc_delete(this)
        class(meanacc), intent(inout) :: this

        deallocate(this%buffer)
    end subroutine meanacc_delete

    subroutine meanacc_add(this)
        class(meanacc), intent(inout) :: this

        this%n = this%n + 1
    end subroutine meanacc_add

    subroutine meanacc_mean(this, out_)
        class(meanacc), intent(in) :: this
        real(c_double), intent(out) :: out_(:)

        if (size(out_) /= size(this%buffer)) &
            stop '[meanacc_mean] Invalid size of out_ array'

        out_(:) = this%buffer(:) / this%n
    end subroutine meanacc_mean

    ! ######## METHODS FOR AUTOCORRELATION ACCUMULATOR ########

    function autocorracc_create(ncomp, nlevels) result(this)
        type(autocorracc) :: this
        integer(c_int64_t), intent(in) :: ncomp, nlevels

        allocate(this%buffer(ncomp))
        allocate(this%dat(ncomp, nlevels))
        call autocorracc_reset(this)
    end function autocorracc_create

    subroutine autocorracc_delete(this)
        class(autocorracc), intent(inout) :: this

        deallocate(this%buffer)
        deallocate(this%dat)
    end subroutine autocorracc_delete

    subroutine autocorracc_reset(this)
        class(autocorracc), intent(inout) :: this

        this%n = 0
        this%buffer(:) = 0
        this%dat(:,:)%spill = 0
        this%dat(:,:)%mean = 0
        this%dat(:,:)%sqdiff = 0
    end subroutine autocorracc_reset

    function autocorracc_has_mean(this) result(answer)
        class(autocorracc), intent(in) :: this
        logical :: answer

        answer = .true.
    end function autocorracc_has_mean

    function autocorracc_has_blocks(this) result(answer)
        class(autocorracc), intent(in) :: this
        logical :: answer

        answer = .true.
    end function autocorracc_has_blocks

    function autocorracc_num_blocks(this) result(num_blocks)
        class(autocorracc), intent(in) :: this
        integer(c_int64_t) :: num_blocks

        num_blocks = size(this%dat, 2)
    end function autocorracc_num_blocks

    subroutine autocorracc_resize(this, nlevels)
        class(autocorracc), intent(inout) :: this
        integer(c_int64_t), intent(in) :: nlevels

        integer(c_int64_t) :: ncomp, nlevels_old
        type(varatom), allocatable :: tmp(:,:)

        ncomp = size(this%dat, 1)
        nlevels_old = size(this%dat, 2)

        allocate(tmp(ncomp, nlevels))
        if (nlevels < nlevels_old) then
            tmp(:,:) = this%dat(:,:nlevels)
        else
            tmp(:,:nlevels_old) = this%dat
            tmp(:,nlevels_old+1:)%spill = 0
            tmp(:,nlevels_old+1:)%mean = 0
            tmp(:,nlevels_old+1:)%sqdiff = 0
        endif
        call move_alloc(tmp, this%dat)
    end subroutine autocorracc_resize

    subroutine autocorracc_add(this)
        class(autocorracc), intent(inout) :: this

        integer :: lvl
        integer(c_int64_t) :: n, i
        real(c_double) :: val_i, incr1, incr2
        !real(c_double) :: tmp

        n = this%n
        lvl = 1

        ! Handle level 0 first using Welford, West, Hanson update formula
        ! NOTE: we use a do loop here rather than vectorized notation, because
        !       otherwise ifort cannot work out aliasing and starts creating
        !       temporaries.
        do i = 1, size(this%dat, 1)
            val_i = this%buffer(i)
            this%dat(i,lvl)%spill = this%dat(i,lvl)%spill + val_i
            incr1 = val_i - this%dat(i,lvl)%mean
            this%dat(i,lvl)%mean = this%dat(i,lvl)%mean + incr1 / (n + 1)
            incr2 = val_i - this%dat(i,lvl)%mean
            this%dat(i,lvl)%sqdiff = this%dat(i,lvl)%sqdiff + incr1 * incr2
        end do

        ! Whenever the lowest n bits are all ones, we need to roll over results
        ! to the n'th level.  See arxiv:1810.05079, Sec III.D for details
        do while(iand(n, 1_c_int64_t) /= 0)
            if (lvl == size(this%dat, 2)) exit
            n = n / 2
            lvl = lvl + 1

            ! add value of lower-level bin to current level
            do i = 1, size(this%dat, 1)
                ! First, empty lower level
                val_i = this%dat(i,lvl-1)%spill / 2
                this%dat(i,lvl-1)%spill = 0

                ! Add to current level
                this%dat(i,lvl)%spill = this%dat(i,lvl)%spill + val_i
                incr1 = val_i - this%dat(i,lvl)%mean
                this%dat(i,lvl)%mean = this%dat(i,lvl)%mean + incr1 / (n + 1)
                incr2 = val_i - this%dat(i,lvl)%mean
                this%dat(i,lvl)%sqdiff = this%dat(i,lvl)%sqdiff + incr1 * incr2
            end do
        end do

        this%n = this%n + 1
        this%buffer(:) = 0
    end subroutine autocorracc_add

    subroutine autocorracc_mean(this, out_)
        class(autocorracc), intent(in) :: this
        real(c_double), intent(out) :: out_(:)

        if (size(out_) /= size(this%dat, 1)) &
            stop '[autocorracc_mean] Invalid size of out_ array'

        out_(:) = this%dat(:,1)%mean
    end subroutine autocorracc_mean

    subroutine autocorracc_blocks(this, out_)
        class(autocorracc), intent(in) :: this
        real(c_double), intent(out) :: out_(:,:)

        integer(c_int64_t) :: n
        integer :: level

        if (size(out_,1) /= size(this%dat,1)) &
            stop '[autocorracc_mean] Invalid size of out_ array'

        n = this%n
        do level = 1, min(size(this%dat, 2), size(out_, 2))
            if (n > 1) then
                out_(:,level) = this%dat(:,level)%sqdiff / (n - 1)
            else
                out_(:,level) = this%dat(:,level)%sqdiff / 0.
            endif
            n = n / 2
        end do
    end subroutine
end module accumulators

#ifdef ACC_TEST_PROGRAM

program test_program
    use accumulators
    use iso_c_binding
    implicit none

    class(accumulator), allocatable, target :: any_acc

    integer(c_int64_t), parameter :: n = 23, ncomp = 1
    logical :: use_autocorr = .true.

    real(c_double), pointer :: ptr(:,:)

    integer :: i
    real(c_double) :: vars(ncomp,n), mean(ncomp)

    if (use_autocorr) then
        allocate(any_acc, source=autocorracc_create(ncomp, n))
    else
        allocate(any_acc, source=meanacc_create(ncomp))
    endif

    ptr(1:1,1:1) => any_acc%buffer

    write (*,*) 'ADDING', 2**n, 'MEASUREMENTS'
    do i = 0, 2**n - 1
        any_acc%buffer(1) = i
        !ptr(1,1) = ptr(1,1) + i
        call any_acc%add()
        !ptr(1,1) = 0
    end do

    if (any_acc%has_mean()) then
        write (*,*) 'MEAN:'
        call any_acc%get_mean(mean)
        write (*,'(ES20.10)') mean(1)
    endif
    if (any_acc%has_blocks()) then
        write (*,*) 'VARIANCES IN BLOCKING ANALYSIS:'
        call any_acc%get_blocks(vars)
        do i = 1, n
            write (*,'(I3, I10, ES20.10)') i-1, 2**(i-1), vars(1,i)
        enddo
    endif

    call any_acc%delete()
    deallocate(any_acc)
end program

#endif
