!> Module providing a basic progress indicator
module type_progress
    use iso_c_binding, only: c_int64_t
    integer, parameter :: PINT = c_int64_t

    integer, parameter :: TITLELEN = 15

    !> The basic type for the progress indicator
    type progress
        character(len=TITLELEN) :: title = 'Processing... '
        integer(PINT) :: minskip   = 1        !< Minimum count skips
        integer(PINT) :: maxskip   = HUGE(1_PINT)  !< Maximum count skips
        integer(PINT) :: unit      = 0        !< Output unit
        real          :: updateint = 60.      !< Desired update interval [s]
        real          :: adaptrate = 0.5      !< How fast the system adapts (0, 1)
        logical(PINT) :: fancy     = .false.  !< Interactive display (walltime)

        integer(PINT) :: ticks     = 0        !< Current value
        integer(PINT) :: target    = -1       !< Target count
        integer(PINT) :: uptick    = 1        !< Next checking value
        integer(PINT) :: skip      = 1        !< Current desired skip rate
        real          :: startrun  = 0        !< Progress run start [s]
        real          :: startintv = 0        !< Current interval start [s]
        integer(PINT) :: updates   = 0        !< Also direction of the cursor
    end type

    !> Parameters for "ramping up" the update interval in cases where the
    !! average tick rate is small and possibly of the order of magnitude
    !! of the cpu_time() discretisation.
    real, parameter :: RAMP_FACTOR = 1.7, RAMP_THRESHOLD = .1

contains
    !> Initialises and starts a new progress run
    subroutine pstart(p, total, fancy, title)
        type(progress), intent(inout) :: p
        integer(PINT), intent(in) :: total
        logical, intent(in), optional :: fancy
        character(len=*), intent(in), optional :: title

        if(present(fancy)) then
            p%fancy = fancy
            if(fancy) p%updateint = 0.3
        endif
        if(present(title)) then
            p%title = repeat(' ',TITLELEN)
            write(p%title,'(A)') title
        endif
        if(total < 0) then
            write(p%unit,"(A,' WARNING: negative target: ',I15)") &
                p%title, total
        endif
        p%target = total
        p%ticks = 0
        p%updates = 0
        p%skip = 1
        p%uptick = 1
        if(p%fancy) then
            write(p%unit,"(A,'   0% [',40X,'] START (no tick)')", advance="no") p%title
        else
            write(p%unit,"(A,' starting:',I12,' ticks')") p%title, p%target
        endif
        call cpu_time(p%startrun)
        p%startintv = p%startrun
    end subroutine

    !> Increments the counter by one
    subroutine ptick(p)
        type(progress), intent(inout) :: p
        character(len=1), parameter :: backsl = char(92)
        character, parameter :: DIRCHARS(4) = (/ '|', '/', '-', backsl /)
        real :: t, dt, trem
        integer :: eqs
        integer(PINT), parameter :: modarg = 4

        ! Increment the ticks and check for the next update
        p%ticks = p%ticks + 1
        if(p%ticks < p%uptick) return
        ! Check for a completed cycle
        if(p%ticks >= p%target) then
            call cpu_time(p%startintv)
            write(p%unit,"(A,A,' complete:',I12,' ticks in ',A6,28(' '))")&
                achar(13), p%title, p%ticks, tstr(p%startintv-p%startrun)
            p%target = huge(1_PINT)    ! Flooding prevention
            return
        endif

        ! Get CPU time and check for zero difference (since, e.g., gfortran has
        ! a 10ms discretisation in the returned seconds)
        call cpu_time(t)
        dt = t - p%startintv
        if(dt < RAMP_THRESHOLD) then
            p%uptick = min(p%uptick + p%skip, p%target)
            p%skip = ceiling(RAMP_FACTOR * p%skip)
            return
        endif
        t = t - p%startrun
        p%updates = p%updates + 1

        ! Dynamically adapt the skipping rate
        p%skip = p%skip + int(p%skip*(p%updateint-dt)/dt*p%adaptrate)
        if(p%skip > p%maxskip) p%skip = p%maxskip
        if(p%skip < p%minskip) p%skip = p%minskip
        p%uptick = min(p%ticks + p%skip, p%target)

        ! Print
        trem = (t/p%ticks)*(p%target-p%ticks)
        if(p%fancy) then
            eqs = floor(p%ticks*40.0/p%target)   ! otherwise risk of overflow
            write(p%unit,"(A,A,1X,I3,'% [',A40,'] ',A6,' (T-',A6,')')", advance="no")&
                achar(13), p%title, floor(p%ticks*100.0/p%target), &
                repeat('=',eqs) // DIRCHARS(mod(p%updates,modarg)+1) // repeat(' ',39-eqs),&
                tstr(t), tstr(trem)
        else
            write(p%unit,"(A,I12,' of',I12,' (',F5.1,'%) ',A6,'  T-',A6)") &
                p%title, p%ticks, p%target, p%ticks*100.0/p%target, tstr(t), tstr(trem)
        endif
        call cpu_time(p%startintv)
    contains
        function tstr(seconds) result(res)
            character(len=6) :: res
            real, intent(in) :: seconds
            integer(PINT) :: secs, mins, hrs, days
            if(seconds < 0) then
                res = ' *****'
                return
            elseif(seconds < 60) then
                write(res,"(F5.2,'s')") seconds
                return
            endif
            secs = floor(seconds)
            if(secs < 0) secs = HUGE(1) ! floor overflow
            mins = secs/60
            secs = secs - mins*60
            if(mins < 60) then
                write(res,"(I2,'m',I2.2,'s')") mins, secs
                return
            endif
            hrs = mins/60
            mins = mins - hrs*60
            if(hrs < 48) then
                write(res,"(I2,'h',I2.2,'m')") hrs, mins
                return
            endif
            days = hrs/24
            hrs = hrs - days*24
            if(days < 100) then
                write(res,"(I2,'d',I2,'h')") days, hrs
            elseif(hrs < 100000) then
                write(res,"(I5,'d')") days
            else
                res = '******'
            endif
        end function
    end subroutine

end module

#ifdef Progress_Test

program Test_Progress
    use type_progress
    integer(PINT) :: i, k
    real :: q
    type(progress) :: p
    interface
        subroutine usleep(useconds) bind(C)
            use iso_c_binding
            integer(c_int32_t), value :: useconds
        end subroutine
    end interface

    call pstart(p,10000000_PINT,.true.)
    do i = 1, 10000000_PINT
        call random_number(q)
        do k = 1, int(1000*q); enddo
        call ptick(p)
    end do
end program Test_Progress

#endif
