module signals
    !< Allows the more or less portable handling of UNIX signals in Fortran 90.
    !! Provides constants and functions.
    use iso_c_binding
#ifdef __INTEL_COMPILER
    use iflport
#endif
    integer, parameter :: SIGNAL_INT  =  2    !< Interrupt program signal (^C)
    integer, parameter :: SIGNAL_USR1 = 10    !< User-defined signal 1
    integer, parameter :: SIGNAL_USR2 = 12    !< User-defined signal 2
    integer, parameter :: SIGNAL_TERM = 15    !< Terminate signal (kill)
    integer, parameter :: SIGNAL_CONT = 18    !< Continue program (fg)
    integer, parameter :: SIGNAL_STOP = 19    !< Stop program (^Z)
    integer, parameter :: SIGNAL_MAX  = 64    !< Highest signal number

    integer, parameter :: HANDLER_DFL = 0     !< Default handler for signals
    integer, parameter :: HANDLER_ERR = -1    !< Error return code

    !> \{ Variable set by trap_notify() to notify the program of any caught
    !!    signals for which it was registered. Note that these values are
    !!    changed asynchronously, so be especially careful with last_signal.
    logical, volatile :: any_signal_fired = .false. !< Any signal caught
    logical, volatile :: signals_fired(SIGNAL_MAX) = .false. !< Signals caught
    integer, volatile :: last_signal = -1           !< Last signal caught
    logical, volatile :: write_stderr = .true.      !< Write message (unsafe)
    !! \}

    logical :: repeat_aborts = .true.

    logical, parameter :: sigdbg = .false.
contains
    subroutine register_trap(sig_num, trap, prev_trap)
        !< \brief Registers a signal handler for a certain signal.
        !!
        !! There are two predefined signal handlers: trap_stop(), which stops
        !! the program and trap_notify(), which sets the signal_fired modules
        !! variables.
        !!
        !! Note that the usual restrictions for signal handlers apply (you should
        !! not do any I/O or other asynchronous calls). Also note that due to a
        !! GCC bug, you should never use the sig argument of your handler directly,
        !! but pass it first through signal_number().
        integer, intent(in) :: sig_num
        interface
            integer function trap(sig)
                integer, intent(in) :: sig
            end function
        end interface
        integer, intent(out), optional :: prev_trap
#ifdef __INTEL_COMPILER
        integer :: prev
        if(sigdbg) write(*,*) 'Registering handler (Intel)'
        prev = signal(sig_num, trap, -1)
        if(prev == HANDLER_ERR) stop '[register_trap] Error registering handler'
        if(present(prev_trap)) prev_trap = prev
#else
#ifdef __xlc__
        ! untested path for IBM compiler
        if(sigdbg) write(*,*) 'Registering handler (IBM)'
        call signal(sig_num, trap)
#else
        ! GNU and derivatives (hopefully)
        if(sigdbg) write(*,*) 'Registering handler (GNU)'
        call signal(sig_num, trap, prev_trap)
#endif
#endif
    end subroutine

    subroutine clear_trap(sig_num, prev_trap)
        !< Resets the signal handler for a certain signal to the default action.
        integer, intent(in) :: sig_num
        integer, intent(out), optional :: prev_trap
#if defined(__INTEL_COMPILER)
        integer :: prev
        if(sigdbg) write(*,*) 'Clearing handler (Intel)'
        prev = signal(sig_num, trap_stop, 0)
        if(prev == HANDLER_ERR) stop '[clear_trap] Error clearing handler'
        if(present(prev_trap)) prev_trap = prev
#else
#ifdef __xlc__
        ! untested path for IBM compiler. There seems to be no facility to clear a signal handler.
        if(sigdbg) write(*,*) 'Clearing handler (IBM) not implemented!'
!        call signal(sig_num, SIG_DFL)
#else
        ! GNU and derivatives (hopefully)
        if(sigdbg) write(*,*) 'Clearing handler (GNU)'
        call signal(sig_num, HANDLER_DFL, prev_trap)
#endif
#endif
    end subroutine

    recursive integer function signal_number(sig)
        !< \brief Extract "true" signal number out of the handler´s sig parameter
        !!
        !! This function is necessary because there is a bug in gfortran until
        !! v4.7 which causes it to pass the signal number by value instead of by
        !! pointer, which would be Fortran standard.
        integer, intent(in) :: sig
        integer(c_intptr_t) :: addr

        addr = loc(sig)
        if (addr < 0 .or. addr > SIGNAL_MAX) then
            signal_number = sig
        else
            signal_number = int(addr)
        endif
    end function

    recursive integer function trap_stop(sig)
        !< Predefined signal handler, which causes the program to halt
        integer, intent(in) :: sig
        write (0,'(A)') '[trap_stop] Stopping program on signal'
        call abort()
        trap_stop = sig
    end function

    recursive integer function trap_notify(sig)
        !< \brief Predefined signal handler, which records caught signals
        !!
        !! This handler updates the module variables any_signal_fired,
        !! signals_fired and last_signal on every call based on the signal number
        !! caught. Use reset_notify() to clear the statistics.
        !!
        !! Note that this function changes these values asynchronously, so be
        !! especially careful with last_signal.
        integer, intent(in) :: sig
        any_signal_fired = .true.
        last_signal = signal_number(sig)
        if(sigdbg .or. write_stderr) &
            write(0,'(A,I3)') '[trap_notify] Trapped signal', last_signal
        if(last_signal > 0 .and. last_signal <= SIGNAL_MAX) then
            if(repeat_aborts .and. signals_fired(last_signal)) &
                stop 'Aborting: Caught signal second time.'
            signals_fired(last_signal) = .true.
        endif
        trap_notify = 0
    end function

    subroutine reset_notify
        !< \brief Clears all statistics maintained by trap_notify().
        any_signal_fired = .false.
        signals_fired = .false.
        last_signal = -1
    end subroutine
end module

#ifdef Signals_Test

program signals_test
    use signals
    integer :: i
    call register_trap(SIGNAL_INT, trap_notify)
    call register_trap(SIGNAL_USR1, trap_notify)
    do i = 1,100
        write(*,'(A)',advance='no') '.'
        call sleep(1)
        if(any_signal_fired) then
            write (*,*) 'TRAPPED:', last_signal
            call clear_trap(last_signal)
            call reset_notify()
        endif
    enddo
end program signals_test

#endif
