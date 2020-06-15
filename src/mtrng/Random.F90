module MRandom
    use ISO_C_BINDING

    ! Binds into random.cc PRNG
    interface
        function new_random() bind(C)
            use ISO_C_BINDING
            type(c_ptr) :: new_random
        end function

        subroutine delete_random(rng) bind(C)
            use ISO_C_BINDING
            type(c_ptr), intent(in), value :: rng
        end subroutine

        subroutine seed_random(rng, seed) bind(C)
            use ISO_C_BINDING
            type(c_ptr), intent(in), value :: rng
            integer(c_int), intent(in), value :: seed
        end subroutine

        function get_random(rng) bind(C)
            use ISO_C_BINDING
            real(c_double) :: get_random
            type(c_ptr), intent(in), value :: rng
        end function

        function get_randint(rng, min, max) bind(C)
            use ISO_C_BINDING
            integer(c_int) :: get_randint
            type(c_ptr), intent(in), value :: rng
            integer(c_int), intent(in), value :: min, max
        end function
    end interface

    type(c_ptr) :: singleton = C_NULL_PTR

contains
    ! Initialise the PRNG
    subroutine sgrnd(seed)
        integer, intent(in) :: seed

        if (C_ASSOCIATED(singleton)) call delete_random(singleton)
        singleton = new_random()
        call seed_random(singleton, seed)
    end subroutine

    ! Get random number in the half-open interval [0,1)
    real(kind=kind(0.D0)) function grnd()
        grnd = get_random(singleton)
    end function

    ! Get random number in the half-open interval [0,1)
    integer(c_int) function randint(min, max)
        integer(c_int), intent(in) :: min, max
        randint = get_randint(singleton, min, max)
    end function
end module
