! ISO C binding wrapper around LIDORT_getPlanck

module planck_m

    use iso_c_binding
    use LIDORT_getPlanck, only :get_planckfunction, get_planckfunction_plus  

contains

    subroutine planckfunction( &
        wavenumber, temperature,   &
        bbfunc, smallv, fail_out, message_len, message_out) bind(c)

        ! Arguments
        real(kind=c_double), intent(in)       :: wavenumber
        real(kind=c_double), intent(in)       :: temperature

        real(kind=c_double), intent(out)      :: bbfunc

        integer(kind=c_int),  intent(out)     :: smallv
        logical(kind=c_bool), intent(out)     :: fail_out

        ! A little bit longer than the maximum message emitted by lidort_getplanck.f90
        integer(kind=c_int),  intent(in)     :: message_len
        character(kind=c_char), intent(inout) :: message_out(message_len)
  
        ! Local variables
        logical :: fail_lcl
        character(kind=c_char, len=message_len) :: message_lcl

        call get_planckfunction(wavenumber - 0.5d0, wavenumber + 0.5d0, temperature, &
            bbfunc, smallv, fail_lcl, message_lcl)

        ! Copy logical value manually
        fail_out = fail_lcl

        ! Copy contents of message and add a null character at the end for C
        message_out = trim(message_lcl) // c_null_char

    end subroutine

    subroutine planckfunction_plus( &
        wavenumber, temperature,   &
        bbfunc, deriv_bbfunc, smallv, fail_out, message_len, message_out) bind(c)

        ! Arguments
        real(kind=c_double), intent(in)       :: wavenumber
        real(kind=c_double), intent(in)       :: temperature

        real(kind=c_double), intent(out)      :: bbfunc
        real(kind=c_double), intent(out)      :: deriv_bbfunc

        integer(kind=c_int),  intent(out)     :: smallv
        logical(kind=c_bool), intent(out)     :: fail_out

        ! A little bit longer than the maximum message emitted by lidort_getplanck.f90
        integer(kind=c_int),  intent(in)     :: message_len
        character(kind=c_char), intent(inout) :: message_out(message_len)
  
        ! Local variables
        logical :: fail_lcl
        character(kind=c_char, len=message_len) :: message_lcl

        call get_planckfunction_plus(wavenumber - 0.5d0, wavenumber + 0.5d0, temperature, &
            bbfunc, deriv_bbfunc, smallv, fail_lcl, message_lcl)

        ! Copy logical value manually
        fail_out = fail_lcl

        ! Copy contents of message and add a null character at the end for C
        message_out = trim(message_lcl) // c_null_char

    end subroutine

end module planck_m
