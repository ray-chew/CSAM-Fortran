module utils_mod
    use, intrinsic :: iso_fortran_env, only : error_unit
    use error_status, only : COMMAND_LINE_ERR
    implicit none

    private
    real, parameter :: PI = acos(-1.0)

    public :: get_fn, rad_to_deg

contains

    subroutine get_fn(fn_grid, fn_topo)
        implicit none
        character(len=*), intent(out) :: fn_grid
        character(len=*), intent(out) :: fn_topo

        if (command_argument_count() /= 2) then
            write(unit=error_unit, fmt='(A)') "Argument error: Expected 2 input arguments."
            stop COMMAND_LINE_ERR
        end if

        call get_command_argument(1, fn_grid)
        call get_command_argument(2, fn_topo)

    end subroutine get_fn

    elemental subroutine rad_to_deg(value)
        implicit none
        real, intent(inout) :: value

        value = value * (180.0 / PI)
        
    end subroutine rad_to_deg


end module utils_mod