module utils_mod
    use, intrinsic :: iso_fortran_env, only : error_unit
    use :: error_status, only : COMMAND_LINE_ERR
    use :: const_mod, only : PI
    ! use :: stdlib_sorting, only : ord_sort
    implicit none

    private

    public :: get_fn, rad_to_deg, get_N_unique

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


    function get_N_unique(arr) result(cnt)
        implicit none
        integer, dimension(:), intent(in) :: arr
        integer :: old_val = 0, cnt, i = 1
        ! integer, dimension(:), allocatable :: work

        ! the first value in the array is always unique.
        old_val = arr(1)
        cnt = 1

        do i=2,size(arr)
            if (arr(i) /= old_val) then
                cnt = cnt + 1
                old_val = arr(i)
            end if
        end do
    end function get_N_unique

end module utils_mod