module utils_mod
    use, intrinsic :: iso_fortran_env, only : error_unit
    use :: error_status, only : COMMAND_LINE_ERR
    use :: fourier_trans_mod, only : llgrid_t
    implicit none

    private
    real, parameter :: PI = acos(-1.0)

    public :: get_fn, rad_to_deg, points_in_triangle

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

    ! translated from Niraj's code
    elemental function points_in_triangle(i, j, ll_grid) result(mask)
        implicit none
        real, intent(in) :: i, j
        type(llgrid_t), intent(in) :: ll_grid
        real, dimension(3) :: tmp_vi, tmp_vj, ei, ej, p2ei, p2ej
        real :: r1, r2, r3
        logical :: mask

        associate (vi => ll_grid%vi, vj => ll_grid%vj)

        tmp_vi = (/vi(2), vi(3), vi(1)/)
        tmp_vj = (/vj(2), vj(3), vj(1)/)
        ei = tmp_vi - vi
        ej = tmp_vj - vj

        p2ei = vi - i
        p2ej = vj - j

        end associate

        r1 = cross_2D(ei(1), ej(1), p2ei(1), p2ej(1))
        r2 = cross_2D(ei(2), ej(2), p2ei(2), p2ej(2))
        r3 = cross_2D(ei(3), ej(3), p2ei(3), p2ej(3))
        
        if (r1 .gt. 0 .and. r2.gt. 0 .and. r3 .gt. 0) then
            mask = .true.
        else if (r1 .lt. 0 .and. r2 .lt. 0 .and. r3 .lt. 0) then
            mask = .true.
        else
            mask = .false.
        end if
    end function points_in_triangle


    pure function cross_2D(x1,x2,y1,y2) result(res)
        implicit none
        real, intent(in) :: x1, x2, y1, y2
        real :: res

        res = x1 * y2 - (y1 * x2)

    end function cross_2D

end module utils_mod