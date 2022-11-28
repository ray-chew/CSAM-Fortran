module triangle_mod
    implicit none

    private

    public :: llgrid_t, set_triangle_verts, get_box_width, points_in_triangle

    type :: llgrid_t
        real :: lat_res, lon_res
        real, dimension(3) :: vi, vj
        real :: lat_min, lat_max, lon_min, lon_max
    end type llgrid_t

contains 

    subroutine set_triangle_verts(llgrid, vi, vj)
        real, dimension(3), intent(in):: vi, vj
        type(llgrid_t), intent(inout) :: llgrid

        llgrid%vi = vi
        llgrid%vj = vj
    end subroutine set_triangle_verts


    subroutine get_box_width(lat_vert, lon_vert, llgrid, pad)
        implicit none
        real, dimension(:), intent(in) :: lat_vert
        real, dimension(:), intent(in) :: lon_vert
        type(llgrid_t), intent(out) :: llgrid
        real, optional, intent(in) :: pad
        real :: pad_

        if (present(pad)) then
            pad_ = pad
        else
            pad_ = 4.0
        end if

        llgrid%lat_max = maxval(lat_vert) + pad_
        llgrid%lat_min = minval(lat_vert) - pad_
        llgrid%lon_max = maxval(lon_vert) + pad_
        llgrid%lon_min = minval(lon_vert) - pad_

    end subroutine get_box_width


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

end module triangle_mod