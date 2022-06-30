module fourier_trans_mod
    implicit none

    private
    integer, parameter ::   nhar_i = 15, &
                            nhar_j = 15

    public :: llgrid_t, set_triangle_verts

    type :: llgrid_t
        real :: lat_res, lon_res
        real, dimension(3) :: vi, vj
        integer, dimension(:), allocatable :: Ni, Nj

    end type llgrid_t

contains

    subroutine set_triangle_verts(llgrid, vi, vj)
        real, dimension(3), intent(in):: vi, vj
        type(llgrid_t), intent(inout) :: llgrid

        llgrid%vi = vi
        llgrid%vj = vj
    end subroutine set_triangle_verts

end module fourier_trans_mod