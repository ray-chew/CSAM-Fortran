module fourier_mod
    use :: topo_mod, only : topo_t
    use :: utils_mod, only : get_N_unique
    use :: const_mod, only : PI
    use :: stdlib_sorting, only : ord_sort
    implicit none

    private

    public :: get_coeffs, recover_coeffs, nhar_i, nhar_j

    integer, parameter ::   nhar_i = 12, &
                            nhar_j = 12

contains

    subroutine get_coeffs(topo_obj, mask, coeffs)
        implicit none
        type(topo_t), intent(inout) :: topo_obj
        logical, dimension(:,:), intent(in) :: mask
        real, dimension(count(mask)) :: lat_tri, lon_tri, topo_tri
        integer, dimension(count(mask)) :: II, JJ, JJ_tmp
        real :: d_lat, d_lon
        integer :: Ni, Nj, i, j, k, l, m, n, N_cos, N_sin

        real, dimension(:), allocatable :: tmp, c_cos, c_sin
        real, dimension(:,:), allocatable, intent(out) :: coeffs

        ! we get lat, lon and topo in the triangle using the mask
        lat_tri = pack(topo_obj%lat_grid, mask=mask)
        lon_tri = pack(topo_obj%lon_grid, mask=mask)
        topo_tri = pack(topo_obj%topo, mask=mask)

        ! get grid spacing assuming equidistant grid.
        d_lat = topo_obj%lat(2) - topo_obj%lat(1)
        d_lon = topo_obj%lon(2) - topo_obj%lon(1)

        II = ceiling((lat_tri - minval(lat_tri)) / d_lat)
        JJ = ceiling((lon_tri - minval(lon_tri)) / d_lon)
        
        Ni = get_N_unique(II)
        JJ_tmp = JJ
        call ord_sort(JJ_tmp)
        Nj = get_N_unique(JJ_tmp)

        ! N_cos = nhar_i * nhar_j
        ! N_sin = nhar_i * nhar_j - 1
        N_cos = nhar_i * nhar_j - nhar_j/2
        N_sin = nhar_i * nhar_j - nhar_j/2 - 1

        allocate (tmp(nhar_i * nhar_j))
        allocate (coeffs(N_cos + N_sin, size(topo_tri)))

        allocate (c_cos(N_cos))
        allocate (c_sin(N_sin))

        ! print *, "Entering loop..."

        do k=1,size(topo_tri)
            l = 1
            m = 1
            n = 1
            do i=0,nhar_i-1
                ! do j=0,nhar_j-1
                do j=-nhar_j/2,nhar_j/2-1
                    tmp(l) = 2.0 * PI * (i*II(k)/real(Ni) + j*JJ(k)/real(Nj))

                    if (.not.(i == 0 .and. j < 0)) then
                        c_cos(n) = cos(tmp(l))
                        n = n + 1
                    end if

                    if (.not.(i == 0 .and. j <= 0)) then
                        c_sin(m) = sin(tmp(l))
                        m = m + 1
                    end if


                    ! if (cos(tmp(l)) == 0) then
                    !     print *, l, k, i, j
                    ! end if
                    ! if (sin(tmp(l)) == 0) then
                    !     print *, l, k, i, j
                    ! end if

                    l = l + 1
                end do
            end do

            coeffs(1:N_cos,k) = c_cos
            coeffs(N_cos+1:N_cos+N_sin,k) = c_sin

            ! do i=1,size(coeffs,dim=1)
            !     if (coeffs(i,k) == 0) then
            !         print *, i, k
            !     end if
            ! end do
        end do

        topo_obj%lat_tri = lat_tri
        topo_obj%lon_tri = lon_tri
        topo_obj%topo_tri = topo_tri

    end subroutine get_coeffs


    subroutine recover_coeffs(topo_obj, sol)
        implicit none
        type(topo_t), intent(inout) :: topo_obj
        real, dimension(:), intent(in) :: sol

        complex, allocatable :: recov_coeffs(:), recov_coeffs_ij(:,:)
        integer :: mid = nhar_i * nhar_i

        allocate (recov_coeffs(mid))
 
        recov_coeffs(2:mid) = cmplx(sol(2:mid), sol(mid+1:size(sol))) / 2.0
        recov_coeffs(1) = cmplx(sol(1))

        recov_coeffs_ij = reshape(recov_coeffs,(/nhar_i,nhar_j/))

        topo_obj%nhar_i = nhar_i
        topo_obj%nhar_j = nhar_j
        topo_obj%fcoeffs = recov_coeffs_ij
    end subroutine recover_coeffs

end module fourier_mod