module fourier_mod
    use :: topo_mod, only : topo_t
    use :: utils_mod, only : get_N_unique, deg_to_rad
    use :: const_mod, only : PI
    use :: stdlib_sorting, only : ord_sort
    implicit none

    private

    public :: get_coeffs, get_axial_coeffs, get_full_coeffs, recover_coeffs, nhar_i, nhar_j

    integer, parameter ::   nhar_i = 12, &
                            nhar_j = 12

    type :: ftmp_t
        integer :: Ni, Nj
        integer, dimension(:), allocatable :: m_i, m_j
        integer, dimension(:), allocatable :: II, JJ
        real, dimension(:,:), allocatable :: terms_i, terms_j

    contains

        procedure, private :: prepare_terms

    end type ftmp_t

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

        N_cos = nhar_i * nhar_j - nhar_j/2
        N_sin = nhar_i * nhar_j - nhar_j/2 - 1

        allocate (tmp(nhar_i * nhar_j))
        allocate (coeffs(N_cos + N_sin, size(topo_tri)))

        allocate (c_cos(N_cos))
        allocate (c_sin(N_sin))

        do k=1,size(topo_tri)
            l = 1
            m = 1
            n = 1
            do i=0,nhar_i-1
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

                    l = l + 1
                end do
            end do

            coeffs(1:N_cos,k) = c_cos
            coeffs(N_cos+1:N_cos+N_sin,k) = c_sin

        end do

        topo_obj%lat_tri = lat_tri
        topo_obj%lon_tri = lon_tri
        topo_obj%topo_tri = topo_tri

    end subroutine get_coeffs


    subroutine get_full_coeffs(topo_obj, mask, coeffs)
        implicit none
        type(topo_t), intent(inout) :: topo_obj
        real, dimension(:,:), allocatable, intent(out) :: coeffs
        logical, dimension(:,:), intent(in) :: mask
        real, dimension(:,:), allocatable :: c_cos, c_sin, terms_sum
        real, dimension(:,:,:), allocatable :: tmp_i, tmp_j

        type(ftmp_t) :: f_obj
        integer :: ncells, dense_sz, nhj_0

        call get_IJ(mask, f_obj, topo_obj)
        call f_obj%prepare_terms()

        ! Note, I am storing terms_i and terms_j as row major here...
        f_obj%terms_i = spread(real(f_obj%m_i), dim=1, ncopies=size(f_obj%II)) * & 
        spread(real(f_obj%II), dim=2, ncopies=size(f_obj%m_i)) / real(f_obj%Ni)
        f_obj%terms_j = spread(real(f_obj%m_j), dim=1, ncopies=size(f_obj%JJ)) * &
        spread(real(f_obj%JJ), dim=2, ncopies=size(f_obj%m_j)) / real(f_obj%Nj)

        ncells = size(f_obj%terms_i, dim=1)
        dense_sz = size(f_obj%terms_i, dim=2) * size(f_obj%terms_j, dim=2)

        f_obj%terms_i = reshape(spread(f_obj%terms_i, dim=2, ncopies=size(f_obj%terms_j, dim=2)), &
        (/ncells, dense_sz/), order=(/1,2/))
        f_obj%terms_j = reshape(spread(f_obj%terms_j, dim=3, ncopies=size(f_obj%terms_i, dim=2)), &
        (/ncells, dense_sz/), order=(/1,2/))

        ! tmp_i = spread(f_obj%terms_i, dim=2, ncopies=size(f_obj%terms_j, dim=2))
        ! tmp_j = spread(f_obj%terms_j, dim=3, ncopies=size(f_obj%terms_i, dim=2))

        ! f_obj%terms_i = reshape(tmp_i, (/ncells, dense_sz/), order=(/1,2/))
        ! f_obj%terms_j = reshape(tmp_j, (/ncells, dense_sz/), order=(/1,2/))

        nhj_0 = nhar_j / 2 + 1 ! index position of the zeroth j-wavenumber
        terms_sum = f_obj%terms_i + f_obj%terms_j

        c_cos = cos(2.0 * PI * terms_sum(:,nhj_0:dense_sz))
        c_sin = sin(2.0 * PI * terms_sum(:,nhj_0+1:dense_sz))

        ! here, we do a transpose of coeffs to return a column-major result.
        allocate (coeffs(size(c_cos, dim=2) + size(c_sin, dim=2), ncells))

        coeffs(1 : size(c_cos, dim=2),:) = transpose(c_cos)
        coeffs(size(c_cos, dim=2)+1 : size(c_cos, dim=2) + size(c_sin, dim=2), :) = transpose(c_sin)

    end subroutine get_full_coeffs


    subroutine get_axial_coeffs(topo_obj, mask, deg_alpha, coeffs)
        implicit none
        type(topo_t), intent(inout) :: topo_obj
        real, intent(inout) :: deg_alpha
        real, dimension(:,:), allocatable, intent(out) :: coeffs
        logical, dimension(:,:), intent(in) :: mask
        real, dimension(:), allocatable :: ktil, ltil, khat, lhat
        real, dimension(:,:), allocatable :: c_cos, c_sin
        real :: alpha
        ! real, dimension(:,:), allocatable :: tmp

        type(ftmp_t) :: f_obj
        integer :: nhj_0

        ! do rotation of topo here...

        call get_IJ(mask, f_obj, topo_obj)
        call f_obj%prepare_terms()

        alpha = deg_to_rad(deg_alpha)

        ! Note that I am initialising the arrays as row-majors here.
        ktil = real(f_obj%m_i) * cos(alpha)
        ltil = real(f_obj%m_i) * sin(alpha)

        f_obj%terms_i = spread(ktil, dim=1, ncopies=size(f_obj%II)) * spread(f_obj%II, dim=2, ncopies=size(ktil)) / f_obj%Ni &
        + spread(ltil, dim=1, ncopies=size(f_obj%JJ)) * spread(f_obj%JJ, dim=2, ncopies=size(ltil)) / f_obj%Nj


        khat = real(f_obj%m_j) * cos(alpha + PI/2.0)
        lhat = real(f_obj%m_j) * sin(alpha + PI/2.0)

        f_obj%terms_j = spread(khat, dim=1, ncopies=size(f_obj%II)) * spread(f_obj%II, dim=2, ncopies=size(khat)) / f_obj%Ni &
        + spread(lhat, dim=1, ncopies=size(f_obj%JJ)) * spread(f_obj%JJ, dim=2, ncopies=size(lhat)) / f_obj%Nj

        nhj_0 = nhar_j / 2 + 1 ! position of zeroth j-wavenumber

        allocate (c_cos(size(f_obj%terms_i, dim=1), nhar_i + int(nhar_j - nhj_0)))

        c_cos(:, 1:nhar_i) = f_obj%terms_i(:,:)
        c_cos(:, int(nhar_i+1):int(nhar_i+1) + int(nhar_j - nhj_0 - 1)) = f_obj%terms_j(:, int(nhj_0+1):nhar_j)

        c_cos = 2.0 * cos(2.0 * PI * c_cos)


        allocate (c_sin(size(f_obj%terms_i, dim=1), int(nhar_i-1) + int(nhar_j - nhj_0)))

        c_sin(:, 1:nhar_i-1) = f_obj%terms_i(:,2:nhar_i)
        c_sin(:, nhar_i:nhar_i + int(nhar_j - nhj_0 - 1)) = f_obj%terms_j(:, int(nhj_0+1):nhar_j)

        c_sin = 2.0 * sin(2.0 * PI * c_sin)


        allocate (coeffs(size(c_cos, dim=2) + size(c_sin,dim=2), size(f_obj%terms_i, dim=1)))

        ! We do a transpose here to make our coefficients column-major again
        coeffs(1:size(c_cos,dim=2),:) = transpose(c_cos)
        coeffs(size(c_cos,dim=2)+1:size(c_cos,dim=2)+size(c_sin,dim=2),:) = transpose(c_sin)

    end subroutine

    subroutine get_IJ(mask, f_obj, topo_obj)
        implicit none
        type(topo_t), intent(inout) :: topo_obj
        type(ftmp_t), intent(inout) :: f_obj
        logical, dimension(:,:), intent(in) :: mask

        ! integer, dimension(count(mask)), intent(out) :: II, JJ
        real, dimension(count(mask)) :: lat_tri, lon_tri, topo_tri
        real :: d_lat, d_lon

        ! we get lat, lon and topo in the triangle using the mask
        lat_tri = pack(topo_obj%lat_grid, mask=mask)
        lon_tri = pack(topo_obj%lon_grid, mask=mask)
        topo_tri = pack(topo_obj%topo, mask=mask)

        topo_obj%lat_tri = lat_tri
        topo_obj%lon_tri = lon_tri
        topo_obj%topo_tri = topo_tri

        ! get grid spacing assuming equidistant grid.
        d_lat = topo_obj%lat(2) - topo_obj%lat(1)
        d_lon = topo_obj%lon(2) - topo_obj%lon(1)

        f_obj%II = ceiling((lat_tri - minval(lat_tri)) / d_lat)
        f_obj%JJ = ceiling((lon_tri - minval(lon_tri)) / d_lon)

    end subroutine get_IJ


    subroutine prepare_terms(self)
        implicit none
        ! type(topo_t), intent(inout) :: topo_obj
        class(ftmp_t), intent(inout) :: self
        ! logical, intent(in) :: compute_terms
        ! integer, dimension(count(mask)), intent(out) :: II, JJ
        integer, dimension(size(self%JJ)) :: JJ_tmp
        integer :: i

        self%Ni = get_N_unique(self%II)
        JJ_tmp = self%JJ
        call ord_sort(JJ_tmp)
        self%Nj = get_N_unique(JJ_tmp)

        self%m_i = (/(i, i=0, nhar_i-1, 1)/)

        if (mod(nhar_j, 2) == 0) then
            self%m_j = (/(i, i=-nhar_j/2, nhar_j/2-1, 1)/)
        else
            self%m_j = (/(i, i=-(nhar_j-1)/2, nhar_j/2-1, 1)/)
        end if

        ! if (compute_terms) then
        ! end if
        
    end subroutine prepare_terms


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