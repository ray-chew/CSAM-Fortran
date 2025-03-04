module lin_reg_mod
    use :: topo_mod, only : topo_t
    use, intrinsic :: iso_fortran_env, only : error_unit
    use :: error_status, only : LINALG_ERR
    use :: fourier_mod, only : recover_coeffs, recover_sol
    use :: utils_mod, only : debug_t, tol_t

    implicit none

    private

    public :: do_lin_reg

contains

    subroutine do_lin_reg(coeffs, topo_obj, mask, ncell, full_spectrum, debug, tol)
        implicit none
        real, dimension(:,:), intent(in) :: coeffs
        type(topo_t), intent(inout) :: topo_obj
        integer, intent(in) :: ncell

        real, dimension(size(coeffs,dim=1)) :: h_hat, sol
        real, dimension(size(coeffs,dim=1),size(coeffs,dim=1)) :: M, Minv
        integer :: nc, nd, istat, nwork
        integer, dimension(size(coeffs,dim=1)) :: ipiv
        real, dimension(size(coeffs,dim=1)) :: work
        real, dimension(size(coeffs,dim=2)) :: z_recon

        logical, intent(in) :: mask(:,:)
        type(debug_t), intent(in) :: debug
        type(tol_t), intent(in) :: tol
        logical, intent(in) :: full_spectrum

        nc = size(coeffs,dim=1) ! wavenumber grid
        nd = size(coeffs,dim=2) ! number of cells
        nwork = size(coeffs,dim=1)

        ! h_hat = matmul(coeffs, topo_obj%topo_tri)
        ! M = matmul(coeffs, transpose(coeffs))

        call dgemm('n','n', nc, 1, nd, 1.0, coeffs, nc, topo_obj%topo_tri, nd, 0.0, h_hat, nc)
        call dgemm('n','t', nc, nc, nd, 1.0, coeffs, nc, coeffs, nc, 0.0, M, nc)

        ! could just replace M...
        Minv = M

        call dgetrf(nc,nc,Minv,nc,ipiv,istat)

        if (istat /= 0) then
            write(unit=error_unit, fmt=*) "Triangulation of matrix unsuccessful for cell: ", ncell, "! Error code: ", istat
            stop LINALG_ERR
        end if

        call dgetri(nc,Minv,nc,ipiv,work,nwork,istat)

        if (istat /= 0) then
            write(unit=error_unit, fmt=*) "Inversion of matrix unsuccessful for cell: ", ncell, "! Error code: ", istat 
            stop LINALG_ERR
        end if

        call dgemm('n','n', nc, 1, nc, 1.0, Minv, nc, h_hat, nc, 0.0, sol, nc)
        call recover_coeffs(topo_obj, sol, full_spectrum)

        if (debug%recover_topo) then
            if (debug%recover_sol) then
                if (debug%verbose) print *, "recovering sol from fcoeffs..."
                call recover_sol(topo_obj, sol, tol, full_spectrum)
            end if

            if (debug%verbose) print *, "recovering topography from spectrum"
            call dgemm('t','n', nd, 1, nc, 1.0, coeffs, nc, sol, nc, 0.0, z_recon, nd)
            call recon_2D(topo_obj, z_recon, mask)
        end if

    end subroutine

    subroutine recon_2D(topo_obj, z_recon, mask)
        implicit none
        type(topo_t), intent(inout) :: topo_obj
        real, intent(in) :: z_recon(:)
        logical, intent(in) :: mask(:,:)
        real, allocatable :: z_recon_2D(:,:)

        integer :: i, j, k

        allocate (z_recon_2D(size(topo_obj%lon),size(topo_obj%lat)))
        z_recon_2D = 0.0

        k = 1
        do j=1,size(mask,2)
            do i=1,size(mask,1)
                if (mask(i,j) .eqv. .true.) then
                    z_recon_2D(i,j) = z_recon(k)
                    k = k+1
                end if
            end do
        end do
        topo_obj%topo_recon_2D = z_recon_2D

        ! print *, shape(z_recon_2D)
    end subroutine recon_2D

end module lin_reg_mod