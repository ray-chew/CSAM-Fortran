program orog_source
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, DP => real64
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use read_data_mod
    use write_data_mod
    use utils_mod
    use triangle_mod
    use topo_mod, only : topo_t, get_topo, dealloc_topo_obj, get_l2_err
    use fourier_mod, only : get_coeffs, nhar_i, nhar_j
    use lin_reg_mod, only : do_lin_reg
    use error_status, only : ALLOCATION_ERR
    use omp_lib

    implicit none
    character(len=128) :: fn, fn_linker, fn_topo, fn_output
    type(debug_t) :: debug_flags
    type(tol_t) :: tol_flags
    type(run_t) :: run_flags
    real, dimension(:), allocatable :: lat_center, lon_center, alphas, wls_lat, wls_lon, alphas_search, opt_deg, err_val
    integer, dimension(:,:), allocatable :: link_map
    real, dimension(:,:), allocatable :: lat_vert, lon_vert, topo_lat, topo_lon, coeffs, errs
    real(kind=DP), dimension(:,:,:), allocatable :: topo_dat
    real :: start, finish, wt_start, wt_finish
    integer :: i, j, Ncells, stat, ncid, nhi_dim_id, nhj_dim_id, ncell_dim_id, lat_dim_id = 0, lon_dim_id = 0
    integer :: ndegrees_dim_id, Ndegrees
    real :: clat, clon, nan, alpha
    complex, allocatable :: fcoeffs(:,:,:)

    type(topo_t) :: topo_obj
    type(llgrid_t) :: llgrid_obj
    logical, dimension(:,:), allocatable :: mask
    logical :: tmp_switch = .false.
    logical, dimension(:,:), allocatable :: cond

    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    Ndegrees = 90

    print *, "Reading grid and linker data..."
    call cpu_time(start)
    call get_fn(fn)
    call get_namelist(fn, fn_linker, fn_topo, fn_output, tol_flags, debug_flags, run_flags)

    call read_data(fn_linker, "clat", lat_center)
    call read_data(fn_linker, "clon", lon_center)
    call read_data(fn_linker, "clat_vertices", lat_vert)
    call read_data(fn_linker, "clon_vertices", lon_vert)
    call read_data(fn_linker, "links", link_map)

    call rad_to_deg(lat_center)
    call rad_to_deg(lon_center)
    call rad_to_deg(lat_vert)
    call rad_to_deg(lon_vert)

    print *, ""
    print *, "Reading topo data..."
    call read_data(fn_topo, "lat", topo_lat)
    call read_data(fn_topo, "lon", topo_lon)
    call read_data(fn_topo, "topo", topo_dat)

    call translate_deg_longitude(lon_center)
    call translate_deg_longitude(lon_vert)
    call translate_deg_longitude(topo_lon)

    print *, "Read topo_lat with shape: ", shape(topo_lat)
    print *, "Read topo_lon with shape: ", shape(topo_lon)
    print *, "Read topo_dat with shape: ", shape(topo_dat)

    call cpu_time(finish)
    print *, "Read time = ", finish-start ," seconds."
    print *, ""


    Ncells = size(lat_center)

    ! START I/O handling
    ncid = create_dataset(fn_output)
    nhi_dim_id = create_dim(ncid, 'nhar_i', nhar_i)
    nhj_dim_id = create_dim(ncid, 'nhat_j', nhar_j)
    ncell_dim_id = create_dim(ncid, 'ncell', Ncells)
    ndegrees_dim_id = create_dim(ncid, 'ndegrees', Ndegrees+1)

    ! if (debug_flags%output) then
        ! ...
    ! end if
    call close_dataset(ncid)
    ! END I/O handling

    if (run_flags%rotation == 1) then
        allocate (err_val(Ndegrees + 1))
        allocate (errs(Ndegrees + 1, Ncells))
        allocate (opt_deg(Ncells))
    end if


    allocate (fcoeffs(nhar_j,nhar_i,Ncells), stat=stat)
    if (stat /= 0) then
        write(unit=error_unit, fmt='(A)') "Error allocating nrecs array for topo_obj"
        stop ALLOCATION_ERR
    end if

    allocate (wls_lat(Ncells))
    allocate (wls_lon(Ncells))

    !$OMP PARALLEL
    !$OMP SINGLE
    print *, "Number of OMP threads used: ", omp_get_num_threads()
    call omp_set_num_threads(6)
    print *, ""
    !$OMP END SINGLE
    !$OMP END PARALLEL


    call cpu_time(start)
    wt_start = omp_get_wtime()
    print *, "Entering meaty loop..."

    !$OMP  PARALLEL DO DEFAULT(PRIVATE)                                     &
    !$OMP& SHARED(topo_lat, topo_lon, topo_dat, lat_center, lon_center,     &
    !$OMP& lat_vert, lon_vert, link_map, fcoeffs, fn_output,                &
    !$OMP& errs, opt_deg, wls_lat, wls_lon)                                 &
    !$OMP& FIRSTPRIVATE(run_flags, tol_flags, debug_flags, Ndegrees, Ncells, err_val, nan) &
    !$OMP& SCHEDULE(DYNAMIC)
    do i = 1, Ncells
        if (debug_flags%verbose) print *, "Starting cell: ", i

        clat = lat_center(i)
        clon = lon_center(i)
        call get_box_width(lat_vert(:,i), lon_vert(:,i), llgrid_obj, tol_flags%box_padding)

        if (debug_flags%verbose) print *, "Getting topo..."

        call get_topo(topo_lat, topo_lon, llgrid_obj, topo_dat, topo_obj, link_map(:,i), i, tol_flags%sp_real)

        !$OMP CRITICAL
        if (debug_flags%output) then
                ncid = open_dataset(fn_output)
                lat_dim_id = create_dim(ncid, 'lat_' // trim(str(i)), size(topo_obj%lat))
                lon_dim_id = create_dim(ncid, 'lon_' // trim(str(i)), size(topo_obj%lon))
                stat = write_data(ncid, 'lat_grid_' // trim(str(i)), topo_obj%lat_grid, (/lon_dim_id,lat_dim_id/))
                stat = write_data(ncid, 'lon_grid_' // trim(str(i)), topo_obj%lon_grid, (/lon_dim_id,lat_dim_id/))
                stat = write_data(ncid, 'topo_' // trim(str(i)), topo_obj%topo, (/lon_dim_id,lat_dim_id/))
                call close_dataset(ncid)
        end if
        !$OMP END CRITICAL

        ! if ((maxval(topo_obj%topo)) < 1.0) then

        !     print *, "Skipping (H < 1.0) cell: ", i
        !     fcoeffs(:,:,i) = nan
        !     wls_lat(i) = 0.0
        !     wls_lon(i) = 0.0
        !     ! print *, "Below sea level cell: ", i

        !     if (run_flags%rotation == 1) then
        !         opt_deg(i) = nan
        !         err_val = nan
        !         errs(:,i) = err_val
        !     end if

        !     call dealloc_topo_obj(topo_obj, .false.)

        if (debug_flags%skip_four) then

            print *, "Skipping fourier computations for cell: ", i
            fcoeffs(:,:,i) = nan 
            wls_lat(i) = 0.0
            wls_lon(i) = 0.0

        else

            cond = (topo_obj%topo < 1.0)
            cond = (topo_obj%topo > 0.0) .and. cond

            if ((sum(pack(topo_obj%topo, cond)) < 10.0) .and. (maxval(topo_obj%topo) < 10.0) ) then
                topo_obj%topo = 0.0
            end if                

            if (debug_flags%verbose) print *, "Setting triangular vertices"
            call set_triangle_verts(llgrid_obj, lat_vert(:,i), lon_vert(:,i))

            if (debug_flags%verbose) print *, "Masking points in triangle"
            mask = points_in_triangle(topo_obj%lat_grid,topo_obj%lon_grid, llgrid_obj)

            if (debug_flags%verbose) print *, "Getting coefficients"
            ! call get_coeffs(topo_obj, mask, coeffs)
            if (run_flags%full_spectrum) then
                call get_coeffs(topo_obj, mask, coeffs)
            else
                ! handle how to rotate axial wavenumbers
                if (run_flags%rotation == 0) then

                    alphas = (/0.0/)

                else if (run_flags%rotation == 1) then

                    alphas_search = (/(j, j=0, Ndegrees, 1)/)

                    tmp_switch = debug_flags%recover_topo
                    debug_flags%recover_topo = .true.

                    do j = 1, size(alphas_search)
                        alpha = alphas_search(j)
                        call get_coeffs(topo_obj, mask, alpha, coeffs)
                        call do_lin_reg(coeffs, topo_obj, mask, i, .true., debug_flags, tol_flags)
                        err_val(j) = get_l2_err(topo_obj)
                    end do

                    debug_flags%recover_topo = tmp_switch
                    alphas = real(minloc(err_val)) - 1.0

                    opt_deg(i) = alphas(1)
                    errs(:,i) = err_val

                else if (run_flags%rotation == 2) then

                    ! alphas = get_pca_angle()
                    alphas = (/0.0/)
                    print *, "Method not implemented"

                else if (run_flags%rotation == 3) then

                    alphas = (/0.0/)
                    print *, "Method not implemented"

                end if
                do j = 1, size(alphas)
                    alpha = alphas(j)
                    call get_coeffs(topo_obj, mask, alpha, coeffs)
                end do
                
            end if

            if (debug_flags%verbose) print *, "Doing linear regression"
            call do_lin_reg(coeffs, topo_obj, mask, i, run_flags%full_spectrum, debug_flags, tol_flags)
            fcoeffs(:,:,i) = topo_obj%fcoeffs
            wls_lat(i) = topo_obj%wavelength_lat
            wls_lon(i) = topo_obj%wavelength_lon

            print *, "Completed cell: ", i

            !$OMP CRITICAL
            if ((debug_flags%output) .and. (debug_flags%recover_topo)) then
                ncid = open_dataset(fn_output)
                stat = write_data(ncid, 'topo_recon_' // trim(str(i)), topo_obj%topo_recon_2D, (/lon_dim_id,lat_dim_id/))
                stat = write_data(ncid, 'mask_' // trim(str(i)), bool2int(mask), (/lon_dim_id,lat_dim_id/))
                call close_dataset(ncid)
            end if
            !$OMP END CRITICAL

            call dealloc_topo_obj(topo_obj, .true.)
            deallocate(coeffs)
            deallocate(mask)
        end if
    end do
    !$OMP END PARALLEL DO

    call cpu_time(finish)
    wt_finish = omp_get_wtime()
    print *, "Done with meaty loop..."
    print *, ""
    print *, "CPU time taken = ", (finish - start)
    print *, "Wall time taken = ", (wt_finish - wt_start)

    print *, ""
    print *, "Writing data output..."

    ncid = open_dataset(fn_output)
    stat = write_data(ncid, 'fcoeffs', fcoeffs, (/nhj_dim_id,nhi_dim_id,ncell_dim_id/))
    stat = write_data(ncid, 'wavelength_lat', wls_lat, (/ncell_dim_id/))
    stat = write_data(ncid, 'wavelength_lon', wls_lon, (/ncell_dim_id/))
    if (run_flags%rotation == 1) then
        stat = write_data(ncid, 'opt_deg', opt_deg, (/ncell_dim_id/))
        stat = write_data(ncid, 'errs', errs, (/ndegrees_dim_id,ncell_dim_id/))
    end if
    call close_dataset(ncid)

    deallocate(lat_center)
    deallocate(lon_center)
    deallocate(lat_vert)
    deallocate(lon_vert)
    deallocate(link_map)

    deallocate(topo_lat)
    deallocate(topo_lon)
    deallocate(topo_dat)

    deallocate(fcoeffs)
    if (run_flags%rotation == 1) then
        deallocate(opt_deg)
        deallocate(errs)
    end if

end program orog_source
