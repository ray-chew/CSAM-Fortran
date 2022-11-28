program orog_source
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, DP => real64
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use read_data_mod
    use write_data_mod
    use utils_mod
    use triangle_mod
    use topo_mod, only : topo_t, get_topo, dealloc_topo_obj
    use fourier_mod, only : get_coeffs, nhar_i, nhar_j
    use lin_reg_mod, only : do_lin_reg
    use error_status, only : ALLOCATION_ERR
    use omp_lib

    implicit none
    character(len=128) :: fn, fn_linker, fn_topo, fn_output
    type(debug_t) :: debug_flags
    type(tol_t) :: tol_flags
    real, dimension(:), allocatable :: lat_center, lon_center
    integer, dimension(:,:), allocatable :: link_map
    real, dimension(:,:), allocatable :: lat_vert, lon_vert, topo_lat, topo_lon, coeffs
    real(kind=DP), dimension(:,:,:), allocatable :: topo_dat
    real :: start, finish, wt_start, wt_finish
    integer :: ref_idx, i, Ncells, stat, ncid, nhi_dim_id, nhj_dim_id, ncell_dim_id, lat_dim_id, lon_dim_id
    real :: clat, clon, nan
    complex, allocatable :: fcoeffs(:,:,:)

    type(topo_t) :: topo_obj
    type(llgrid_t) :: llgrid_obj
    logical, dimension(:,:), allocatable :: mask
    integer, parameter :: chunk=100

    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

    print *, "Reading grid and linker data..."
    call cpu_time(start)
    call get_fn(fn)
    call get_namelist(fn, fn_linker, fn_topo, fn_output, tol_flags, debug_flags)

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

    print *, "Read topo_lat with shape: ", shape(topo_lat)
    print *, "Read topo_lon with shape: ", shape(topo_lon)
    print *, "Read topo_dat with shape: ", shape(topo_dat)

    call cpu_time(finish)
    print '("Read time = ",f6.3," seconds.")',finish-start

    ! Snippet of code for tests...topo_obj
    ref_idx = 1000

    clat = lat_center(ref_idx)
    clon = lon_center(ref_idx)
    ! End Snippet

    Ncells = size(lat_center)
    ! Ncells = 1
    ! width = 2.0
    i = 1

    ! START I/O handling
    ncid = create_dataset(fn_output)
    nhi_dim_id = create_dim(ncid, 'nhar_i', nhar_i)
    nhj_dim_id = create_dim(ncid, 'nhat_j', nhar_j)
    ncell_dim_id = create_dim(ncid, 'ncell', Ncells)

    ! if (debug_flags%output) then
        ! ...
    ! end if
    call close_dataset(ncid)
    ! END I/O handling


    allocate (fcoeffs(nhar_j,nhar_i,Ncells), stat=stat)
    if (stat /= 0) then
        write(unit=error_unit, fmt='(A)') "Error allocating nrecs array for topo_obj"
        stop ALLOCATION_ERR
    end if

    !$OMP PARALLEL
    !$OMP SINGLE
    print *, "Number of OMP threads used: ", omp_get_num_threads()
    !$OMP END SINGLE
    !$OMP END PARALLEL


    call cpu_time(start)
    wt_start = omp_get_wtime()
    print *, "Entering meaty loop..."

    !$OMP  PARALLEL DO SHARED(topo_lat, topo_lon, topo_dat, lat_center, lon_center, lat_vert, lon_vert, link_map,       &
    !$OMP& fcoeffs, fn_output)          &
    !$OMP& DEFAULT(PRIVATE)             &
    !$OMP& FIRSTPRIVATE(tol_flags, debug_flags)
    !$OMP& SCHEDULE(DYNAMIC, chunk)
    do i = 1, Ncells
        ! print *, "Starting cell: ", i
        clat = lat_center(i)
        clon = lon_center(i)
        call get_box_width(lat_vert(:,i), lon_vert(:,i), llgrid_obj, tol_flags%box_padding)

        ! print *, "Getting topo..."
        call get_topo(topo_lat, topo_lon, clat, clon, llgrid_obj, topo_dat, topo_obj, link_map(:,i), i, tol_flags%sp_real)

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

        if ((maxval(topo_obj%topo)) < 1.0) then

            print *, "Skipping ocean cell: ", i
            !OMP CRITICAL
            fcoeffs(:,:,i) = nan
            !OMP END CRITICAL

        else if (debug_flags%skip_four) then

            print *, "Skipping fourier computations cell: ", i
            !OMP CRITICAL
            fcoeffs(:,:,i) = nan 

        else

            ! print *, "Setting triangular vertices"
            call set_triangle_verts(llgrid_obj, lat_vert(:,i), lon_vert(:,i))
            ! print *, "Masking points in triangle"
            mask = points_in_triangle(topo_obj%lat_grid,topo_obj%lon_grid, llgrid_obj)

            ! print *, "Getting coefficients"
            call get_coeffs(topo_obj, mask, coeffs)
            ! print *, "Doing linear regression"
            call do_lin_reg(coeffs, topo_obj, mask, i, debug_flags%recover_topo)
            !OMP CRITICAL
            fcoeffs(:,:,i) = topo_obj%fcoeffs
            !OMP END CRITICAL

            print *, "Completed cell: ", i

            !$OMP CRITICAL
            if ((debug_flags%output) .and. (debug_flags%recover_topo)) then
                ncid = open_dataset(fn_output)
                stat = write_data(ncid, 'topo_recon_' // trim(str(i)), topo_obj%topo_recon_2D, (/lon_dim_id,lat_dim_id/))
                stat = write_data(ncid, 'mask_' // trim(str(i)), bool2int(mask), (/lon_dim_id,lat_dim_id/))
                call close_dataset(ncid)
            end if
            !$OMP END CRITICAL

        end if
    end do
    !$OMP END PARALLEL DO

    call cpu_time(finish)
    wt_finish = omp_get_wtime()
    print *, "Done with meaty loop..."
    print *, "CPU time taken = ", (finish - start)
    print *, "Wall time taken = ", (wt_finish - wt_start)

    print *, ""
    print *, "Writing data output..."

    ncid = open_dataset(fn_output)    
    stat = write_data(ncid, 'fcoeffs', fcoeffs, (/nhj_dim_id,nhi_dim_id,ncell_dim_id/))
    call close_dataset(ncid)

end program orog_source
