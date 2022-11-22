program grid_linker
    use read_data_mod
    use write_data_mod
    use utils_mod
    use topo_mod
    use omp_lib
    use, intrinsic :: iso_fortran_env, only : error_unit

    implicit none
    character(len=1024)                        :: fn
    character(len=1024)                        :: fn_grid
    character(len=1024)                        :: fn_topo
    type(flags_t)                           :: flags
    integer                                 :: Ncells, Nrecs
    integer                                 :: i
    integer, dimension(:,:), allocatable    :: icon_topo_links

    real, dimension(:), allocatable :: lat_center, lon_center
    real, dimension(:,:), allocatable :: topo_lat, topo_lon
    real :: clat, clon
    real :: start, finish, wt_start, wt_finish

    integer :: ncid, nrec_dim_id, ncell_dim_id, link_var_id

    call omp_set_num_threads(4)

    print *, "Start preprocessing grid data..."
    call get_fn(fn)
    fn = trim(fn)
    call get_namelist(fn, fn_grid, fn_topo, flags)
    fn_grid = trim(fn_grid)
    fn_topo = trim(fn_topo)

    call read_data(fn_grid, "clat", lat_center)
    call read_data(fn_grid, "clon", lon_center)
    ! call read_data(fn_grid, "clat_vertices", lat_vert)
    ! call read_data(fn_grid, "clon_vertices", lon_vert)

    call rad_to_deg(lat_center)
    call rad_to_deg(lon_center)
    ! call rad_to_deg(lat_vert)
    ! call rad_to_deg(lon_vert)

    if (size(lat_center) /= size(lon_center)) then
        write(unit = error_unit, fmt='(2A)') "Size of lat-lon entries for ICON grid are unequal."
    end if

    print *, "Reading topo data..."
    call read_data(fn_topo, "lat", topo_lat)
    call read_data(fn_topo, "lon", topo_lon)
    ! call read_data(fn_topo, "topo", topo_dat)

    print *, "Read topo_lat with shape: ", shape(topo_lat)
    print *, "Read topo_lon with shape: ", shape(topo_lon)
    ! print *, "Read topo_dat with shape: ", shape(topo_dat)

    Ncells = size(lat_center)
    Ncells = 4
    Nrecs = size(topo_lat, dim=2)
    !!! IN FORTRAN, OUTERMOST INDEX IS THE FASTEST !!!
    allocate (icon_topo_links(Nrecs,Ncells))
    icon_topo_links = 0

    call cpu_time(start)
    wt_start = omp_get_wtime()

    !$OMP PARALLEL DO
    do i = 1, Ncells
        ! print *, "Ncell = ", i
        clat = lat_center(i)
        clon = lon_center(i)
        call get_topo_idx(topo_lat, topo_lon, clat, clon, 2.0, i, icon_topo_links)
    end do
    !$OMP END PARALLEL DO

    call cpu_time(finish)
    wt_finish = omp_get_wtime()
    print '("CPU time taken = ",f6.3," seconds.")', finish-start
    print '("Wall time taken = ",f6.3," seconds.")', wt_finish-wt_start

    print *, "End preprocessing grid data..."
    print *, "Writing preprocessed grid data..."
    ncid = open_dataset(fn_grid)
    nrec_dim_id = create_dim(ncid, 'nrecs', Nrecs)
    ncell_dim_id = get_dim_id(ncid, 'cell')
    link_var_id = write_data(ncid, 'links', icon_topo_links, (/nrec_dim_id,ncell_dim_id/))
    call close_dataset(ncid)


    print *, "Peach is the best cat ever"

end program grid_linker