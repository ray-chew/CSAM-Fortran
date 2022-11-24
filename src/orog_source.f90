program orog_source
    use, intrinsic :: iso_fortran_env, only : error_unit, DP => real64
    use read_data_mod
    use write_data_mod
    use utils_mod
    use topo_mod, only : topo_t, get_topo, dealloc_topo_obj
    use fourier_mod, only : llgrid_t, set_triangle_verts, get_coeffs, points_in_triangle
    use lin_reg_mod, only : do_lin_reg

    implicit none
    character(len=1024) :: fn, fn_linker, fn_topo
    type(flags_t) :: flags
    real, dimension(:), allocatable :: lat_center, lon_center
    integer, dimension(:,:), allocatable :: link_map
    real, dimension(:,:), allocatable :: lat_vert, lon_vert, topo_lat, topo_lon, coeffs
    real(kind=DP), dimension(:,:,:), allocatable :: topo_dat
    real :: start, finish
    integer :: ref_idx, i, Ncells
    real :: clat, clon, width

    type(topo_t) :: topo_obj
    type(llgrid_t) :: llgrid_obj
    logical, dimension(:,:), allocatable :: mask

    print *, "Reading grid and linker data..."
    call cpu_time(start)
    call get_fn(fn)
    call get_namelist(fn, fn_linker, fn_topo, flags)

    call read_data(fn_linker, "clat", lat_center)
    call read_data(fn_linker, "clon", lon_center)
    call read_data(fn_linker, "clat_vertices", lat_vert)
    call read_data(fn_linker, "clon_vertices", lon_vert)
    call read_data(fn_linker, "links", link_map)

    call rad_to_deg(lat_center)
    call rad_to_deg(lon_center)
    call rad_to_deg(lat_vert)
    call rad_to_deg(lon_vert)

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
    ncells = 4
    width = 2.0

    call cpu_time(start)
    print *, "Entering meaty loop..."


    ! SHARED(topo_lat, topo_lon, topo_dat, topo_obj, link_map) PRIVATE(i, clat, clon, mask, coeffs) FIRSTPRIVATE(width)
    do i = 1, Ncells
        clat = lat_center(i)
        clon = lon_center(i)
        call get_topo(topo_lat, topo_lon, clat, clon, width, topo_dat, topo_obj, link_map(:,i), i)
        call set_triangle_verts(llgrid_obj, lat_vert(:,ref_idx), lon_vert(:,ref_idx))
        mask = points_in_triangle(topo_obj%lat_grid,topo_obj%lon_grid, llgrid_obj)
        call get_coeffs(topo_obj, mask, coeffs)
        call do_lin_reg(coeffs, topo_obj, mask)
    end do



end program orog_source