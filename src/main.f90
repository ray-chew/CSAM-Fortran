program main
    use read_grid_mod
    use utils_mod
    use topo_mod

    implicit none
    character(len=1024) :: fn_grid, fn_topo
    real, dimension(:), allocatable :: lat_center, lon_center
    real, dimension(:,:), allocatable :: lat_vert, lon_vert, topo_lat, topo_lon
    real, dimension(:,:,:), allocatable :: topo_dat
    real :: lat_ref, lon_ref

    real :: start, finish

    type(topo_t) :: topo_obj

    ! hard-coded index for test.
    integer, parameter :: ref_idx = 441 !13680

    print *, "Reading grid data..."

    call get_fn(fn_grid, fn_topo)
    fn_grid = trim(fn_grid)
    fn_topo = trim(fn_topo)

    call get_data(fn_grid, "clat", lat_center)
    call get_data(fn_grid, "clon", lon_center)
    call get_data(fn_grid, "clat_vertices", lat_vert)
    call get_data(fn_grid, "clon_vertices", lon_vert)

    ! lat_vert = transpose(lat_vert)
    ! lon_vert = transpose(lon_vert)

    call rad_to_deg(lat_center)
    call rad_to_deg(lon_center)

    lat_ref = lat_center(ref_idx)
    lon_ref = lon_center(ref_idx)

    print *, "Reference (lat, lon): ", lat_ref, lon_ref

    print *, "Reading topo data..."
    call get_data(fn_topo, "lat", topo_lat)
    call get_data(fn_topo, "lon", topo_lon)
    call get_data(fn_topo, "topo", topo_dat)

    ! topo_lat = transpose(topo_lat)
    ! topo_lon = transpose(topo_lon)
    ! topo_dat = transpose(topo_dat)

    print *, "Read topo_lat with shape: ", shape(topo_lat)
    print *, "Read topo_lon with shape: ", shape(topo_lon)
    print *, "Read topo_dat with shape: ", shape(topo_dat)

    print *, "Gathering subpoints..."
    call cpu_time(start)
    call set_attrs(topo_lat, topo_lon, lat_ref, lon_ref, 2.0, topo_dat, topo_obj)
    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start

    call dealloc_all()
    call dealloc_obj(topo_obj)

contains

    subroutine dealloc_all()
    implicit none
        deallocate(lat_center)
        deallocate(lon_center)
        deallocate(lat_vert)
        deallocate(lon_vert)
    
        deallocate(topo_lat)
        deallocate(topo_lon)
        deallocate(topo_dat)
    end subroutine dealloc_all

end program main