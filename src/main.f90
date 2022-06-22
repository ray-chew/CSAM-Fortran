program main
    use read_grid_mod
    use utils_mod

    implicit none
    character(len=1024) :: fn
    real, dimension(:), allocatable :: lat, lon
    real, dimension(:,:), allocatable :: latv, lonv
    real :: lat_ref, lon_ref

    ! hard-coded index for test.
    integer, parameter :: ref_idx = 13680

    call get_fn(fn)
    fn = trim(fn)

    call get_grid_data(fn, "clat", lat)
    call get_grid_data(fn, "clon", lon)
    call get_grid_data(fn, "clat_vertices", latv)
    call get_grid_data(fn, "clon_vertices", lonv)

    latv = transpose(latv)
    lonv = transpose(lonv)

    call rad_to_deg(lat)
    call rad_to_deg(lon)

    lat_ref = lat(ref_idx)
    lon_ref = lon(ref_idx)

    print *, lat_ref, lon_ref

    deallocate(lat)
    deallocate(lon)
    deallocate(latv)
    deallocate(lonv)

end program main