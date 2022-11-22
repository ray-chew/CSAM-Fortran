module topo_mod
    use, intrinsic :: iso_fortran_env, only : error_unit, DP => real64
    use :: error_status, only : ALLOCATION_ERR
    implicit none

    private

    public :: topo_t, get_topo, get_topo_idx, dealloc_topo_obj

    type :: topo_t
        real, dimension(:), allocatable :: lat, lon, lat_tri, lon_tri, topo_tri
        real, dimension(:,:), allocatable :: topo, lat_grid, lon_grid, topo_recon_2D

        integer :: nhar_i, nhar_j
        complex, dimension(:,:), allocatable :: fcoeffs
    end type topo_t

contains

    subroutine get_topo_idx(lat, lon, clat, clon, width, ncell, icon_topo_links)
        implicit none
        real, dimension(:,:), intent(in) :: lat
        real, dimension(:,:), intent(in) :: lon
        real, intent(in) :: clat
        real, intent(in) :: clon
        real, intent(in) :: width
        integer, intent(in) :: ncell
        integer, dimension(:,:), intent(out) :: icon_topo_links

        integer :: nrecs, nlat, nlon
        integer :: i, j, stat
        logical, dimension(:,:), allocatable :: tmp_indices
        logical :: lat_indices(size(lat,dim=1), size(lat,dim=2)), lon_indices(size(lon,dim=1), size(lon,dim=2))
        real, dimension(:,:), allocatable :: lat_grid, lon_grid

        ! get the outer-most axis of the lat-lon grid
        nlat = size(lat, dim=1)
        nlon = size(lon, dim=1)

        ! get the number of records
        nrecs = size(lat,dim=2)

        lat_indices = .false.
        lon_indices = .false.

        ! print *, count(lat_indices)
        ! print *, count(lon_indices)

        allocate (tmp_indices(nlon, nlat))

        j = 1
        do i=1,nrecs
            ! print *, "nrec = ", i
            lat_grid = spread(lat(:,i), 1, nlon)
            lon_grid = spread(lon(:,i), 2, nlat)

            tmp_indices = (abs(lat_grid - clat) <= width)
            tmp_indices = tmp_indices .and. (abs(lon_grid - clon) <= width)

            ! print *, count(tmp_indices)

            if (count(tmp_indices) > 0) then
                icon_topo_links(j,ncell) = i
                j = j + 1
            end if
        end do


    end subroutine get_topo_idx

    subroutine get_topo(lat, lon, clat, clon, width, topo, obj)
        implicit none
        type(topo_t), intent(out) :: obj
        real, dimension(:,:), intent(in) :: lat, lon      
        real, intent(in) :: clat, clon, width
        real(kind=DP), dimension(:,:,:), intent(inout) :: topo

        integer :: nrecs, nlat, nlon
        integer :: i, j, stat
        logical, dimension(:), allocatable :: tmp_lat, tmp_lon
        logical, dimension(:,:), allocatable :: tmp_indices
        logical, dimension(:,:,:), allocatable :: tmp_indices_3D
        logical :: lat_indices(size(lat,dim=1), size(lat,dim=2)), lon_indices(size(lon,dim=1), size(lon,dim=2))
        real, dimension(:), allocatable :: tmp_topo
        real, dimension(:,:), allocatable :: lat_grid, lon_grid, topo_grid

        ! get the outer-most axis of the lat-lon grid
        nlat = size(lat, dim=1)
        nlon = size(lon, dim=1)

        lat_indices = .false.
        lon_indices = .false.

        ! get the number of records
        nrecs = size(topo,dim=3)
        ! nrecs = 27

        lat_indices = .false.
        lon_indices = .false.

        ! print *, count(lat_indices)
        ! print *, count(lon_indices)

        ! we suppose that the allocation of the temporary 2D arrays will not throw an error.
        allocate (tmp_lat(nlat))
        allocate (tmp_lon(nlon))
        allocate (tmp_indices(nlon, nlat))

        allocate (tmp_indices_3D(nlon, nlat, nrecs), stat=stat)
        tmp_indices_3D = .false.
        if (stat /= 0) then
            write(unit=error_unit, fmt='(A)') "Error allocating indices array for topo_obj"
            stop ALLOCATION_ERR
        end if

        ! do i=1,nrecs
        do i=23,26
            print *, "nrec = ", i
            lat_grid = spread(lat(:,i), 1, nlon)
            lon_grid = spread(lon(:,i), 2, nlat)
            topo_grid = topo(:,:,i)

            tmp_indices = (abs(lat_grid - clat) <= width)
            tmp_indices = tmp_indices .and. (abs(lon_grid - clon) <= width)

            print *, count(tmp_indices)
            
            tmp_lat = .false.
            do j=1,nlon
                tmp_lat = tmp_lat .or. tmp_indices(j,:)
            end do
            lat_indices(:,i) = tmp_lat

            print *, count(tmp_lat)

            tmp_lon = .false.
            do j=1,nlat
                tmp_lon = tmp_lon .or. tmp_indices(:,j)
            end do
            lon_indices(:,i) = tmp_lon

            print *, count(tmp_lon)

            ! Why do we need to flip the latitude axis?
            tmp_indices_3D(:,:,i) = tmp_indices(:,ubound(tmp_indices,dim=2):lbound(tmp_indices,dim=2):-1)
        end do

        obj%lat  = pack(lat,  lat_indices)
        obj%lon  = pack(lon,  lon_indices)

        print *, "obj%lat size: ", size(obj%lat)
        print *, "obj%lon size: ", size(obj%lon)

        allocate (tmp_topo(size(obj%lat) * size(obj%lon)))
        
        tmp_topo = pack(topo, tmp_indices_3D)
        
        obj%topo = reshape(tmp_topo, shape=(/size(obj%lon),size(obj%lat)/), order=(/1,2/))
        obj%topo = obj%topo(:,ubound(obj%topo,dim=2):lbound(obj%topo,dim=2):-1)

        if ((size(obj%lat) * size(obj%lon)) /= size(obj%topo)) then
            write(unit=error_unit, fmt='(A)') "Error: Gathered subpoints shapes do not match."
        end if

        obj%lat_grid = spread(obj%lat, 1, size(obj%lon))
        obj%lon_grid = spread(obj%lon, 2, size(obj%lat))
    end subroutine get_topo

    subroutine dealloc_topo_obj(obj)
        implicit none
        type(topo_t), intent(inout) :: obj

        deallocate(obj%lat)
        deallocate(obj%lon)
        deallocate(obj%lat_grid)
        deallocate(obj%lon_grid)
        deallocate(obj%topo)

        deallocate(obj%lat_tri)
        deallocate(obj%lon_tri)
        deallocate(obj%topo_tri)

    end subroutine dealloc_topo_obj

end module topo_mod