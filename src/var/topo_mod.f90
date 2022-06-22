module topo_mod
    use, intrinsic :: iso_fortran_env, only : error_unit
    use :: error_status, only : ALLOCATION_ERR
    implicit none

    private

    public :: topo_t, set_attrs, dealloc_obj

    type :: topo_t
        private
        real, dimension(:), allocatable :: lat, lon, topo
        logical, dimension(:,:,:), allocatable :: indices

    end type topo_t

contains

    subroutine set_attrs(lat, lon, clat, clon, width, topo, obj)
        implicit none
        type(topo_t), intent(out) :: obj
        real, dimension(:,:), intent(in) :: lat, lon
        real, dimension(:,:), allocatable :: latg, long
        logical :: lat_indices(size(lat,dim=1), size(lat,dim=2)), lon_indices(size(lon,dim=1), size(lon,dim=2))
        logical, dimension(:), allocatable :: tmp_lat, tmp_lon
        logical, dimension(:,:), allocatable :: tmp_indices
        real, intent(in) :: clat, clon, width
        real, dimension(:,:,:), intent(inout) :: topo
        real, dimension(1000000), allocatable :: tmp_lat_f, tmp_lon_f
        integer :: nrecs, nlat, nlon
        integer :: i, j, stat

        ! get the outer-most axis of the lat-lon grid
        nlat = size(lat, dim=1)
        nlon = size(lon, dim=1)

        ! get the number of records
        nrecs = size(topo,dim=3)

        allocate (obj%indices(nlon, nlat, nrecs), stat=stat)
        if (stat /= 0) then
            write(unit=error_unit, fmt='(A)') "Error allocating indices array for topo_obj"
            stop ALLOCATION_ERR
        end if

        ! we suppose that the allocation of the temporary arrays will not throw an error.
        allocate (tmp_lat(nlat))
        allocate (tmp_lon(nlon))
        allocate (tmp_indices(nlon, nlat))

        do i=1,nrecs
            print *, "nrec = ", i
            latg = spread(lat(:,i), 1, nlon)
            long = spread(lon(:,i), 2, nlat)
            tmp_indices = (abs(latg - clat) <= width)
            tmp_indices = tmp_indices .and. (abs(long - clon) <= width)
            
            tmp_lat = .false.
            do j=1,nlon
                tmp_lat = tmp_lat .or. tmp_indices(j,:)
            end do
            lat_indices(:,i) = tmp_lat

            tmp_lon = .false.
            do j=1,nlat
                tmp_lon = tmp_lon .or. tmp_indices(:,j)
            end do
            lon_indices(:,i) = tmp_lon

            ! Why do we need to flip the latitude axis?
            obj%indices(:,:,i) = tmp_indices(:,ubound(tmp_indices,dim=2):lbound(tmp_indices,dim=2):-1)
        end do
 
        deallocate(tmp_lat)
        deallocate(tmp_lon)
        deallocate(tmp_indices)

        obj%topo = pack(topo, obj%indices)
        obj%lat  = pack(lat,  lat_indices)
        obj%lon  = pack(lon,  lon_indices)

        print *, minval(obj%topo), maxval(obj%topo)

        print *, "obj%topo size: ", size(obj%topo)
        print *, "obj%lat size: ", size(obj%lat)
        print *, "obj%lon size: ", size(obj%lon)
        if ((size(obj%lat) * size(obj%lon)) /= size(obj%topo)) then
            write(unit=error_unit, fmt='(A)') "Error: Gathered subpoints shapes do not match."
        end if

    end subroutine set_attrs

    subroutine dealloc_obj(obj)
        implicit none
        type(topo_t), intent(inout) :: obj

        deallocate(obj%lat)
        deallocate(obj%lon)
        deallocate(obj%topo)
        deallocate(obj%indices)

    end subroutine dealloc_obj

end module topo_mod