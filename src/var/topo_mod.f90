module topo_mod
    use, intrinsic :: iso_fortran_env, only : error_unit
    use :: error_status, only : ALLOCATION_ERR
    use :: stdlib_sorting, only : sort_index, int_size
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
        real, dimension(:,:), allocatable :: latg, long, topog, tmp_lon_f2D, tmp_topo_f2D, tmp_topo_test
        real, dimension(:), allocatable :: latf, lonf, topof
        logical :: lat_indices(size(lat,dim=1), size(lat,dim=2)), lon_indices(size(lon,dim=1), size(lon,dim=2))
        logical, dimension(:), allocatable :: tmp_lat, tmp_lon, tmp_indices_j
        logical, dimension(:,:), allocatable :: tmp_indices
        real, intent(in) :: clat, clon, width
        real, dimension(:,:,:), intent(inout) :: topo
        real, dimension(500000) :: tmp_lat_f, tmp_lon_f, tmp_topo_f
        real, dimension(:), allocatable :: tmp_lat_g, tmp_lon_g, tmp_topo_g
        integer :: nrecs, nlat, nlon
        integer :: i, j, k, stat

        integer(kind=int_size), dimension(:), allocatable :: lat_idx_sorted, lon_idx_sorted, iwork
        real, dimension(:), allocatable :: dummy, work

        ! get the outer-most axis of the lat-lon grid
        nlat = size(lat, dim=1)
        nlon = size(lon, dim=1)

        ! get the number of records
        ! nrecs = size(topo,dim=3)
        nrecs = 27

        allocate (obj%indices(nlon, nlat, nrecs), stat=stat)
        if (stat /= 0) then
            write(unit=error_unit, fmt='(A)') "Error allocating indices array for topo_obj"
            stop ALLOCATION_ERR
        end if

        ! we suppose that the allocation of the temporary arrays will not throw an error.
        allocate (tmp_lat(nlat))
        allocate (tmp_lon(nlon))
        allocate (tmp_indices(nlon, nlat))

        k = 1
        do i=26,nrecs
            print *, "nrec = ", i
            latg = spread(lat(:,i), 1, nlon)
            long = spread(lon(:,i), 2, nlat)
            topog = topo(:,:,i)

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

            tmp_indices = tmp_indices(:,ubound(tmp_indices,dim=2):lbound(tmp_indices,dim=2):-1)
            latg = latg(:,ubound(latg,dim=2):lbound(latg,dim=2):-1)
            tmp_indices_j = pack(tmp_indices, .true.)
            latf = pack(latg, .true.)
            lonf = pack(long, .true.)
            topof = pack(topog, .true.)
            do j=1,size(tmp_indices_j)
                if (tmp_indices_j(j) .eqv. .true.) then
                    tmp_lat_f(k) = latf(j)
                    tmp_lon_f(k) = lonf(j)
                    tmp_topo_f(k)= topof(j)
                    k = k + 1
                end if
            end do
        end do

        print *, k

        tmp_lat_g = tmp_lat_f(1:k-1)
        tmp_lon_g = tmp_lon_f(1:k-1)
        tmp_topo_g = tmp_topo_f(1:k-1)
    
        deallocate(tmp_lat)
        deallocate(tmp_lon)
        deallocate(tmp_indices)

        print *, "Before:"
        print *, minval(tmp_lat_g), maxval(tmp_lat_g)
        print *, minval(tmp_lon_g), maxval(tmp_lon_g)
        print *, minval(tmp_topo_g), maxval(tmp_topo_g)

        obj%lat  = pack(lat,  lat_indices)
        obj%lon  = pack(lon,  lon_indices)
        obj%topo = pack(topo, obj%indices)
        allocate (tmp_topo_test(size(obj%lon),size(obj%lat)))
        tmp_topo_test = reshape(obj%topo, shape=(/size(obj%lon),size(obj%lat)/), order=(/1,2/))
        tmp_topo_test = tmp_topo_test(:,ubound(tmp_topo_test,dim=2):lbound(tmp_topo_test,dim=2):-1)

        allocate (lat_idx_sorted(size(tmp_lat_g)))
        allocate (work(size(lat_idx_sorted)/2))
        allocate (iwork(size(lat_idx_sorted)/2))
        call sort_index(tmp_lat_g, lat_idx_sorted, work, iwork)

        tmp_topo_g = tmp_topo_g(lat_idx_sorted)
        tmp_lon_g = tmp_lon_g(lat_idx_sorted)
        deallocate(lat_idx_sorted)

        print *, "After:"
        print *, minval(tmp_lat_g), maxval(tmp_lat_g)
        print *, minval(tmp_lon_g), maxval(tmp_lon_g)
        print *, minval(tmp_topo_g), maxval(tmp_topo_g)

        tmp_lon_f2D = reshape(tmp_lon_g, shape=(/size(obj%lon),size(obj%lat)/))

        allocate (tmp_topo_f2D(size(obj%lon),size(obj%lat)))
        tmp_topo_f2D = reshape(tmp_topo_g, shape=(/size(obj%lon),size(obj%lat)/), order=(/1,2/))

        print *, "tmp_topo_f2D:"
        print *, minval(tmp_topo_f2D), maxval(tmp_topo_f2D)

        dummy = tmp_lon_f2D(1,:)

        ! call sort_index(dummy, lon_idx_sorted, work, iwork)
        ! do j = 1, size(tmp_topo_f2D), size(tmp_topo_f2D, dim=1)
        !     print *, j
        !     tmp_topo_f2D(j, :) = tmp_topo_f2D(j, lon_idx_sorted)
        ! end do


        ! obj%topo = pack(topo, obj%indices)

        ! print *, minval(obj%topo), maxval(obj%topo)

        ! print *, "obj%topo size: ", size(obj%topo)
        print *, "obj%lat size: ", size(obj%lat)
        print *, "obj%lon size: ", size(obj%lon)
        ! if ((size(obj%lat) * size(obj%lon)) /= size(obj%topo)) then
        !     write(unit=error_unit, fmt='(A)') "Error: Gathered subpoints shapes do not match."
        ! end if

    end subroutine set_attrs

    subroutine dealloc_obj(obj)
        implicit none
        type(topo_t), intent(inout) :: obj

        deallocate(obj%lat)
        deallocate(obj%lon)
        ! deallocate(obj%topo)
        ! deallocate(obj%indices)

    end subroutine dealloc_obj

end module topo_mod