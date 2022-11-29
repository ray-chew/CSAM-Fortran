module topo_mod
    use, intrinsic :: iso_fortran_env, only : error_unit, DP => real64
    use :: error_status, only : ALLOCATION_ERR
    use :: triangle_mod, only : llgrid_t
    use omp_lib
    implicit none

    private

    public :: topo_t, get_topo, get_topo_idx, dealloc_topo_obj

    type :: topo_t
        real, dimension(:), allocatable :: lat, lon, lat_tri, lon_tri, topo_tri
        real, dimension(:,:), allocatable :: topo, lat_grid, lon_grid, topo_recon_2D

        integer :: nhar_i, nhar_j
        complex, dimension(:,:), allocatable :: fcoeffs
    end type topo_t

    interface get_topo
        module procedure get_topo_by_search
        module procedure get_topo_by_index
    end interface get_topo

contains

    subroutine get_topo_idx(lat, lon, llgrid, ncell, icon_topo_links)
        implicit none
        real, dimension(:,:), intent(in) :: lat
        real, dimension(:,:), intent(in) :: lon
        ! real, intent(inout) :: width
        type(llgrid_t), intent(out) :: llgrid
        integer, intent(in) :: ncell
        integer, dimension(:,:), intent(out) :: icon_topo_links

        integer :: nrecs, nlat, nlon
        integer :: i, j
        logical, dimension(:,:), allocatable :: tmp_indices, cond_lat, cond_lon
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

            ! tmp_indices = (abs(lat_grid - clat) <= width)
            cond_lat = ((llgrid%lat_min < lat_grid) .and. (lat_grid < llgrid%lat_max))
            cond_lon = ((llgrid%lon_min < lon_grid) .and. (lon_grid < llgrid%lon_max))
            tmp_indices = cond_lat .and. cond_lon
            ! tmp_indices = tmp_indices .and. (abs(lon_grid - clon) <= width)

            ! print *, count(tmp_indices)

            if (count(tmp_indices) > 0) then
                icon_topo_links(j,ncell) = i
                j = j + 1
            end if
        end do

        ! if (j == 1) then
        !     print *, "Doubling width for: thread num ", omp_get_thread_num(), "on grid cell", ncell, "with width", width
        !     width = 8.0 * width
        !     call get_topo_idx(lat, lon, clat, clon, width, ncell, icon_topo_links)
        ! end if

    end subroutine get_topo_idx

    subroutine get_topo_by_search(lat, lon, llgrid, topo, obj)
        implicit none
        type(topo_t), intent(out) :: obj
        real, dimension(:,:), intent(in) :: lat, lon
        real(kind=DP), dimension(:,:,:), intent(inout) :: topo
        type(llgrid_t), intent(out) :: llgrid

        integer :: nrecs, nlat, nlon
        integer :: i, j, stat
        logical, dimension(:), allocatable :: tmp_lat, tmp_lon
        logical, dimension(:,:), allocatable :: tmp_indices, cond_lat, cond_lon
        logical, dimension(:,:,:), allocatable :: tmp_indices_3D
        logical :: lat_indices(size(lat,dim=1), size(lat,dim=2)), lon_indices(size(lon,dim=1), size(lon,dim=2))
        real, dimension(:), allocatable :: tmp_topo
        real, dimension(:,:), allocatable :: lat_grid, lon_grid

        ! get the outer-most axis of the lat-lon grid
        nlat = size(lat, dim=1)
        nlon = size(lon, dim=1)

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

        do i=1,nrecs
        ! do i=23,26
            print *, "nrec = ", i
            lat_grid = spread(lat(:,i), 1, nlon)
            lon_grid = spread(lon(:,i), 2, nlat)
            ! topo_grid = topo(:,:,i)

            ! tmp_indices = (abs(lat_grid - clat) <= width)
            ! tmp_indices = tmp_indices .and. (abs(lon_grid - clon) <= width)
            cond_lat = ((llgrid%lat_min < lat_grid) .and. (lat_grid < llgrid%lat_max))
            cond_lon = ((llgrid%lon_min < lon_grid) .and. (lon_grid < llgrid%lon_max))
            tmp_indices = cond_lat .and. cond_lon

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

            ! We need to flip the latitude axis. Why? Eye-power.
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
    end subroutine get_topo_by_search


    subroutine get_topo_by_index(lat, lon, llgrid, topo, obj, link_map, ncell, tol)
        implicit none
        type(topo_t), intent(out) :: obj
        real, dimension(:,:), intent(in) :: lat, lon
        real(kind=DP), dimension(:,:,:), intent(in) :: topo
        integer, dimension(:), intent(in) :: link_map
        integer, intent(in) :: ncell
        type(llgrid_t), intent(out) :: llgrid

        integer :: nrecs, nlat, nlon, nrec
        integer :: i, j, stat
        logical, dimension(:), allocatable :: tmp_lat, tmp_lon
        logical, dimension(:,:), allocatable :: tmp_indices, cond_lat, cond_lon
        logical, dimension(:,:,:), allocatable :: tmp_indices_3D
        ! logical :: lat_indices(size(lat,dim=1)), lon_indices(size(lon,dim=1))
        ! logical :: lat_indices(size(lat,dim=1), size(lat,dim=2)), lon_indices(size(lon,dim=1), size(lon,dim=2))
        logical, dimension(:,:), allocatable :: lat_indices, lon_indices
        real, dimension(:), allocatable :: tmp_topo
        real, dimension(:,:), allocatable :: lat_grid, lon_grid, lat_nrecs, lon_nrecs, tmp_topo_grid
        real, dimension(:,:,:), allocatable :: topo_nrecs
        real, optional, intent(in) :: tol
        real :: tol_
        if (present(tol)) then
            tol_ = tol
        else
            tol_ = 1e-5
        end if
        
        ! get the outer-most axis of the lat-lon grid
        nlat = size(lat, dim=1)
        nlon = size(lon, dim=1)

        ! get the number of records
        ! nrecs = size(topo,dim=3)
        nrecs = count(link_map > 0)
        ! nrecs = 27

        ! print *, count(lat_indices)
        ! print *, count(lon_indices)

        ! we suppose that the allocation of the temporary 2D arrays will not throw an error.
        allocate (tmp_lat(nlat))
        allocate (tmp_lon(nlon))
        allocate (cond_lat(nlon, nlat))
        allocate (cond_lon(nlon, nlat))
        allocate (tmp_indices(nlon, nlat))
        allocate (lat_indices(nlat,nrecs))
        allocate (lon_indices(nlon,nrecs))
        allocate (lat_nrecs(nlat,nrecs))
        allocate (lon_nrecs(nlon,nrecs))

        lat_indices = .false.
        lon_indices = .false.

        allocate (tmp_indices_3D(nlon, nlat, nrecs), stat=stat)
        tmp_indices_3D = .false.
        if (stat /= 0) then
            write(unit=error_unit, fmt='(A)') "Error allocating indices array for topo_obj"
            stop ALLOCATION_ERR
        end if

        allocate (topo_nrecs(nlon,nlat,nrecs), stat=stat)
        if (stat /= 0) then
            write(unit=error_unit, fmt='(A)') "Error allocating nrecs array for topo_obj"
            stop ALLOCATION_ERR
        end if

        do i=1,nrecs
            if ((i == 1) .and. (link_map(i) == 0)) then
                write(unit=error_unit, fmt='(A)') "ICON grid cell ", ncell, "has no linked topography."
            end if

            if (link_map(i) == 0) exit

            nrec = link_map(i)
            
            lat_grid = spread(lat(:,nrec), 1, nlon)
            lon_grid = spread(lon(:,nrec), 2, nlat)
            ! topo_grid = topo(:,:,nrec)

            lat_nrecs(:,i) = lat(:,nrec)
            lon_nrecs(:,i) = lon(:,nrec)
            tmp_topo_grid = topo(:,:,nrec)
            tmp_topo_grid = tmp_topo_grid(:,ubound(tmp_topo_grid,dim=2):lbound(tmp_topo_grid,dim=2):-1)
            topo_nrecs(:,:,i) = tmp_topo_grid !topo(:,:,nrec)
            ! topo_nrecs(:,:,i) = topo(:,:,nrec)

            cond_lat = .false.
            cond_lon = .false.
            tmp_indices = .false.

            ! tmp_indices = (abs(lat_grid - clat) <= width)
            ! tmp_indices = tmp_indices .and. (abs(lon_grid - clon) <= width)
            cond_lat = ((llgrid%lat_min <= lat_grid) .and. (lat_grid <= llgrid%lat_max))
            cond_lon = ((llgrid%lon_min <= lon_grid) .and. (lon_grid <= llgrid%lon_max))
            tmp_indices = cond_lat .and. cond_lon
            ! print *, count(tmp_indices)
            
            tmp_lat = .false.
            do j=1,nlon
                tmp_lat = tmp_lat .or. tmp_indices(j,:)
            end do
            lat_indices(:,i) = tmp_lat

            ! print *, "count(tmp_lat) = ", count(tmp_lat)

            tmp_lon = .false.
            do j=1,nlat
                tmp_lon = tmp_lon .or. tmp_indices(:,j)
            end do
            lon_indices(:,i) = tmp_lon

            ! print *, "count(tmp_lon) = ", count(tmp_lon)

            ! Why do we need to flip the latitude axis?
            tmp_indices_3D(:,:,i) = tmp_indices !tmp_indices(:,ubound(tmp_indices,dim=2):lbound(tmp_indices,dim=2):-1)
        end do

        ! print *, "nloop done"

        ! tmp_lon = .false.
        ! do j = 1, nrecs
        !     tmp_lon = tmp_lon .or. lon_indices(:,j)
        ! end do

        ! print *, "count(tmp_lon) = ", count(tmp_lon)

        ! this do loop takes care of the corner cases where topography data are in separate nfiles/nrecs but with the same lat or lon underlying grid.
        ! put tolerance into nml flags
        ! if (nrecs > 2) then
        !     do i = 1, nrecs-1
        !         if (all(abs(lat_nrecs(:,i) - lat_nrecs(:,i+1)) < tol)) then
        !             lat_indices(:,i+1) = .false.
        !         end if
        !         if (all(abs(lon_nrecs(:,i) - lon_nrecs(:,i+1)) < tol)) then
        !             lon_indices(:,i+1) = .false.
        !         end if
        !     end do
        !     if (all(abs(lat_nrecs(:,1) - lat_nrecs(:,nrecs)) < tol)) then
        !         lat_indices(:,nrecs) = .false.
        !     end if
        !     if (all(abs(lon_nrecs(:,1) - lon_nrecs(:,nrecs)) < tol)) then
        !         lon_indices(:,nrecs) = .false.
        !     end if
        ! else if (nrecs == 2 ) then
        !     if (all(abs(lat_nrecs(:,1) - lat_nrecs(:,nrecs)) < tol)) then
        !         lat_indices(:,nrecs) = .false.
        !     end if
        !     if (all(abs(lon_nrecs(:,1) - lon_nrecs(:,nrecs)) < tol)) then
        !         lon_indices(:,nrecs) = .false.
        !     end if
        ! end if

        if (nrecs > 1) then
            do i = 0, nrecs-2
                call check_dupl_arrs(nrecs-i, lat_nrecs, tol, lat_indices)
                call check_dupl_arrs(nrecs-i, lon_nrecs, tol, lon_indices)
            end do
        end if

        ! print *, count(lat_indices)
        ! print *, count(lon_indices)
        ! print *, count(tmp_indices_3D)

        obj%lat  = pack(lat_nrecs,  lat_indices)
        obj%lon  = pack(lon_nrecs,  lon_indices)

        ! print *, "obj%lat size: ", size(obj%lat)
        ! print *, "obj%lon size: ", size(obj%lon)

        allocate (tmp_topo(size(obj%lat) * size(obj%lon)))
        
        tmp_topo = pack(topo_nrecs, tmp_indices_3D)
        
        obj%topo = reshape(tmp_topo, shape=(/size(obj%lon),size(obj%lat)/), order=(/1,2/))
        ! obj%topo = obj%topo(:,ubound(obj%topo,dim=2):lbound(obj%topo,dim=2):-1)

        if ((size(obj%lat) * size(obj%lon)) /= size(obj%topo)) then
            write(unit=error_unit, fmt='(A)') "Error: Gathered subpoints shapes do not match."
        end if

        obj%lat_grid = spread(obj%lat, 1, size(obj%lon))
        obj%lon_grid = spread(obj%lon, 2, size(obj%lat))
    end subroutine get_topo_by_index


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

        deallocate(obj%fcoeffs)

    end subroutine dealloc_topo_obj 


    subroutine check_dupl_arrs(end, arr_nrecs, tol, arr_indices)
        implicit none
        integer, intent(in) :: end
        ! integer, intent(in) :: nrecs
        real, dimension(:,:), intent(in) :: arr_nrecs
        real, intent(in) :: tol
        logical, dimension(:,:), intent(inout) :: arr_indices

        integer :: i

        do i = 1, end-1
            if (all(abs(arr_nrecs(:,end) - arr_nrecs(:,end-i)) < tol)) then
                arr_indices(:,end-i) = .false.
            end if
        end do

        ! if ((start == 1) .and. (nrecs > 2)) then
        !     if (all(abs(arr_nrecs(:,start) - arr_nrecs(:,nrecs)) < tol)) then
        !         arr_indices(:,start) = .false.
        !     end if           
        ! end if
        
    end subroutine check_dupl_arrs

end module topo_mod
