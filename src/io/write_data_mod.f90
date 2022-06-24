module write_data_mod
    use :: netcdf
    use :: netcdf_check, only : nc_check
    use, intrinsic :: iso_fortran_env, only : error_unit
    use :: error_status, only : NOT_ALLOCATED_ERR

    implicit none

    ! private
    ! integer, parameter :: nlat, nlon

    ! public :: write_data
    
    interface write_data
        module procedure write_1D
        module procedure write_2D
        ! we do not need to write 3D datasets for now...
        ! module procedure write_3D 
    end interface write_data

contains
    function create_dataset(fname) result(ncid)
        implicit none
        character(len=*), intent(in) :: fname
        integer :: ncid

        call nc_check(nf90_create(fname, cmode=nf90_noclobber, ncid=ncid))

    end function create_dataset


    function open_dataset(fname) result(ncid)
        implicit none
        character(len=*), intent(in) :: fname
        integer :: ncid

        call nc_check(nf90_open(fname, mode=nf90_write, ncid=ncid))
        call nc_check(nf90_redef(ncid))

    end function open_dataset


    subroutine close_dataset(ncid)
        implicit none
        integer, intent(in) :: ncid

        call nc_check(nf90_close(ncid))

    end subroutine close_dataset


    function create_dim(ncid, dimname, dimsz) result(dimid)
        implicit none
        integer, intent(in) :: ncid, dimsz
        character(len=*), intent(in) :: dimname
        integer :: dimid

        ! define dimensions of datasets.
        call nc_check(nf90_def_dim(ncid, dimname, dimsz, dimid))

    end function create_dim


    function write_1D(ncid, varname, array, dimids) result(varid)
        implicit none

        character(len=*), intent(in) :: varname
        integer, intent(in) :: ncid
        integer, dimension(:), intent(in) :: dimids
        real, dimension(:), intent(in) :: array
        integer :: varid

        ! if .not. (allocated(nlat) and allocated(nlon)) then
        !     write(unit=error_unit, fmt=('A')) "Error: nlat and/or nlot not defined."
        !     stop NOT_ALLOCATED_ERR
        ! end if

        ! We store variables as single precision!
        call nc_check(nf90_def_var(ncid, varname, nf90_float, dimids, varid))
        call nc_check(nf90_enddef(ncid))
        call nc_check(nf90_put_var(ncid, varid, array))
        call nc_check(nf90_redef(ncid))

    end function write_1D

    function write_2D(ncid, varname, array, dimids) result(varid)
        implicit none

        character(len=*), intent(in) :: varname
        integer, intent(in) :: ncid
        integer, dimension(:), intent(in) :: dimids
        real, dimension(:,:), intent(in) :: array
        integer :: varid

        ! We store variables as single precision!
        call nc_check(nf90_def_var(ncid, varname, nf90_float, dimids, varid))
        call nc_check(nf90_enddef(ncid))
        call nc_check(nf90_put_var(ncid, varid, array))
        call nc_check(nf90_redef(ncid))

    end function write_2D

    subroutine write_attrs(ncid, varid, name, values)
        implicit none

        integer, intent(in) :: ncid, varid
        character(len=*), intent(in) :: name, values

        call nc_check(nf90_put_att(ncid, varid, name, values))

    end subroutine write_attrs

end module write_data_mod