module write_data_mod
    use :: netcdf
    use :: netcdf_check, only : nc_check
    use, intrinsic :: iso_fortran_env, only : error_unit, DP => real64
    use :: error_status, only : NOT_ALLOCATED_ERR

    implicit none

    ! private
    ! integer, parameter :: nlat, nlon

    ! public :: write_data
    
    interface write_data
        module procedure write_1D
        module procedure write_2D
        module procedure write_2D_int
        module procedure write_3D
        module procedure write_3D_cmplx
    end interface write_data

contains
    function create_dataset(fname) result(ncid)
        implicit none
        character(len=*), intent(in) :: fname
        integer :: ncid

        call nc_check(nf90_create(fname, cmode=nf90_clobber, ncid=ncid))

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


    function get_dim_id(ncid, dimname) result(dimid)
        implicit none
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: dimname
        integer :: dimid

        call nc_check(nf90_inq_dimid(ncid, dimname, dimid))

    end function get_dim_id


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

    function write_2D_int(ncid, varname, array, dimids) result(varid)
        implicit none

        character(len=*), intent(in) :: varname
        integer, intent(in) :: ncid
        integer, dimension(:), intent(in) :: dimids
        integer, dimension(:,:), intent(in) :: array
        integer :: varid

        ! We store variables as single precision!
        call nc_check(nf90_def_var(ncid, varname, nf90_int, dimids, varid))
        call nc_check(nf90_enddef(ncid))
        call nc_check(nf90_put_var(ncid, varid, array))
        call nc_check(nf90_redef(ncid))

    end function write_2D_int

    function write_3D(ncid, varname, array, dimids) result(varid)
        implicit none

        character(len=*), intent(in) :: varname
        integer, intent(in) :: ncid
        integer, dimension(:), intent(in) :: dimids
        real(kind=DP), dimension(:,:,:), intent(in) :: array
        integer :: varid

        ! varname = trim(varname)
        call nc_check(nf90_def_var(ncid, varname, nf90_double, dimids, varid))
        call nc_check(nf90_enddef(ncid))
        call nc_check(nf90_put_var(ncid, varid, array))
        call nc_check(nf90_redef(ncid))

    end function write_3D

    function write_3D_cmplx(ncid, varname, array, dimids) result(stat)
        implicit none

        character(len=*), intent(in) :: varname
        character(len=128) :: r_varname, re_ = "re_"
        character(len=128) :: i_varname, im_ = "im_"
        integer, intent(in) :: ncid
        integer, dimension(:), intent(in) :: dimids
        complex, dimension(:,:,:), intent(in) :: array
        real(kind=DP), dimension(:,:,:), allocatable :: r_arr, i_arr
        integer :: stat, re_id, im_id
        integer :: d1, d2, d3

        d1 = size(array, dim = 1)
        d2 = size(array, dim = 2)
        d3 = size(array, dim = 3)

        allocate (r_arr(d1,d2,d3))
        allocate (i_arr(d1,d2,d3))

        r_arr = real(array)
        i_arr = real(aimag(array))

        r_varname = trim(re_) // trim(varname)
        i_varname = trim(im_) // trim(varname)

        re_id = write_3D(ncid, trim(r_varname), r_arr, dimids)
        call write_attrs(ncid, re_id, 'long_name', 're(fcoeffs)')

        im_id = write_3D(ncid, trim(i_varname), i_arr, dimids)
        call write_attrs(ncid, im_id, 'long_name', 'im(fcoeffs)')

        stat = 1

    end function write_3D_cmplx

    subroutine write_attrs(ncid, varid, name, values)
        implicit none

        integer, intent(in) :: ncid, varid
        character(len=*), intent(in) :: name, values

        call nc_check(nf90_put_att(ncid, varid, name, values))

    end subroutine write_attrs

end module write_data_mod