module netcdf_check
    use netcdf
    use error_status, only : NETCDF_READ_ERR
    use, intrinsic :: iso_fortran_env, only : error_unit
    implicit none

    private

        public :: nc_check

contains

    subroutine nc_check(stat)
        implicit none
        integer, intent(in) :: stat

        if (stat /= nf90_noerr) then
                write(unit=error_unit, fmt='(2A)') 'NetCDF read error: ', trim(nf90_strerror(stat))
            stop NETCDF_READ_ERR
        end if

    end subroutine nc_check

end module netcdf_check