module utils_mod
    use, intrinsic :: iso_fortran_env, only : error_unit, output_unit
    use :: error_status, only : COMMAND_LINE_ERR
    use :: const_mod, only : PI
    ! use :: stdlib_sorting, only : ord_sort
    implicit none

    interface get_fn
        module procedure get_nml_fn
        module procedure get_grid_fn
    end interface get_fn

    private
    public :: get_fn, get_namelist, rad_to_deg, get_N_unique, flags_t

    type :: flags_t
        logical :: debug
    end type flags_t

contains

    ! ref: https://cyber.dabamos.de/programming/modernfortran/namelists.html
    subroutine get_namelist(fn, fn_grid,fn_topo, fn_output, flags)
        !! Reads Namelist from given file.
        character(len=*),  intent(in)    :: fn
        character(len=*),  intent(out)   :: fn_grid
        character(len=*),  intent(out)   :: fn_topo
        character(len=*),  intent(out)   :: fn_output
        type(flags_t),     intent(out)   :: flags
        integer                          :: stat
        integer                          :: unit

        ! Namelist definition.
        namelist /USERDATA/ fn_grid, fn_topo, fn_output, flags

        ! Check whether file exists.
        inquire (file=fn, iostat=stat)

        if (stat /= 0) then
            write (error_unit, '("Error: input file ", a, " does not exist")') fn
            return
        end if

        ! Open and read Namelist file.
        open (action='read', file=fn, iostat=stat, newunit=unit)
        read (nml=USERDATA, iostat=stat, unit=unit)
        ! print *, "Stat = ", stat
        ! if (stat /= 0) write (error_unit, '("Error: invalid Namelist format")')

        close (unit)
        
        fn_grid = trim(fn_grid)
        fn_topo = trim(fn_topo)
        fn_output = trim(fn_output)

        write(unit=output_unit, nml=USERDATA)

    end subroutine get_namelist


    subroutine get_nml_fn(fn)
        implicit none
        character(len=*), intent(out) :: fn

        if (command_argument_count() /= 1) then
            write(unit=error_unit, fmt='(A)') "Argument error: Expected 1 input arguments."
            stop COMMAND_LINE_ERR
        end if

        call get_command_argument(1, fn)
        fn = trim(fn)

    end subroutine get_nml_fn

    subroutine get_grid_fn(fn_grid, fn_topo)
        implicit none
        character(len=*), intent(out) :: fn_grid
        character(len=*), intent(out) :: fn_topo

        if (command_argument_count() /= 2) then
            write(unit=error_unit, fmt='(A)') "Argument error: Expected 2 input arguments."
            stop COMMAND_LINE_ERR
        end if

        call get_command_argument(1, fn_grid)
        call get_command_argument(2, fn_topo)

    end subroutine get_grid_fn


    elemental subroutine rad_to_deg(value)
        implicit none
        real, intent(inout) :: value

        value = value * (180.0 / PI)
        
    end subroutine rad_to_deg


    function get_N_unique(arr) result(cnt)
        implicit none
        integer, dimension(:), intent(in) :: arr
        integer :: old_val = 0, cnt, i = 1
        ! integer, dimension(:), allocatable :: work

        ! the first value in the array is always unique.
        old_val = arr(1)
        cnt = 1

        do i=2,size(arr)
            if (arr(i) /= old_val) then
                cnt = cnt + 1
                old_val = arr(i)
            end if
        end do
    end function get_N_unique

end module utils_mod