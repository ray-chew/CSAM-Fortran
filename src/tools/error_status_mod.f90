module error_status
    implicit none
    integer, parameter ::   NETCDF_READ_ERR   = 1,      &
                            COMMAND_LINE_ERR  = 2,      &
                            ALLOCATION_ERR    = 3,      &
                            NOT_ALLOCATED_ERR = 4,      &    
                            LINALG_ERR = 5

end module error_status