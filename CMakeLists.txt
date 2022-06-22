cmake_minimum_required(VERSION 3.0)

project(four_topo LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_VERBOSE_MAKEFILE ON)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

add_compile_options(-g -Wall -Wextra)

add_executable(${PROJECT_NAME} 
                src/main.f90 
                src/io/read_grid_mod.f90
                src/io/netcdf_check_mod.f90
                src/tools/error_status_mod.f90
                src/tools/utils_mod.f90
		        src/var/topo_mod.f90
                )
target_include_directories(${PROJECT_NAME}
                            PUBLIC
                            "/usr/include/" 
                            )
target_link_libraries(${PROJECT_NAME} "netcdff")
