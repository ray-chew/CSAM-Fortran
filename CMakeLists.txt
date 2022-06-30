cmake_minimum_required(VERSION 3.0)

project(four_topo LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(fortran_stdlib REQUIRED)

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_VERBOSE_MAKEFILE ON)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

add_compile_options(-g -Wall -Wextra -fdefault-real-8)

add_executable(${PROJECT_NAME} 
                src/main.f90 
                src/io/read_data_mod.f90
                src/io/write_data_mod.f90
                src/io/netcdf_check_mod.f90
                src/tools/error_status_mod.f90
                src/tools/utils_mod.f90
	            src/var/topo_mod.f90
                src/math/fourier_mod.f90
                src/math/lin_reg_mod.f90
                )

target_include_directories(${PROJECT_NAME}
                            PUBLIC
                            "/usr/include/"
                            )
target_link_libraries(  ${PROJECT_NAME} 
                        PRIVATE
                        "netcdff"
                        fortran_stdlib::fortran_stdlib
                        blas
                        lapack
                        )
