cmake_minimum_required(VERSION 3.0)

project(four_topo LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(fortran_stdlib REQUIRED)
find_package(OpenMP REQUIRED)

find_path (NETCDFF_INCLUDES_F90 netcdf.mod HINTS NETCDFF_DIR NETCDFF_DIR)
message(" inc ff =  ${NETCDFF_INCLUDES_F90}  ")
find_library (NETCDFF_LIBRARIES_F90  NAMES netcdff HINTS NETCDFF_DIR )
message(" lib ff = ${NETCDFF_LIBRARIES_F90}  ")
mark_as_advanced(NETCDF_LIBRARIES_F90)

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_VERBOSE_MAKEFILE ON)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_compile_options(-fopenmp -g -Wall -Wextra -fdefault-real-8)

link_libraries("-fopenmp")

add_executable(${PROJECT_NAME} 
                src/orog_source.f90
                src/io/read_data_mod.f90
                src/io/write_data_mod.f90
                src/io/netcdf_check_mod.f90
                src/tools/const_mod.f90
                src/tools/error_status_mod.f90
                src/tools/utils_mod.f90
                src/tools/triangle_mod.f90
                src/var/topo_mod.f90
                src/math/fourier_mod.f90
                src/math/lin_reg_mod.f90
                )

target_include_directories(${PROJECT_NAME}
                            PUBLIC
                            "/usr/include"
                            ${NETCDFF_INCLUDES_F90}
                            )

target_link_libraries(  ${PROJECT_NAME} 
                        PRIVATE
                        ${NETCDFF_LIBRARIES_F90}
                        fortran_stdlib::fortran_stdlib
                        blas
                        lapack
                        )
