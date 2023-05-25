cmake_minimum_required(VERSION 3.16)
project(gmres)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "-Ofast -march=native -fsanitize=address -lblas")

find_package(BLAS REQUIRED)

add_executable(gmres main.f90)
target_link_libraries(gmres ${CMAKE_SOURCE_DIR}/libskit.a)
