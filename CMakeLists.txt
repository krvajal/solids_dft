cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

enable_language(Fortran)

project(freeparticle)
set(SRC_FILES freeparticle.f90)

include_directories(/home/miguel/fortran-utils/src)
link_directories(/home/miguel/fortran-utils/src)
add_executable(free.x ${SRC_FILES})

target_link_libraries(free.x fortran_utils)
target_link_libraries(free.x lapack)