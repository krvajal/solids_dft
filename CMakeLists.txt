cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

enable_language(Fortran)

project(freeparticle)






include_directories(~/fortran-utils/src)
include_directories(/usr/local/include)
link_directories(~/fortran-utils/src)

add_definitions(-Wno-implicit-interface)

add_definitions(-fcheck=all -fbacktrace)
add_definitions(-ffpe-trap=zero,overflow,underflow)

set(SRC_FILES freeparticle.f90)

link_directories(/usr/local/lib)



add_executable(gth.x main.f90 pseudopot.f90 projectors.f90 gvect.f90 fft.f90 density.f90 xc.f90)


target_link_libraries(gth.x fortran_utils)
target_link_libraries(gth.x lapack)
target_link_libraries(gth.x fftw3)




