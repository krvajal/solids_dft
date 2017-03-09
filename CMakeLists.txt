cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

enable_language(Fortran)

project(freeparticle)



<<<<<<< HEAD
=======



>>>>>>> 49974c06db85d3f26820264cf0cfc47e90f0ea14
include_directories(~/fortran-utils/src)
include_directories(/usr/local/include)
link_directories(~/fortran-utils/src)

add_definitions(-Wno-implicit-interface)
<<<<<<< HEAD
add_definitions(-fcheck=all -fbacktrace)
add_definitions(-ffpe-trap=zero,overflow,underflow)

=======
add_definitions(-fcheck=all)
>>>>>>> 49974c06db85d3f26820264cf0cfc47e90f0ea14

set(SRC_FILES freeparticle.f90)

link_directories(/usr/local/lib)





<<<<<<< HEAD
add_executable(gth.x main.f90 gth_potential.f90  gvect.f90 fft.f90 density.f90 xc.f90 energy.f90)
=======
add_executable(gth.x main.f90 gth_potential.f90 projectors.f90 gvect.f90 fft.f90 density.f90 xc.f90)
>>>>>>> 49974c06db85d3f26820264cf0cfc47e90f0ea14

target_link_libraries(gth.x fortran_utils)
target_link_libraries(gth.x lapack)
target_link_libraries(gth.x fftw3)




