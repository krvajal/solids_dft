cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

enable_language(Fortran)

project(freeparticle)



set(SRC_FILES freeparticle.f90)

include_directories(~/fortran-utils/src)
include_directories(/usr/local/include)
link_directories(~/fortran-utils/src)
link_directories(/usr/local/lib)
add_executable(free.x ${SRC_FILES})

target_link_libraries(free.x fortran_utils)
target_link_libraries(free.x lapack)

# test targets

add_executable(test.x projectors.f90 test_projectors.f90)
target_link_libraries(test.x fortran_utils)
target_link_libraries(test.x lapack)

# particle with gth
add_executable(gth.x particle_gth.f90 gth_potential.f90 projectors.f90 gvect.f90 xc.f90 density.f90)
target_link_libraries(gth.x fortran_utils)
target_link_libraries(gth.x lapack)

add_executable(testfft.x fft.f90 test_fft.f90 )
target_link_libraries(testfft.x fftw3)
target_link_libraries(testfft.x fortran_utils)
target_link_libraries(testfft.x lapack)




