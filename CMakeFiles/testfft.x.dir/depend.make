# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

CMakeFiles/testfft.x.dir/fft.f90.o: /usr/local/include/fftw3.f03

CMakeFiles/testfft.x.dir/fft.mod.proxy: CMakeFiles/testfft.x.dir/fft.f90.o.provides
CMakeFiles/testfft.x.dir/fft.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod fft CMakeFiles/testfft.x.dir/fft.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/testfft.x.dir/fft.f90.o.provides.build
CMakeFiles/testfft.x.dir/build: CMakeFiles/testfft.x.dir/fft.f90.o.provides.build

CMakeFiles/testfft.x.dir/test_fft.f90.o.requires: CMakeFiles/testfft.x.dir/fft.mod.proxy
CMakeFiles/testfft.x.dir/test_fft.f90.o: CMakeFiles/testfft.x.dir/fft.mod.stamp
CMakeFiles/testfft.x.dir/test_fft.f90.o: /home/miguel/fortran-utils/src/types.mod