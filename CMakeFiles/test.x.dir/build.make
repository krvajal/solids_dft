# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.5.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.5.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle

# Include any dependencies generated for this target.
include CMakeFiles/test.x.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test.x.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test.x.dir/flags.make

CMakeFiles/test.x.dir/projectors.f90.o: CMakeFiles/test.x.dir/flags.make
CMakeFiles/test.x.dir/projectors.f90.o: projectors.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/test.x.dir/projectors.f90.o"
	/usr/local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/projectors.f90 -o CMakeFiles/test.x.dir/projectors.f90.o

CMakeFiles/test.x.dir/projectors.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/test.x.dir/projectors.f90.i"
	/usr/local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/projectors.f90 > CMakeFiles/test.x.dir/projectors.f90.i

CMakeFiles/test.x.dir/projectors.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/test.x.dir/projectors.f90.s"
	/usr/local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/projectors.f90 -o CMakeFiles/test.x.dir/projectors.f90.s

CMakeFiles/test.x.dir/projectors.f90.o.requires:

.PHONY : CMakeFiles/test.x.dir/projectors.f90.o.requires

CMakeFiles/test.x.dir/projectors.f90.o.provides: CMakeFiles/test.x.dir/projectors.f90.o.requires
	$(MAKE) -f CMakeFiles/test.x.dir/build.make CMakeFiles/test.x.dir/projectors.f90.o.provides.build
.PHONY : CMakeFiles/test.x.dir/projectors.f90.o.provides

CMakeFiles/test.x.dir/projectors.f90.o.provides.build: CMakeFiles/test.x.dir/projectors.f90.o


CMakeFiles/test.x.dir/test_projectors.f90.o: CMakeFiles/test.x.dir/flags.make
CMakeFiles/test.x.dir/test_projectors.f90.o: test_projectors.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/test.x.dir/test_projectors.f90.o"
	/usr/local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/test_projectors.f90 -o CMakeFiles/test.x.dir/test_projectors.f90.o

CMakeFiles/test.x.dir/test_projectors.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/test.x.dir/test_projectors.f90.i"
	/usr/local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/test_projectors.f90 > CMakeFiles/test.x.dir/test_projectors.f90.i

CMakeFiles/test.x.dir/test_projectors.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/test.x.dir/test_projectors.f90.s"
	/usr/local/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/test_projectors.f90 -o CMakeFiles/test.x.dir/test_projectors.f90.s

CMakeFiles/test.x.dir/test_projectors.f90.o.requires:

.PHONY : CMakeFiles/test.x.dir/test_projectors.f90.o.requires

CMakeFiles/test.x.dir/test_projectors.f90.o.provides: CMakeFiles/test.x.dir/test_projectors.f90.o.requires
	$(MAKE) -f CMakeFiles/test.x.dir/build.make CMakeFiles/test.x.dir/test_projectors.f90.o.provides.build
.PHONY : CMakeFiles/test.x.dir/test_projectors.f90.o.provides

CMakeFiles/test.x.dir/test_projectors.f90.o.provides.build: CMakeFiles/test.x.dir/test_projectors.f90.o


# Object files for target test.x
test_x_OBJECTS = \
"CMakeFiles/test.x.dir/projectors.f90.o" \
"CMakeFiles/test.x.dir/test_projectors.f90.o"

# External object files for target test.x
test_x_EXTERNAL_OBJECTS =

test.x: CMakeFiles/test.x.dir/projectors.f90.o
test.x: CMakeFiles/test.x.dir/test_projectors.f90.o
test.x: CMakeFiles/test.x.dir/build.make
test.x: CMakeFiles/test.x.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking Fortran executable test.x"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test.x.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test.x.dir/build: test.x

.PHONY : CMakeFiles/test.x.dir/build

CMakeFiles/test.x.dir/requires: CMakeFiles/test.x.dir/projectors.f90.o.requires
CMakeFiles/test.x.dir/requires: CMakeFiles/test.x.dir/test_projectors.f90.o.requires

.PHONY : CMakeFiles/test.x.dir/requires

CMakeFiles/test.x.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test.x.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test.x.dir/clean

CMakeFiles/test.x.dir/depend:
	cd /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle /Users/carvajal/Dropbox/Maestria/codigos/dft/freeparticle/CMakeFiles/test.x.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test.x.dir/depend

