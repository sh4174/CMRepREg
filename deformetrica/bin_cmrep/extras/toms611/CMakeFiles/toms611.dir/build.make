# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/cmrep

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep

# Include any dependencies generated for this target.
include extras/toms611/CMakeFiles/toms611.dir/depend.make

# Include the progress variables for this target.
include extras/toms611/CMakeFiles/toms611.dir/progress.make

# Include the compile flags for this target's objects.
include extras/toms611/CMakeFiles/toms611.dir/flags.make

extras/toms611/CMakeFiles/toms611.dir/toms611.c.o: extras/toms611/CMakeFiles/toms611.dir/flags.make
extras/toms611/CMakeFiles/toms611.dir/toms611.c.o: /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/cmrep/extras/toms611/toms611.c
	$(CMAKE_COMMAND) -E cmake_progress_report /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object extras/toms611/CMakeFiles/toms611.dir/toms611.c.o"
	cd /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep/extras/toms611 && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/toms611.dir/toms611.c.o   -c /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/cmrep/extras/toms611/toms611.c

extras/toms611/CMakeFiles/toms611.dir/toms611.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/toms611.dir/toms611.c.i"
	cd /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep/extras/toms611 && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/cmrep/extras/toms611/toms611.c > CMakeFiles/toms611.dir/toms611.c.i

extras/toms611/CMakeFiles/toms611.dir/toms611.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/toms611.dir/toms611.c.s"
	cd /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep/extras/toms611 && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/cmrep/extras/toms611/toms611.c -o CMakeFiles/toms611.dir/toms611.c.s

extras/toms611/CMakeFiles/toms611.dir/toms611.c.o.requires:
.PHONY : extras/toms611/CMakeFiles/toms611.dir/toms611.c.o.requires

extras/toms611/CMakeFiles/toms611.dir/toms611.c.o.provides: extras/toms611/CMakeFiles/toms611.dir/toms611.c.o.requires
	$(MAKE) -f extras/toms611/CMakeFiles/toms611.dir/build.make extras/toms611/CMakeFiles/toms611.dir/toms611.c.o.provides.build
.PHONY : extras/toms611/CMakeFiles/toms611.dir/toms611.c.o.provides

extras/toms611/CMakeFiles/toms611.dir/toms611.c.o.provides.build: extras/toms611/CMakeFiles/toms611.dir/toms611.c.o

# Object files for target toms611
toms611_OBJECTS = \
"CMakeFiles/toms611.dir/toms611.c.o"

# External object files for target toms611
toms611_EXTERNAL_OBJECTS =

extras/toms611/libtoms611.a: extras/toms611/CMakeFiles/toms611.dir/toms611.c.o
extras/toms611/libtoms611.a: extras/toms611/CMakeFiles/toms611.dir/build.make
extras/toms611/libtoms611.a: extras/toms611/CMakeFiles/toms611.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libtoms611.a"
	cd /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep/extras/toms611 && $(CMAKE_COMMAND) -P CMakeFiles/toms611.dir/cmake_clean_target.cmake
	cd /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep/extras/toms611 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/toms611.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
extras/toms611/CMakeFiles/toms611.dir/build: extras/toms611/libtoms611.a
.PHONY : extras/toms611/CMakeFiles/toms611.dir/build

extras/toms611/CMakeFiles/toms611.dir/requires: extras/toms611/CMakeFiles/toms611.dir/toms611.c.o.requires
.PHONY : extras/toms611/CMakeFiles/toms611.dir/requires

extras/toms611/CMakeFiles/toms611.dir/clean:
	cd /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep/extras/toms611 && $(CMAKE_COMMAND) -P CMakeFiles/toms611.dir/cmake_clean.cmake
.PHONY : extras/toms611/CMakeFiles/toms611.dir/clean

extras/toms611/CMakeFiles/toms611.dir/depend:
	cd /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/cmrep /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/cmrep/extras/toms611 /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep/extras/toms611 /media/shong/IntHard1/Projects/4DShapeAnalysis/Code/SPTSkeleton/CMRepReg/deformetrica/bin_cmrep/extras/toms611/CMakeFiles/toms611.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : extras/toms611/CMakeFiles/toms611.dir/depend
