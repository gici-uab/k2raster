# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_SOURCE_DIR = /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5

# Include any dependencies generated for this target.
include test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/flags.make

test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o: test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/flags.make
test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o: test/mvrtree/MVRTreeLoad.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o"
	cd /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o -c /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test/mvrtree/MVRTreeLoad.cc

test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.i"
	cd /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test/mvrtree/MVRTreeLoad.cc > CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.i

test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.s"
	cd /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test/mvrtree/MVRTreeLoad.cc -o CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.s

test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o.requires:
.PHONY : test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o.requires

test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o.provides: test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o.requires
	$(MAKE) -f test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/build.make test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o.provides.build
.PHONY : test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o.provides

test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o.provides.build: test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o

# Object files for target test-mvrtree-MVRTreeLoad
test__mvrtree__MVRTreeLoad_OBJECTS = \
"CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o"

# External object files for target test-mvrtree-MVRTreeLoad
test__mvrtree__MVRTreeLoad_EXTERNAL_OBJECTS =

bin/test-mvrtree-MVRTreeLoad: test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o
bin/test-mvrtree-MVRTreeLoad: test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/build.make
bin/test-mvrtree-MVRTreeLoad: bin/libspatialindex.a
bin/test-mvrtree-MVRTreeLoad: test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../bin/test-mvrtree-MVRTreeLoad"
	cd /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test-mvrtree-MVRTreeLoad.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/build: bin/test-mvrtree-MVRTreeLoad
.PHONY : test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/build

test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/requires: test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/mvrtree/MVRTreeLoad.cc.o.requires
.PHONY : test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/requires

test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/clean:
	cd /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test && $(CMAKE_COMMAND) -P CMakeFiles/test-mvrtree-MVRTreeLoad.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/clean

test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/depend:
	cd /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5 /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5 /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test /home/fernando/Doutoramento/Software/RTree/spatialindex-src-1.8.5/test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/test-mvrtree-MVRTreeLoad.dir/depend
