# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/qingshuo/DensestSubgraph

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qingshuo/DensestSubgraph

# Include any dependencies generated for this target.
include allocation/CMakeFiles/allocation.dir/depend.make

# Include the progress variables for this target.
include allocation/CMakeFiles/allocation.dir/progress.make

# Include the compile flags for this target's objects.
include allocation/CMakeFiles/allocation.dir/flags.make

allocation/CMakeFiles/allocation.dir/allocation.cpp.o: allocation/CMakeFiles/allocation.dir/flags.make
allocation/CMakeFiles/allocation.dir/allocation.cpp.o: allocation/allocation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qingshuo/DensestSubgraph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object allocation/CMakeFiles/allocation.dir/allocation.cpp.o"
	cd /home/qingshuo/DensestSubgraph/allocation && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/allocation.dir/allocation.cpp.o -c /home/qingshuo/DensestSubgraph/allocation/allocation.cpp

allocation/CMakeFiles/allocation.dir/allocation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/allocation.dir/allocation.cpp.i"
	cd /home/qingshuo/DensestSubgraph/allocation && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qingshuo/DensestSubgraph/allocation/allocation.cpp > CMakeFiles/allocation.dir/allocation.cpp.i

allocation/CMakeFiles/allocation.dir/allocation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/allocation.dir/allocation.cpp.s"
	cd /home/qingshuo/DensestSubgraph/allocation && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qingshuo/DensestSubgraph/allocation/allocation.cpp -o CMakeFiles/allocation.dir/allocation.cpp.s

# Object files for target allocation
allocation_OBJECTS = \
"CMakeFiles/allocation.dir/allocation.cpp.o"

# External object files for target allocation
allocation_EXTERNAL_OBJECTS =

allocation/liballocation.so: allocation/CMakeFiles/allocation.dir/allocation.cpp.o
allocation/liballocation.so: allocation/CMakeFiles/allocation.dir/build.make
allocation/liballocation.so: allocation/CMakeFiles/allocation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qingshuo/DensestSubgraph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library liballocation.so"
	cd /home/qingshuo/DensestSubgraph/allocation && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/allocation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
allocation/CMakeFiles/allocation.dir/build: allocation/liballocation.so

.PHONY : allocation/CMakeFiles/allocation.dir/build

allocation/CMakeFiles/allocation.dir/clean:
	cd /home/qingshuo/DensestSubgraph/allocation && $(CMAKE_COMMAND) -P CMakeFiles/allocation.dir/cmake_clean.cmake
.PHONY : allocation/CMakeFiles/allocation.dir/clean

allocation/CMakeFiles/allocation.dir/depend:
	cd /home/qingshuo/DensestSubgraph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qingshuo/DensestSubgraph /home/qingshuo/DensestSubgraph/allocation /home/qingshuo/DensestSubgraph /home/qingshuo/DensestSubgraph/allocation /home/qingshuo/DensestSubgraph/allocation/CMakeFiles/allocation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : allocation/CMakeFiles/allocation.dir/depend

