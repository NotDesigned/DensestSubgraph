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
include CMakeFiles/debug.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/debug.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/debug.dir/flags.make

CMakeFiles/debug.dir/test.cpp.o: CMakeFiles/debug.dir/flags.make
CMakeFiles/debug.dir/test.cpp.o: test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qingshuo/DensestSubgraph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/debug.dir/test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/debug.dir/test.cpp.o -c /home/qingshuo/DensestSubgraph/test.cpp

CMakeFiles/debug.dir/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/debug.dir/test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qingshuo/DensestSubgraph/test.cpp > CMakeFiles/debug.dir/test.cpp.i

CMakeFiles/debug.dir/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/debug.dir/test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qingshuo/DensestSubgraph/test.cpp -o CMakeFiles/debug.dir/test.cpp.s

# Object files for target debug
debug_OBJECTS = \
"CMakeFiles/debug.dir/test.cpp.o"

# External object files for target debug
debug_EXTERNAL_OBJECTS =

debug: CMakeFiles/debug.dir/test.cpp.o
debug: CMakeFiles/debug.dir/build.make
debug: utility/libutility.so
debug: reduction/libreduction.so
debug: allocation/liballocation.so
debug: extraction/libextraction.so
debug: verification/libverification.so
debug: CMakeFiles/debug.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qingshuo/DensestSubgraph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable debug"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/debug.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/debug.dir/build: debug

.PHONY : CMakeFiles/debug.dir/build

CMakeFiles/debug.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/debug.dir/cmake_clean.cmake
.PHONY : CMakeFiles/debug.dir/clean

CMakeFiles/debug.dir/depend:
	cd /home/qingshuo/DensestSubgraph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qingshuo/DensestSubgraph /home/qingshuo/DensestSubgraph /home/qingshuo/DensestSubgraph /home/qingshuo/DensestSubgraph /home/qingshuo/DensestSubgraph/CMakeFiles/debug.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/debug.dir/depend

