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
include utility/CMakeFiles/utility.dir/depend.make

# Include the progress variables for this target.
include utility/CMakeFiles/utility.dir/progress.make

# Include the compile flags for this target's objects.
include utility/CMakeFiles/utility.dir/flags.make

utility/CMakeFiles/utility.dir/graph.cpp.o: utility/CMakeFiles/utility.dir/flags.make
utility/CMakeFiles/utility.dir/graph.cpp.o: utility/graph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qingshuo/DensestSubgraph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object utility/CMakeFiles/utility.dir/graph.cpp.o"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utility.dir/graph.cpp.o -c /home/qingshuo/DensestSubgraph/utility/graph.cpp

utility/CMakeFiles/utility.dir/graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utility.dir/graph.cpp.i"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qingshuo/DensestSubgraph/utility/graph.cpp > CMakeFiles/utility.dir/graph.cpp.i

utility/CMakeFiles/utility.dir/graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utility.dir/graph.cpp.s"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qingshuo/DensestSubgraph/utility/graph.cpp -o CMakeFiles/utility.dir/graph.cpp.s

utility/CMakeFiles/utility.dir/flownetwork.cpp.o: utility/CMakeFiles/utility.dir/flags.make
utility/CMakeFiles/utility.dir/flownetwork.cpp.o: utility/flownetwork.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qingshuo/DensestSubgraph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object utility/CMakeFiles/utility.dir/flownetwork.cpp.o"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utility.dir/flownetwork.cpp.o -c /home/qingshuo/DensestSubgraph/utility/flownetwork.cpp

utility/CMakeFiles/utility.dir/flownetwork.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utility.dir/flownetwork.cpp.i"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qingshuo/DensestSubgraph/utility/flownetwork.cpp > CMakeFiles/utility.dir/flownetwork.cpp.i

utility/CMakeFiles/utility.dir/flownetwork.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utility.dir/flownetwork.cpp.s"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qingshuo/DensestSubgraph/utility/flownetwork.cpp -o CMakeFiles/utility.dir/flownetwork.cpp.s

utility/CMakeFiles/utility.dir/args.cpp.o: utility/CMakeFiles/utility.dir/flags.make
utility/CMakeFiles/utility.dir/args.cpp.o: utility/args.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qingshuo/DensestSubgraph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object utility/CMakeFiles/utility.dir/args.cpp.o"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utility.dir/args.cpp.o -c /home/qingshuo/DensestSubgraph/utility/args.cpp

utility/CMakeFiles/utility.dir/args.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utility.dir/args.cpp.i"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qingshuo/DensestSubgraph/utility/args.cpp > CMakeFiles/utility.dir/args.cpp.i

utility/CMakeFiles/utility.dir/args.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utility.dir/args.cpp.s"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qingshuo/DensestSubgraph/utility/args.cpp -o CMakeFiles/utility.dir/args.cpp.s

utility/CMakeFiles/utility.dir/ratioselection.cpp.o: utility/CMakeFiles/utility.dir/flags.make
utility/CMakeFiles/utility.dir/ratioselection.cpp.o: utility/ratioselection.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qingshuo/DensestSubgraph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object utility/CMakeFiles/utility.dir/ratioselection.cpp.o"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/utility.dir/ratioselection.cpp.o -c /home/qingshuo/DensestSubgraph/utility/ratioselection.cpp

utility/CMakeFiles/utility.dir/ratioselection.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/utility.dir/ratioselection.cpp.i"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qingshuo/DensestSubgraph/utility/ratioselection.cpp > CMakeFiles/utility.dir/ratioselection.cpp.i

utility/CMakeFiles/utility.dir/ratioselection.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/utility.dir/ratioselection.cpp.s"
	cd /home/qingshuo/DensestSubgraph/utility && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qingshuo/DensestSubgraph/utility/ratioselection.cpp -o CMakeFiles/utility.dir/ratioselection.cpp.s

# Object files for target utility
utility_OBJECTS = \
"CMakeFiles/utility.dir/graph.cpp.o" \
"CMakeFiles/utility.dir/flownetwork.cpp.o" \
"CMakeFiles/utility.dir/args.cpp.o" \
"CMakeFiles/utility.dir/ratioselection.cpp.o"

# External object files for target utility
utility_EXTERNAL_OBJECTS =

utility/libutility.so: utility/CMakeFiles/utility.dir/graph.cpp.o
utility/libutility.so: utility/CMakeFiles/utility.dir/flownetwork.cpp.o
utility/libutility.so: utility/CMakeFiles/utility.dir/args.cpp.o
utility/libutility.so: utility/CMakeFiles/utility.dir/ratioselection.cpp.o
utility/libutility.so: utility/CMakeFiles/utility.dir/build.make
utility/libutility.so: utility/CMakeFiles/utility.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qingshuo/DensestSubgraph/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX shared library libutility.so"
	cd /home/qingshuo/DensestSubgraph/utility && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/utility.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
utility/CMakeFiles/utility.dir/build: utility/libutility.so

.PHONY : utility/CMakeFiles/utility.dir/build

utility/CMakeFiles/utility.dir/clean:
	cd /home/qingshuo/DensestSubgraph/utility && $(CMAKE_COMMAND) -P CMakeFiles/utility.dir/cmake_clean.cmake
.PHONY : utility/CMakeFiles/utility.dir/clean

utility/CMakeFiles/utility.dir/depend:
	cd /home/qingshuo/DensestSubgraph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qingshuo/DensestSubgraph /home/qingshuo/DensestSubgraph/utility /home/qingshuo/DensestSubgraph /home/qingshuo/DensestSubgraph/utility /home/qingshuo/DensestSubgraph/utility/CMakeFiles/utility.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : utility/CMakeFiles/utility.dir/depend

