# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/robot/Mycodes

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/robot/Mycodes/build

# Include any dependencies generated for this target.
include CMakeFiles/test_Controller3.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_Controller3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_Controller3.dir/flags.make

CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o: CMakeFiles/test_Controller3.dir/flags.make
CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o: ../src/test_Controller3.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/robot/Mycodes/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o -c /home/robot/Mycodes/src/test_Controller3.cpp

CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/robot/Mycodes/src/test_Controller3.cpp > CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.i

CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/robot/Mycodes/src/test_Controller3.cpp -o CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.s

CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o.requires:
.PHONY : CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o.requires

CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o.provides: CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o.requires
	$(MAKE) -f CMakeFiles/test_Controller3.dir/build.make CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o.provides.build
.PHONY : CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o.provides

CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o.provides.build: CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o
.PHONY : CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o.provides.build

# Object files for target test_Controller3
test_Controller3_OBJECTS = \
"CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o"

# External object files for target test_Controller3
test_Controller3_EXTERNAL_OBJECTS =

test_Controller3: CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o
test_Controller3: /usr/local/lib/libboost_thread.so
test_Controller3: /usr/local/lib/libboost_python.so
test_Controller3: /usr/xenomai/lib/libnative.so
test_Controller3: /usr/xenomai/lib/libxenomai.so
test_Controller3: /usr/xenomai/lib/librtdm.so
test_Controller3: /usr/lib/libpython2.6.so
test_Controller3: libconstants.a
test_Controller3: CMakeFiles/test_Controller3.dir/build.make
test_Controller3: CMakeFiles/test_Controller3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable test_Controller3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_Controller3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_Controller3.dir/build: test_Controller3
.PHONY : CMakeFiles/test_Controller3.dir/build

CMakeFiles/test_Controller3.dir/requires: CMakeFiles/test_Controller3.dir/src/test_Controller3.cpp.o.requires
.PHONY : CMakeFiles/test_Controller3.dir/requires

CMakeFiles/test_Controller3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_Controller3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_Controller3.dir/clean

CMakeFiles/test_Controller3.dir/depend:
	cd /home/robot/Mycodes/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/robot/Mycodes /home/robot/Mycodes /home/robot/Mycodes/build /home/robot/Mycodes/build /home/robot/Mycodes/build/CMakeFiles/test_Controller3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_Controller3.dir/depend
