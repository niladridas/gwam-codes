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
include CMakeFiles/collectData3.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/collectData3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/collectData3.dir/flags.make

CMakeFiles/collectData3.dir/src/collectData3.cpp.o: CMakeFiles/collectData3.dir/flags.make
CMakeFiles/collectData3.dir/src/collectData3.cpp.o: ../src/collectData3.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/robot/Mycodes/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/collectData3.dir/src/collectData3.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/collectData3.dir/src/collectData3.cpp.o -c /home/robot/Mycodes/src/collectData3.cpp

CMakeFiles/collectData3.dir/src/collectData3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/collectData3.dir/src/collectData3.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/robot/Mycodes/src/collectData3.cpp > CMakeFiles/collectData3.dir/src/collectData3.cpp.i

CMakeFiles/collectData3.dir/src/collectData3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/collectData3.dir/src/collectData3.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/robot/Mycodes/src/collectData3.cpp -o CMakeFiles/collectData3.dir/src/collectData3.cpp.s

CMakeFiles/collectData3.dir/src/collectData3.cpp.o.requires:
.PHONY : CMakeFiles/collectData3.dir/src/collectData3.cpp.o.requires

CMakeFiles/collectData3.dir/src/collectData3.cpp.o.provides: CMakeFiles/collectData3.dir/src/collectData3.cpp.o.requires
	$(MAKE) -f CMakeFiles/collectData3.dir/build.make CMakeFiles/collectData3.dir/src/collectData3.cpp.o.provides.build
.PHONY : CMakeFiles/collectData3.dir/src/collectData3.cpp.o.provides

CMakeFiles/collectData3.dir/src/collectData3.cpp.o.provides.build: CMakeFiles/collectData3.dir/src/collectData3.cpp.o
.PHONY : CMakeFiles/collectData3.dir/src/collectData3.cpp.o.provides.build

# Object files for target collectData3
collectData3_OBJECTS = \
"CMakeFiles/collectData3.dir/src/collectData3.cpp.o"

# External object files for target collectData3
collectData3_EXTERNAL_OBJECTS =

collectData3: CMakeFiles/collectData3.dir/src/collectData3.cpp.o
collectData3: /usr/local/lib/libboost_thread.so
collectData3: /usr/local/lib/libboost_python.so
collectData3: /usr/xenomai/lib/libnative.so
collectData3: /usr/xenomai/lib/libxenomai.so
collectData3: /usr/xenomai/lib/librtdm.so
collectData3: /usr/lib/libpython2.6.so
collectData3: libconstants.a
collectData3: libsamlibs.a
collectData3: CMakeFiles/collectData3.dir/build.make
collectData3: CMakeFiles/collectData3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable collectData3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/collectData3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/collectData3.dir/build: collectData3
.PHONY : CMakeFiles/collectData3.dir/build

CMakeFiles/collectData3.dir/requires: CMakeFiles/collectData3.dir/src/collectData3.cpp.o.requires
.PHONY : CMakeFiles/collectData3.dir/requires

CMakeFiles/collectData3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/collectData3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/collectData3.dir/clean

CMakeFiles/collectData3.dir/depend:
	cd /home/robot/Mycodes/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/robot/Mycodes /home/robot/Mycodes /home/robot/Mycodes/build /home/robot/Mycodes/build /home/robot/Mycodes/build/CMakeFiles/collectData3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/collectData3.dir/depend

