# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake.exe

# The command to remove a file.
RM = /usr/bin/cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/build"

# Include any dependencies generated for this target.
include CMakeFiles/simulation.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/simulation.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/simulation.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/simulation.dir/flags.make

CMakeFiles/simulation.dir/simulation.c.o: CMakeFiles/simulation.dir/flags.make
CMakeFiles/simulation.dir/simulation.c.o: ../simulation.c
CMakeFiles/simulation.dir/simulation.c.o: CMakeFiles/simulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/simulation.dir/simulation.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/simulation.dir/simulation.c.o -MF CMakeFiles/simulation.dir/simulation.c.o.d -o CMakeFiles/simulation.dir/simulation.c.o -c "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/simulation.c"

CMakeFiles/simulation.dir/simulation.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/simulation.dir/simulation.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/simulation.c" > CMakeFiles/simulation.dir/simulation.c.i

CMakeFiles/simulation.dir/simulation.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/simulation.dir/simulation.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/simulation.c" -o CMakeFiles/simulation.dir/simulation.c.s

# Object files for target simulation
simulation_OBJECTS = \
"CMakeFiles/simulation.dir/simulation.c.o"

# External object files for target simulation
simulation_EXTERNAL_OBJECTS =

simulation.exe: CMakeFiles/simulation.dir/simulation.c.o
simulation.exe: CMakeFiles/simulation.dir/build.make
simulation.exe: CMakeFiles/simulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable simulation.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/simulation.dir/build: simulation.exe
.PHONY : CMakeFiles/simulation.dir/build

CMakeFiles/simulation.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/simulation.dir/cmake_clean.cmake
.PHONY : CMakeFiles/simulation.dir/clean

CMakeFiles/simulation.dir/depend:
	cd "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2" "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2" "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/build" "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/build" "/cygdrive/c/Users/Tric/Documents/C simulation/Chapter 2/build/CMakeFiles/simulation.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/simulation.dir/depend

