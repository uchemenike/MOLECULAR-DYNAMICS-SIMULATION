cmake_minimum_required(VERSION 3.0)
project(simulation C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

set(TARGET simulation)

# Add the source files to the executable
add_executable(${TARGET} simulation.c)
