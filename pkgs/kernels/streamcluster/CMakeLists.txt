cmake_minimum_required(VERSION 3.16)

# set the project name
project(project)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES src/*.cpp src/*.hpp)

# add the executable
add_executable(project ${SOURCES})