cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(mesh)

# Build target module
add_executable(mesh main.cpp)
target_include_directories(mesh PRIVATE .)

# Link s0s
target_include_directories(mesh PRIVATE ../../include)
target_include_directories(mesh PRIVATE ../../../m0sh/include)
target_include_directories(mesh PRIVATE ../../../../thirdparty/eigen)
