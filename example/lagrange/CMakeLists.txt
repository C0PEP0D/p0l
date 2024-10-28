cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(lagrange)

# Build target module
add_executable(lagrange main.cpp)
target_include_directories(lagrange PRIVATE .)

# Link s0s
target_include_directories(lagrange PRIVATE ../../include)
target_include_directories(lagrange PRIVATE ../../../../thirdparty/eigen)
