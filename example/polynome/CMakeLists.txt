cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(polyome)

# Build target module
add_executable(polyome main.cpp)
target_include_directories(polyome PRIVATE .)

# Link
target_include_directories(polyome PRIVATE ../../include)
