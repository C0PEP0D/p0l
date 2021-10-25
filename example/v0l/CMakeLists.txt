cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(v0l)

# Build target module
add_executable(v0l main.cpp)
target_include_directories(v0l PRIVATE .)

# Link s0s
target_include_directories(v0l PRIVATE ../../include)
target_include_directories(v0l PRIVATE ../../../m0sh/include)
target_include_directories(v0l PRIVATE ../../../v0l/include)
target_include_directories(v0l PRIVATE ../../../../thirdparty/eigen)
