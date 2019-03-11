cmake_minimum_required (VERSION 3.5.1)
project (tinyraytracer)

set(CMAKE_CXX_STANDARD 14)

file(GLOB SOURCES *.h *.cpp)
add_executable(tinyraytracer main.cpp geometry.h)
