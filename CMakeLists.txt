cmake_minimum_required(VERSION 3.5)
project(DiskEvolution)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -msse3 -std=c++11")

set(SOURCE_FILES main.cpp GridGeometry.cpp GridGeometry.h DiskWind.cpp DiskWind.h ParameterParser.cpp ParameterParser.h)
add_executable(DiskEvolution ${SOURCE_FILES} GridGeometry.cpp GridGeometry.h DiskWind.cpp DiskWind.h ParameterParser.cpp ParameterParser.h)