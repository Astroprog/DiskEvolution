cmake_minimum_required(VERSION 3.5)
project(DiskEvolution)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse3 -mfpmath=sse -fargument-noalias-global -ffast-math -march=native -std=c++11")
link_directories(/opt/local/lib/openmpi-gcc5)
include_directories(/opt/local/include/openmpi-gcc5)
set(SOURCE_FILES main.cpp GridGeometry.cpp GridGeometry.h DiskWind.cpp DiskWind.h ParameterParser.cpp ParameterParser.h Simulation.cpp Simulation.h)
add_executable(DiskEvolution ${SOURCE_FILES} GridGeometry.cpp GridGeometry.h DiskWind.cpp DiskWind.h ParameterParser.cpp ParameterParser.h Simulation.cpp Simulation.h)
target_link_libraries(DiskEvolution mpi)