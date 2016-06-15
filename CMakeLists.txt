cmake_minimum_required(VERSION 3.5)
project(DiskEvolution)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse3 -mfpmath=sse -fargument-noalias-global -ffast-math -march=native -std=c++11")

set(Boost_INCLUDE_DIR /opt/local/include)
set(Boost_LIBRARY_DIR /opt/local/lib)
find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

set(SOURCE_FILES main.cpp GridGeometry.cpp GridGeometry.h DiskWind.cpp DiskWind.h ParameterParser.cpp ParameterParser.h Simulation.cpp Simulation.h)
add_executable(DiskEvolution ${SOURCE_FILES})
target_link_libraries( DiskEvolution ${Boost_LIBRARIES} )
