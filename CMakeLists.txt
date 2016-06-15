cmake_minimum_required(VERSION 3.5)
project(DiskEvolution)


if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_C_COMPILER /usr/lib64/openmpi/bin/mpicc)
    set(CMAKE_CXX_COMPILER /usr/lib64/openmpi/bin/mpicxx)
    link_directories(/usr/lib64/openmpi/lib)
    include_directories(/usr/include/openmpi-x86_64)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_C_COMPILER /opt/local/bin/mpicc-openmpi-gcc5)
    set(CMAKE_CXX_COMPILER /opt/local/bin/mpicxx-openmpi-gcc5)
    link_directories(/opt/local/lib/openmpi-gcc5)
    include_directories(/opt/local/include/openmpi-gcc5)
endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse3 -mfpmath=sse -fargument-noalias-global -ffast-math -march=native -std=c++11")

set(SOURCE_FILES main.cpp GridGeometry.cpp GridGeometry.h DiskWind.cpp DiskWind.h ParameterParser.cpp ParameterParser.h Simulation.cpp Simulation.h)
add_executable(DiskEvolution ${SOURCE_FILES} GridGeometry.cpp GridGeometry.h DiskWind.cpp DiskWind.h ParameterParser.cpp ParameterParser.h Simulation.cpp Simulation.h)
target_link_libraries(DiskEvolution mpi)