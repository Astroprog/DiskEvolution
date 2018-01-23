mpicxx-openmpi-mp -O3 -std=c++11 -I/opt/local/include/openmpi-mp -L/opt/local/lib/openmpi-mp -o DiskEvolution -lmpi main.cpp ParameterParser.cpp Simulation.cpp GridGeometry.cpp DiskWind.cpp 
