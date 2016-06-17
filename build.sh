module load mpi/openmpi/1.10-gnu-5.2
mpicxx -O3 -std=c++11 -I/opt/bwhpc/common/mpi/openmpi/1.10.2-gnu-5.2/include/ -L/opt/bwhpc/common/mpi/openmpi/1.10.2-gnu-5.2/lib -o DiskEvolution -lmpi main.cpp ParameterParser.cpp Simulation.cpp GridGeometry.cpp DiskWind.cpp 
