CC = mpicxx-openmpi-mp
LDFLAGS = -O3 -std=c++11
CCLAGS = -O3 -std=c++11 

# -I/opt/local/include/openmpi-mp
# -L/opt/local/lib/openmpi-mp

CCFLAGS += `pkg-config --cflags ompi`
LDFLAGS += `pkg-config --libs ompi`
LDFLAGS += -L/opt/local/lib -I/opt/local/include 

EXECUTABLE = bin/DiskEvolution
SOURCES = $(wildcard */*.cpp)
OBJECTS = $(SOURCES: .cpp=.o)


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(SOURCE_DIR)%.o: $SOURCE_DIR%.cpp
	$(CC) $(CCFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE)
