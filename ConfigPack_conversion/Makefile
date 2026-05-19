CXX      = g++
CXXFLAGS = -std=c++17 -O3 -Wall -fopenmp

# Spack paths (fall back to /usr/local if spack is unavailable)
GSL_ROOT   = $(shell spack location -i gsl   2>/dev/null || echo /usr/local)
BOOST_ROOT = $(shell spack location -i boost 2>/dev/null || echo /usr/local)
HDF5_ROOT      = $(shell spack location -i hdf5     2>/dev/null || echo /usr/local)
HIGHFIVE_ROOT  = $(shell spack location -i highfive 2>/dev/null || echo /usr/local)
MPI_ROOT   = $(shell spack location -i openmpi 2>/dev/null || echo /usr/local)

# Core library includes (ConfigPackData and its dependencies)
CORE_INCLUDE = ../cores

INCLUDES = -I$(CORE_INCLUDE) -I. \
           -I$(BOOST_ROOT)/include \
           -I$(GSL_ROOT)/include \
           -I$(HDF5_ROOT)/include \
           -I$(HIGHFIVE_ROOT)/include \
           -I$(MPI_ROOT)/include

LDFLAGS = -L$(GSL_ROOT)/lib   -Wl,-rpath,$(GSL_ROOT)/lib \
          -L$(BOOST_ROOT)/lib -Wl,-rpath,$(BOOST_ROOT)/lib \
          -L$(HDF5_ROOT)/lib  -Wl,-rpath,$(HDF5_ROOT)/lib \
          -L$(MPI_ROOT)/lib   -Wl,-rpath,$(MPI_ROOT)/lib \
          -fopenmp

LIBS = -lgsl -lgslcblas -lboost_system -lhdf5 -lmpi

TARGET = configpack_to_h5
SOURCE = ConfigPack_To_H5.cpp

CORE_SOURCES = $(CORE_INCLUDE)/etc.cpp \
               $(CORE_INCLUDE)/GeometryVector.cpp \
               $(CORE_INCLUDE)/RandomGenerator.cpp \
               $(CORE_INCLUDE)/PeriodicCellList.cpp

OBJECTS = $(SOURCE:.cpp=.o) $(CORE_SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS)

.PHONY: all clean
