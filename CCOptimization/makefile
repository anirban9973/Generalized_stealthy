#=============================-
# Makefile to compile a simple file format manager
# for converting a ConfigPack file to txt-format files, and vice versa.
#
# Author: Jaeuk Kim
# Date: February 2022
# Email: phy000.kim@gmail.com
#=============================

#==============================
# Spack-based dependency resolution (machine-independent)
#==============================

#----- Source directories (relative paths, unchanged) -----
core = ../cores/
pair = ../PairStat/
pot  = ../Potential/

#----- Locate packages via spack -----
GSL_PREFIX   := $(shell spack location -i gsl)
NLOPT_PREFIX := $(shell spack location -i nlopt)
BOOST_PREFIX := $(shell spack location -i boost)

IDIR = -I$(GSL_PREFIX)/include \
       -I$(NLOPT_PREFIX)/include \
       -I$(BOOST_PREFIX)/include \
       -I$(core) -I$(pair) -I$(pot)

GSLLIB  = -L$(GSL_PREFIX)/lib -lgsl -lgslcblas -Wl,-rpath,$(GSL_PREFIX)/lib
NLOPTLIB = -L$(NLOPT_PREFIX)/lib -lnlopt -Wl,-rpath,$(NLOPT_PREFIX)/lib
LIBS = -lm $(GSLLIB) $(NLOPTLIB)

PREPRO =

#----- Output directories and target -----
ODIR   = ./obj/
EXEdir = ./
TARGET = $(EXEdir)CCO.out

#==============================
#	Common flages
#==============================
CC = g++
DEBUG = -g
CFLAGS = -Wall -Wno-reorder -std=c++14 $(DEBUG) $(IDIR) -c -fopenmp -O3
LFLAGS = -Wall -Wno-reorder -std=c++14 $(DEBUG) $(LIBS) -fopenmp -O3

$(info	$(CFLAGS))
$(info	$(LFLAGS))

objs_core = RandomGenerator.o etc.o GeometryVector.o PeriodicCellList.o 
objs_pair = StructureFactor.o PairCorrelation.o
objs_pot = Potential.o StructureOptimization.o MD_System.o
objs_main = main.o MultiStealthyPotential.o
#objs_MC = MonteCarlo.o

DEPS = $(core)GeometryVector.h $(core)RandomGenerator.h $(core)etc.h $(core)PeriodicCellList.h \
$(pair)StructureFactor.h $(pair)PairCorrelation.h \
$(pot)Potential.h $(pot)CollectiveCoordinatePotential.h $(pot)StructureOptimization.h $(pot)MD_System.h \
RepulsiveCCPotential.h MultiStealthyPotential.h

#  ../GenerateConfigs.h ../HardSphereFluids.h ../RandomSequentialAddition.h

OBJ_core = $(patsubst %.o, $(ODIR)%.o, $(objs_core))
OBJ_pot = $(patsubst %.o, $(ODIR)%.o, $(objs_pot))
OBJ_pair = $(patsubst %.o, $(ODIR)%.o, $(objs_pair))
OBJ_main = $(patsubst %.o, $(ODIR)%.o, $(objs_main))
OBJ = $(OBJ_core) $(OBJ_pot) $(OBJ_main) $(OBJ_pair)

$(info	$(OBJ))

SRC_core = $(patsubst %.o, $(core)%.cpp, $(objs_core))
SRC_pot = $(patsubst %.o, $(pot)%.cpp, $(objs_pot))
SRC_pair = $(patsubst %.o, $(pair)%.cpp, $(objs_pair))
SRC_main = $(patsubst %.o, %.cpp, $(objs_main))
#SRC_misc = $(patsubst %.o, ../%.cpp, $(objs_misc))

$(info	$(SRC_core))
$(info	$(SRC_pot))
$(info	$(SRC_pair))
$(info	$(SRC_main))



all: | $(ODIR)
all: $(OBJ) $(TARGET)

$(ODIR):
	mkdir -p $@

# compilations...

$(OBJ_pair): $(ODIR)%.o: $(addprefix $(pair), %.cpp) $(DEPS) #$(addprefix $(pair), %.cpp) $(DEPS)
	echo $@ $<  
	$(CC) $(PREPRO) -o $@ $< $(CFLAGS)  

$(OBJ_core): $(ODIR)%.o: $(addprefix $(core), %.cpp) $(DEPS)#$(addprefix ../, %.cpp) $(DEPS)
	echo $@ $<  
	$(CC) $(PREPRO) -o $@ $< $(CFLAGS)  

$(OBJ_pot): $(ODIR)%.o: $(addprefix $(pot), %.cpp) $(DEPS) #$(addprefix ../, %.cpp) $(DEPS)
	echo $@ $<  
	$(CC) $(PREPRO) -o $@ $< $(CFLAGS)  

$(OBJ_main): $(ODIR)%.o: $(notdir %.cpp) $(DEPS) #$(notdir %.cpp) $(DEPS)
	echo $@ $<
	$(CC) $(PREPRO) -o $@ $< $(CFLAGS) 

# $(OBJ): $(ODIR)%.o: $(SRC) $(DEPS)
# 	echo $@ $<
# 	$(CC) -o $@ $< $(CFLAGS) 

$(TARGET): $(OBJ)
	@echo "Start linking"
	$(CC) $(PREPRO) -o $@ $(OBJ) $(LFLAGS)


clean:
	rm -f $(ODIR)*.o $(TARGET)

.PHONY: all
all:
