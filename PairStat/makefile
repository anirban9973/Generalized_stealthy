#=============================-
# Makefile to compile a simple file format manager
# for converting a ConfigPack file to txt-format files, and vice versa.
#
# Author: Jaeuk Kim
# Date: February 2022
# Email: phy000.kim@gmail.com
#=============================

MACHINE = HOME
#MACHINE = DELLA
#MACHINE = NURION
core = /home/jaeukk/new_codes/c++/cores/
#==============================
# DELLA
#==============================
ifeq ($(MACHINE), DELLA)
#echo $(MACHINE)
$(info	$(MACHINE))
#----- Directories ------------
IDIR = -I/usr/local/include/ -I/usr/include/ -I$(core) #-I/home/jaeukk/install/include # include directory
ODIR = obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = ./# directory for executables

#	Executable
TARGET = $(EXEdir)PairStat.out

#----- Flags and libraries -----
GSLLIB = -L/usr/lib64 -lgsl -lgslcblas
NLOPTLIB = -L/usr/lib64 -lnlopt #-L/home/jaeukk/install/lib64 -lnlopt
LIBS = -lm $(GSLLIB) $(NLOPTLIB)
PREPRO = -DDELLA
endif

#==============================
# MY LAPTOP
#==============================
ifeq ($(MACHINE), HOME)
#echo $(MACHINE)
$(info	$(MACHINE))
#----- Directories ------------
IDIR = -I/usr/include/ -I/home/jaeukk/lib/gsl/include -I/usr/local/include/ -I$(core) # include directory
ODIR = obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = /home/jaeukk/EXC/# directory for executables

#	Executable
TARGET = $(EXEdir)PairStat.out

#----- Flags and libraries -----
GSLLIB = -L/home/jaeukk/lib/gsl/lib -lgsl -lgslcblas -Wl,-rpath,/home/jaeukk/lib/gsl/lib
NLOPTLIB = -L/usr/local/lib -lnlopt
LIBS = -lm $(GSLLIB) $(NLOPTLIB)
PREPRO = 
endif

#==============================
# CLUSTER 2
#==============================
ifeq ($(MACHINE), NURION)
#echo $(MACHINE)
$(info	$(MACHINE))
#----- Directories ------------
IDIR = -I/usr/include/ -I$(HOME)/install/include -I/usr/local/include/ -I$(core) # include directory
ODIR = obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = ~/EXC/# directory for executables

#	Executable
TARGET = $(EXEdir)PairStat.out

#----- Flags and libraries -----
GSLLIB = -L$(HOME)/install/lib -lgsl -lgslcblas -Wl,-rpath,$(HOME)/install/lib
NLOPTLIB = -L$(HOME)/install/lib -lnlopt
LIBS = -lm $(GSLLIB) $(NLOPTLIB)
PREPRO = 
endif

#==============================
#	Common flages
#==============================
CC = g++
DEBUG = -g
CFLAGS = -Wall -std=c++14 $(DEBUG) $(IDIR) -c -fopenmp -O3
LFLAGS = -Wall -std=c++14 $(DEBUG) $(LIBS) -fopenmp -O3

$(info	$(CFLAGS))
$(info	$(LFLAGS))

objs_core = GeometryVector.o RandomGenerator.o etc.o PeriodicCellList.o 
objs_pot = #Potential.o StructureOptimization.o 
objs_main = main.o Computation.o NumberVariance.o StructureFactor.o PairCorrelation.o ContactStatistics.o

DEPS = $(core)/GeometryVector.h $(core)/RandomGenerator.h $(core)/etc.h $(core)/PeriodicCellList.h StructureFactor.h PairCorrelation.h \
Computation.h NumberVariance.h ContactStatistics.h 
#  ../GenerateConfigs.h ../HardSphereFluids.h ../RandomSequentialAddition.h

OBJ_core = $(patsubst %.o, $(ODIR)%.o, $(objs_core))
OBJ_main = $(patsubst %.o, $(ODIR)%.o, $(objs_main))
OBJ = $(OBJ_core)  $(OBJ_main) 

$(info	$(OBJ))

SRC_core = $(patsubst %.o, $(core)%.cpp, $(objs_core))
#SRC_pot = $(patsubst %.o, ../%.cpp, $(objs_pot))
SRC_main = $(patsubst %.o, %.cpp, $(objs_main))
#SRC_misc = $(patsubst %.o, ../%.cpp, $(objs_misc))

$(info	$(SRC_core))
$(info	$(SRC_main))
#$(info	$(SRC_pot))


all: | $(ODIR)
all: $(OBJ) $(TARGET)

$(ODIR):
	mkdir -p $@

# compilations...
$(OBJ_core): $(ODIR)%.o: $(addprefix $(core), %.cpp) $(DEPS)
	echo $@ $<  
	$(CC) $(PREPRO) -o $@ $< $(CFLAGS)  

$(OBJ_main): $(ODIR)%.o: $(notdir %.cpp) $(DEPS)
	echo $@ $<
	$(CC) $(PREPRO) -o $@ $< $(CFLAGS) 

	echo "Start linking"

$(TARGET): $(OBJ)
	$(CC) $(PREPRO) -o $@ $(OBJ) $(LFLAGS)


# $(OBJ_pot): $(ODIR)%.o: $(addprefix ../, %.cpp) $(DEPS)
# 	echo $@ $<  
# 	$(CC) $(PREPRO) -o $@ $< $(CFLAGS)  


#$(OBJ): $(ODIR)%.o: $(SRC) $(DEPS)
#	echo $@ $<
#	$(CC) -o $@ $< $(CFLAGS) 

clean:
	rm -f $(ODIR)*.o $(TARGET)

.PHONY: all
all:
