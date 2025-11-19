#=============================-
# Makefile to compile a simple file format manager
# for converting a ConfigPack file to txt-format files, and vice versa.
#
# Author: Jaeuk Kim
# Date: September 2020
# Email: phy000.kim@gmail.com
#=============================

MACHINE = HOME
#MACHINE = DELLA

#==============================
# DELLA
#==============================
ifeq ($(MACHINE), DELLA)
#----- Directories ------------
IDIR = -I/usr/local/include/ -I$(core)# include directory
ODIR = obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = ~/EXC/# directory for executables

#	Executable
TARGET = $(EXEdir)manager.out

#----- Flags and libraries -----
GSLLIB = -lgsl -lgslcblas
LIBS = -lm $(GSLLIB)
endif

#==============================
# MY LAPTOP
#==============================
ifeq ($(MACHINE), HOME)
$(info	$(MACHINE))
#----- Directories ------------
core = ../cores/

IDIR = -I/usr/include/ -I/home/jaeukk/lib/gsl/include -I$(core)# include directory
ODIR = obj/# object directory
#SDIR = ./# source directory
#Idir = ./# header directory
EXEdir = ~/EXC/# directory for executables

#	Executable
TARGET = $(EXEdir)manager.out #te.out## output file name

#----- Flags and libraries -----
GSLLIB = -L/home/jaeukk/lib/gsl/lib -lgsl -lgslcblas -lm -Wl,-rpath,/home/jaeukk/lib/gsl/lib
LIBS = -lm $(GSLLIB) 
endif


#==============================
#	Common flages
#==============================
CC = g++
DEBUG = -g
CFLAGS = -Wall -Wno-reorder -std=c++14 $(DEBUG) $(IDIR) -c -O3 -fopenmp
LFLAGS = -Wall -Wno-reorder -std=c++14 $(DEBUG) $(LIBS) -O3 -fopenmp

#==============================
#	Flags and libraries
#==============================
#CC = g++
#DEBUG = -g
#GSLLIB = -L/home/jaeukk/lib/gsl/lib -lgsl -lgslcblas -lm
#LIBS = -lm $(GSLLIB) 

$(info	$(CFLAGS))
$(info	$(LFLAGS))


objs_core = GeometryVector.o PeriodicCellList.o RandomGenerator.o etc.o
objs_netw = main.o

DEPS = $(core)/GeometryVector.h $(core)/PeriodicCellList.h $(core)/RandomGenerator.h $(core)/etc.h

OBJ_core = $(patsubst %.o, $(ODIR)%.o, $(objs_core))
OBJ_netw = $(patsubst %.o, $(ODIR)%.o, $(objs_netw))  
OBJ = $(OBJ_core) $(OBJ_netw)
$(info  $(OBJ))

SRC_core = $(patsubst %.o, $(core)/%.cpp, $(objs_core))
SRC_netw = $(patsubst %.o, %.cpp, $(objs_netw))

$(info  $(SRC_core))
$(info  $(SRC_netw))

all: $(OBJ) $(TARGET)

$(TARGET): $(OBJ)
	$(CC) -o $@ $(OBJ) $(LFLAGS)

# compilations...
$(OBJ_core): $(ODIR)%.o: $(addprefix $(core)/, %.cpp) $(DEPS)
	echo $@ $<  
	$(CC) $(PREPROC) -o $@ $< $(CFLAGS)  

$(OBJ_netw): $(ODIR)%.o: $(notdir %.cpp) $(DEPS)
	echo $@ $<
	$(CC) $(PREPROC) -o $@ $< $(CFLAGS) 


clean:
	rm -f $(ODIR)*.o $(TARGET)

.PHONY: all
