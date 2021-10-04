# MAKEFILE FOR fildsim.f90
# Use \ to concatenate lines and && to separate commands

# DIRECTORIES
BASE_DIR=$(PWD)
SRC_DIR=$(PWD)/src/
BIN_DIR=$(PWD)/bin/

# COMPILER TO USE
COMP=gfortran

# OPTIMIZATION SETTINGS
OPT=-O3 -ftree-vectorize
FFLAGS= -g -fbacktrace -Wall #-fopenmp -g #-Wall #-fbacktrace -O0 -fcheck=all

# FORTRAN MODULES
MODULES= sinpa_module.f90

# TARGET
TARGET= sinpa.f90

# EXECUTABLE NAME
EXE=SINPA.go


all: sinpa_mod sinpa

sinpa_mod:
	cd $(SRC_DIR) && \
	$(COMP) $(OPT) -c $(MODULES)

sinpa:
	cd $(SRC_DIR) && \
	$(COMP) $(OPT) $(FFLAGS) $(MODULES:.f90=.o) $(TARGET) -o $(EXE)

	mv $(SRC_DIR)*.o $(BIN_DIR)
	mv $(SRC_DIR)*.mod $(BIN_DIR)
	mv $(SRC_DIR)*.go $(BIN_DIR)

clean:
	rm -f $(BIN_DIR)*.go
	rm -f $(BIN_DIR)*.o
	rm -f $(BIN_DIR)*.mod
