# MAKEFILE FOR fildsim.f90
# Use \ to concatenate lines and && to separate commands

# DIRECTORIES
BASE_DIR=$(PWD)
SRC_DIR=$(PWD)/src/
BIN_DIR=$(PWD)/bin/

# COMPILER TO USE
COMP=gfortran

# OPTIMIZATION SETTINGS
OPT_TOKs=-O3 -ftree-vectorize -mavx512f -mavx512cd
OPT=-O3 -ftree-vectorize
FFLAGS= -g -fbacktrace -Wall #-fopenmp -g #-Wall #-fbacktrace -O0 -fcheck=all

# FORTRAN MODULES
MODULES= sinpa_module.f90

# TARGET
TARGET= sinpa.f90

# EXECUTABLE NAME
EXE=SINPA.go


all: sinpa_mod sinpa
all_cluster: sinpa_mod_cluster sinpa_cluster
debug: sinpa_mod sinpa_debug

sinpa_mod:
	cd $(SRC_DIR) && \
	$(COMP) $(OPT) -c $(MODULES) 

sinpa_mod_cluster:
	cd $(SRC_DIR) && \
	$(COMP) $(OPT_TOKs) -c $(MODULES)

sinpa:
	cd $(SRC_DIR) && \
	$(COMP) $(OPT) $(MODULES:.f90=.o) $(TARGET) -o $(EXE) 

	mv $(SRC_DIR)*.o $(BIN_DIR)
	mv $(SRC_DIR)*.mod $(BIN_DIR)
	mv $(SRC_DIR)*.go $(BIN_DIR)

sinpa_debug:
	cd $(SRC_DIR) && \
	$(COMP)  $(FFLAGS) $(MODULES:.f90=.o) $(TARGET) -o $(EXE)

	mv $(SRC_DIR)*.o $(BIN_DIR)
	mv $(SRC_DIR)*.mod $(BIN_DIR)
	mv $(SRC_DIR)*.go $(BIN_DIR)

sinpa_cluster:
	cd $(SRC_DIR) && \
	$(COMP) $(OPT_TOKs) $(MODULES:.f90=.o) $(TARGET) -o $(EXE)

	mv $(SRC_DIR)*.o $(BIN_DIR)
	mv $(SRC_DIR)*.mod $(BIN_DIR)
	mv $(SRC_DIR)*.go $(BIN_DIR)

clean:
	rm -f $(BIN_DIR)*.go
	rm -f $(BIN_DIR)*.o
	rm -f $(BIN_DIR)*.mod
