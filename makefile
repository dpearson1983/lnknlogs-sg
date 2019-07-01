CXX = g++
CXXOPTS = -march=native -mtune=native -O3 -Wall
LDFLAGS = -lgsl -lgslcblas -lm -lfftw3 -lfftw3_omp -lchealpix
OMP = -fopenmp
SRC_DIR := source
OBJ_DIR := obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

lnknlogs_sg: $(OBJ_FILES)
	$(CXX) $(OMP) $(LDFLAGS) $^ -o $@
	mkdir -p $(HOME)/bin
	cp $@ $(HOME)/bin/$@
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(OMP) $(CXXOPTS) -c $< -o $@

clean:
	rm obj/*.o
	rm $(HOME)/bin/lnknlogs_sg
