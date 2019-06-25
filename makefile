CXX = g++
CXXOPTS = -march=native -mtune=native -O3
LDFLAGS = -lgsl -lgslcblas -lm -lfftw3 -lfftw3_omp -fopenmp
SRC_DIR := source
OBJ_DIR := obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

lnknlogs_sg: $(OBJ_FILES)
	$(CXX) $(LDFLAGS) $^ -o $@
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXOPTS) -c $< -o $@
