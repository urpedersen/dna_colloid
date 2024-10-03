# Compiler and flags
CXX = g++
CXXFLAGS = -O3
SRC_DIR = src
BIN_DIR = bin
TARGET = $(BIN_DIR)/dna_colloid

# Find all .cpp files in the src directory
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Default target
all: $(TARGET)

# Create the binary directory and compile the project
$(TARGET): $(SRCS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

# Clean the build files
clean:
	rm -rf $(BIN_DIR)

# Phony targets to avoid conflict with files
.PHONY: all clean
