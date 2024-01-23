CC = gcc
CFLAGS = -Wall -I./include -Wno-deprecated-non-prototype

SRC_DIR = src
LIB_DIR = lib
BUILD_DIR = build
BIN_DIR = bin

# List of source files
SRC = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(SRC_DIR)/python/*.c)
LIB_SRC = $(wildcard $(LIB_DIR)/*.c)

# List of object files derived from source files
OBJ = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRC)) $(patsubst $(LIB_DIR)/%.c,$(BUILD_DIR)/%.o,$(LIB_SRC))

# Final executable target
EXECUTABLE = $(BIN_DIR)/num-int

# Libraries to link against
LIBS = -L./lib/ -lm -lgsl

# Default target, build the executable
all: $(BIN_DIR) $(BUILD_DIR) $(EXECUTABLE)

# Rule to build the executable
$(EXECUTABLE): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LIBS)

# Rule to build object files from source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(wildcard $(SRC_DIR)/*.h)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/python/%.c $(wildcard $(SRC_DIR)/python/*.h) | $(BUILD_DIR)/python
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(LIB_DIR)/%.c $(wildcard $(LIB_DIR)/*.h)
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to create build directory
$(BUILD_DIR):
	mkdir -p $@

# Rule to create build/python directory
$(BUILD_DIR)/python:
	mkdir -p $@

# Rule to create bin directory
$(BIN_DIR):
	mkdir -p $@

# Rule to clean up generated files
clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/python/*.o $(EXECUTABLE)

# Marking non-file targets as phony
.PHONY: all clean
