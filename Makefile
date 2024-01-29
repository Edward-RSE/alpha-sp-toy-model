CC = mpicc
CFLAGS = -Wall -O3 -Wno-deprecated-non-prototype -std=gnu99

# List of directories
SOURCE_DIR = ./src
LIB_DIR = ./lib
INCLUDE_DIR = ./inc
OBJ_DIR = ./obj
BIN_DIR = ./bin

# List of source files
PYTHON_SOURCE = $(wildcard $(SOURCE_DIR)/python/*.c)
NUM_INT_SOURCE = $(wildcard $(SOURCE_DIR)/num-int/*.c)
NODE_SHARE_SOURCE = $(wildcard $(SOURCE_DIR)/node-share/*.c)

# List of object files
PYTHON_OBJECTS = $(patsubst $(SOURCE_DIR)/%.c,$(OBJ_DIR)/%.o,$(PYTHON_SOURCE))
NUM_INT_OBJECTS = $(patsubst $(SOURCE_DIR)/%.c,$(OBJ_DIR)/%.o,$(NUM_INT_SOURCE))
NODE_SHARE_OBJECTS = $(patsubst $(SOURCE_DIR)/%.c,$(OBJ_DIR)/%.o,$(NODE_SHARE_SOURCE))

NUM_INT_EXE = $(BIN_DIR)/num-int
NODE_SHARE_EXE = $(BIN_DIR)/node-share
LIBS = -L$(LIB_DIR) -lm -lgsl

all: directories $(NUM_INT_EXE) $(NODE_SHARE_EXE)

# Rule to create objects
$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.c
	@mkdir -p $(@D)  # Ensure the output directory exists
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

# Rule to create num-int executable
$(NUM_INT_EXE): $(NUM_INT_OBJECTS) $(PYTHON_OBJECTS)
	$(CC) $(NUM_INT_OBJECTS) $(PYTHON_OBJECTS) -o $@ $(LIBS)

# Rule to create node-share executable
$(NODE_SHARE_EXE): $(NODE_SHARE_OBJECTS) $(PYTHON_OBJECTS)
	$(CC) $(NODE_SHARE_OBJECTS) $(PYTHON_OBJECTS) -o $@ $(LIBS)

# Rule to create directories
directories:
	mkdir -p $(OBJ_DIR) $(BIN_DIR)

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean directories
