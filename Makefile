TARGET := lj_gas.out
SRC_DIR := src
BIN_DIR := bin

SRCS := $(wildcard $(SRC_DIR)/*.cpp)
CXX := clang++
CXXFLAGS := -Wall -Wextra -std=c++20 -Iinclude

TARGET_BIN := $(BIN_DIR)/$(TARGET)

all: $(TARGET_BIN)

$(TARGET_BIN):
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $@

clean:
	rm -rf $(BIN_DIR)

.PHONY: all clean

