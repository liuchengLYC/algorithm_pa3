CXX       := g++
CXXFLAGS  := -std=c++17 -O2 -Wall -Wextra -pedantic

SRC_DIR   := src
BIN_DIR   := bin
TARGET    := $(BIN_DIR)/router

SRCS      := $(wildcard $(SRC_DIR)/*.cpp)
OBJS      := $(patsubst $(SRC_DIR)/%.cpp,$(SRC_DIR)/%.o,$(SRCS))
DEPS      := $(OBJS:.o=.d)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

-include $(DEPS)

clean:
	rm -f $(SRC_DIR)/*.o $(SRC_DIR)/*.d $(TARGET)
