PWD        = $(shell pwd)

CXX 	  = g++
CFLAGS    = -Wall -Wno-parentheses -O3
LFLAGS    = -Wall -lz --static

.PHONY: all clean


TARGET     = decoder

# Define the source files
CSRCS      = $(wildcard $(PWD)/*.cc) 

# Automatically derive object files from source files
COBJS      = $(CSRCS:.cc=.o)

# Default target: builds the executable
all: $(TARGET)

# Rule to link object files into the executable
$(TARGET): $(COBJS)
	$(CXX) $(COBJS) -o $(TARGET)

# Rule to compile individual C source files into object files
%.o: %.cc
	$(CXX) $(CFLAGS) -c -o $@ $<


# Clean target: removes compiled files
clean:
	rm -f $(TARGET) $(COBJS)