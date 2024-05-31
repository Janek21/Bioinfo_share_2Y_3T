# Makefile with rules to compile C and C++ programs, generating executable
# files with names ending in _c for C programs and _cpp for C++ files.

CC := gcc       # Compiler for C
CXX := g++      # Compiler for C++
CXXFLAGS := -std=c++11  # C++11 standard

.PHONY: all clean

# List of C source files
C_SOURCES := $(wildcard *.c)
C_EXECUTABLES := $(C_SOURCES:.c=_c)

# List of C++ source files
CPP_SOURCES := $(wildcard *.cpp)
CPP_EXECUTABLES := $(CPP_SOURCES:.cpp=_cpp)

all: $(C_EXECUTABLES) $(CPP_EXECUTABLES)

%_c: %.c
	$(CC) -o $@ $<

%_cpp: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(C_EXECUTABLES) $(CPP_EXECUTABLES)

