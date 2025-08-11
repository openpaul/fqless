# Simple Makefile for fqless
# A less-like viewer for FastQ files

# Compiler and flags
CXX = g++
CXXFLAGS = -g -std=c++11 -Wall -pedantic -O2
LDFLAGS = -lncurses -lz -pthread

# Program name
PROGRAM = fqless

# Source files
SOURCES = main.cpp fqless.cpp DNA.cpp fastq.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(PROGRAM)

# Build the program
$(PROGRAM): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $(PROGRAM)

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(PROGRAM)

# Install the program (optional)
install: $(PROGRAM)
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp $(PROGRAM) $(DESTDIR)$(PREFIX)/bin/

# Uninstall the program (optional)
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(PROGRAM)

# Set default prefix
PREFIX ?= /usr/local

# Phony targets
.PHONY: all clean install uninstall

# Dependencies (basic header dependencies)
main.o: main.cpp options.h DNA.h fastq.h fqless.h
fqless.o: fqless.cpp fqless.h options.h DNA.h fastq.h
DNA.o: DNA.cpp DNA.h options.h
fastq.o: fastq.cpp fastq.h DNA.h options.h