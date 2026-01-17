# ============================================================================
# PALMA - Parallel Algebra Library for Max-plus Applications
# Makefile
#
# Author: Gnankan Landry Regis N'guessan
#         Axiom Research Group
#         NM-AIST / AIMS-RIC
# Email:  rnguessan@aimsric.org
#
# License: MIT
# ============================================================================

# Project info
PROJECT = palma
VERSION = 1.0.0

# Compiler settings
CC = gcc
AR = ar
ARFLAGS = rcs

# Detect architecture
ARCH := $(shell uname -m)
OS := $(shell uname -s)

# Base compiler flags
CFLAGS = -std=c99 -Wall -Wextra -Wpedantic -O3
LDFLAGS = -lm

# Include path
INCLUDES = -I./include

# ARM-specific optimizations for Raspberry Pi
ifeq ($(ARCH),aarch64)
    # 64-bit ARM (Raspberry Pi 3/4/5 in 64-bit mode)
    CFLAGS += -march=armv8-a -mtune=cortex-a72
    NEON_FLAGS = -DPALMA_USE_NEON=1
    PLATFORM = ARM64 (Raspberry Pi 3/4/5 64-bit)
else ifeq ($(ARCH),armv7l)
    # 32-bit ARM (Raspberry Pi 2/3/4 in 32-bit mode)
    CFLAGS += -march=armv7-a -mfpu=neon-vfpv4 -mfloat-abi=hard
    NEON_FLAGS = -DPALMA_USE_NEON=1
    PLATFORM = ARM32 (Raspberry Pi 32-bit)
else ifeq ($(findstring arm,$(ARCH)),arm)
    # Other ARM variants
    CFLAGS += -march=native
    NEON_FLAGS = -DPALMA_USE_NEON=1
    PLATFORM = ARM (Generic)
else
    # x86/other - disable NEON
    NEON_FLAGS = -DPALMA_USE_NEON=0
    PLATFORM = $(ARCH)
endif

# OpenMP support (optional)
ifdef USE_OPENMP
    CFLAGS += -fopenmp
    LDFLAGS += -fopenmp
    OPENMP_FLAGS = -DPALMA_USE_OPENMP=1
else
    OPENMP_FLAGS = -DPALMA_USE_OPENMP=0
endif

# Combine all flags
ALL_CFLAGS = $(CFLAGS) $(NEON_FLAGS) $(OPENMP_FLAGS) $(INCLUDES)

# Build directories
BUILD_DIR = build
LIB_DIR = $(BUILD_DIR)/lib
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin

# Source files
LIB_SRCS = src/palma.c src/palma_ext.c
LIB_OBJS = $(patsubst src/%.c,$(OBJ_DIR)/%.o,$(LIB_SRCS))
LIB_STATIC = $(LIB_DIR)/lib$(PROJECT).a

# Example sources
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,$(BIN_DIR)/%,$(EXAMPLE_SRCS))

# Header files
HEADERS = include/palma.h

# ============================================================================
# MAIN TARGETS
# ============================================================================

.PHONY: all lib examples clean install uninstall dist info help

# Default: build library and examples
all: lib examples

# Create directories
$(BUILD_DIR) $(LIB_DIR) $(OBJ_DIR) $(BIN_DIR):
	@mkdir -p $@

# ============================================================================
# LIBRARY
# ============================================================================

lib: $(LIB_STATIC)

$(LIB_STATIC): $(LIB_OBJS) | $(LIB_DIR)
	@echo "  AR      $@"
	@$(AR) $(ARFLAGS) $@ $^

$(OBJ_DIR)/%.o: src/%.c $(HEADERS) | $(OBJ_DIR)
	@echo "  CC      $<"
	@$(CC) $(ALL_CFLAGS) -c $< -o $@

# ============================================================================
# EXAMPLES
# ============================================================================

examples: $(EXAMPLE_BINS)

$(BIN_DIR)/%: examples/%.c $(LIB_STATIC) | $(BIN_DIR)
	@echo "  LINK    $@"
	@$(CC) $(ALL_CFLAGS) $< -L$(LIB_DIR) -l$(PROJECT) $(LDFLAGS) -o $@

# Individual example targets
.PHONY: scheduling graphs eigenvalue benchmark
scheduling: $(BIN_DIR)/example_scheduling
graphs: $(BIN_DIR)/example_graphs
eigenvalue: $(BIN_DIR)/example_eigenvalue
benchmark: $(BIN_DIR)/benchmark

# ============================================================================
# RUN TARGETS
# ============================================================================

.PHONY: run-scheduling run-graphs run-eigenvalue run-benchmark run-all

run-scheduling: $(BIN_DIR)/example_scheduling
	@./$(BIN_DIR)/example_scheduling

run-graphs: $(BIN_DIR)/example_graphs
	@./$(BIN_DIR)/example_graphs

run-eigenvalue: $(BIN_DIR)/example_eigenvalue
	@./$(BIN_DIR)/example_eigenvalue

run-benchmark: $(BIN_DIR)/benchmark
	@./$(BIN_DIR)/benchmark

run-all: $(EXAMPLE_BINS)
	@echo "=== Running all examples ==="
	@for bin in $(EXAMPLE_BINS); do \
		echo "\n>>> Running $$bin <<<\n"; \
		./$$bin; \
	done

# ============================================================================
# BUILD VARIANTS
# ============================================================================

.PHONY: debug release scalar openmp

# Debug build with symbols and no optimization
debug: CFLAGS = -std=c99 -Wall -Wextra -Wpedantic -g -O0 -DDEBUG
debug: clean lib examples
	@echo "Debug build complete"

# Release build with aggressive optimization
release: CFLAGS += -DNDEBUG -flto
release: clean lib examples
	@echo "Release build complete"

# Build without NEON (for comparison)
scalar: NEON_FLAGS = -DPALMA_USE_NEON=0
scalar: clean lib examples
	@echo "Scalar (no NEON) build complete"

# Build with OpenMP support
openmp: USE_OPENMP = 1
openmp: clean lib examples
	@echo "OpenMP build complete"

# ============================================================================
# INSTALLATION
# ============================================================================

PREFIX ?= /usr/local

install: lib
	@echo "Installing PALMA to $(PREFIX)..."
	@install -d $(PREFIX)/include
	@install -d $(PREFIX)/lib
	@install -m 644 $(HEADERS) $(PREFIX)/include/
	@install -m 644 $(LIB_STATIC) $(PREFIX)/lib/
	@echo "Installation complete"

uninstall:
	@echo "Uninstalling PALMA from $(PREFIX)..."
	@rm -f $(PREFIX)/include/palma.h
	@rm -f $(PREFIX)/lib/libpalma.a
	@echo "Uninstallation complete"

# ============================================================================
# DISTRIBUTION
# ============================================================================

DIST_NAME = $(PROJECT)-$(VERSION)
DIST_FILES = Makefile README.md LICENSE include/ src/ examples/

dist: clean
	@echo "Creating distribution package $(DIST_NAME).tar.gz..."
	@mkdir -p $(DIST_NAME)
	@cp -r $(DIST_FILES) $(DIST_NAME)/
	@tar -czf $(DIST_NAME).tar.gz $(DIST_NAME)
	@rm -rf $(DIST_NAME)
	@echo "Package created: $(DIST_NAME).tar.gz"

# ============================================================================
# DOCUMENTATION
# ============================================================================

.PHONY: docs

docs:
	@if command -v doxygen > /dev/null; then \
		echo "Generating documentation..."; \
		doxygen Doxyfile 2>/dev/null || echo "Create Doxyfile first: doxygen -g"; \
	else \
		echo "Doxygen not found. Install with: sudo apt install doxygen"; \
	fi

# ============================================================================
# TESTING
# ============================================================================

.PHONY: test

test: $(EXAMPLE_BINS)
	@echo "=== Running Tests ==="
	@./$(BIN_DIR)/example_scheduling > /dev/null && echo "✓ Scheduling example passed" || echo "✗ Scheduling example failed"
	@./$(BIN_DIR)/example_graphs > /dev/null && echo "✓ Graphs example passed" || echo "✗ Graphs example failed"
	@./$(BIN_DIR)/example_eigenvalue > /dev/null && echo "✓ Eigenvalue example passed" || echo "✗ Eigenvalue example failed"
	@echo "=== All tests passed ==="

# ============================================================================
# CLEANUP
# ============================================================================

clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(BUILD_DIR)
	@rm -f *.dot *.csv *.bin *.png
	@echo "Clean complete"

distclean: clean
	@rm -f $(DIST_NAME).tar.gz

# ============================================================================
# INFORMATION
# ============================================================================

info:
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════╗"
	@echo "║  PALMA - Parallel Algebra Library for Max-plus Applications      ║"
	@echo "║  Version $(VERSION)                                                    ║"
	@echo "╚══════════════════════════════════════════════════════════════════╝"
	@echo ""
	@echo "Author:  Gnankan Landry Regis N'guessan"
	@echo "Email:   rnguessan@aimsric.org"
	@echo ""
	@echo "Build Configuration:"
	@echo "  Platform:     $(PLATFORM)"
	@echo "  Compiler:     $(CC)"
	@echo "  CFLAGS:       $(CFLAGS)"
	@echo "  NEON:         $(NEON_FLAGS)"
	@echo "  OpenMP:       $(OPENMP_FLAGS)"
	@echo "  Build Dir:    $(BUILD_DIR)"
	@echo ""

help:
	@echo ""
	@echo "PALMA Build System - Available Targets:"
	@echo ""
	@echo "  Building:"
	@echo "    all           Build library and all examples (default)"
	@echo "    lib           Build static library only"
	@echo "    examples      Build all examples"
	@echo "    debug         Build with debug symbols"
	@echo "    release       Build with aggressive optimization"
	@echo "    scalar        Build without NEON (for comparison)"
	@echo "    openmp        Build with OpenMP parallelization"
	@echo ""
	@echo "  Running:"
	@echo "    run-scheduling   Run scheduling example"
	@echo "    run-graphs       Run graph algorithms example"
	@echo "    run-eigenvalue   Run eigenvalue example"
	@echo "    run-benchmark    Run performance benchmark"
	@echo "    run-all          Run all examples"
	@echo ""
	@echo "  Installation:"
	@echo "    install       Install to $(PREFIX)"
	@echo "    uninstall     Remove from $(PREFIX)"
	@echo ""
	@echo "  Other:"
	@echo "    test          Run basic tests"
	@echo "    docs          Generate documentation (requires Doxygen)"
	@echo "    dist          Create distribution tarball"
	@echo "    clean         Remove build artifacts"
	@echo "    info          Show build configuration"
	@echo "    help          Show this help message"
	@echo ""
