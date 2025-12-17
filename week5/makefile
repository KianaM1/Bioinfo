# Test Makefile

usage:
	@echo "Running all targets"
	@echo make hello 
	@echo make bye

# Set the name variable.
NAME = Joe

# Include in future Makefiles after defining variables.
SHELL = bash
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY: usage hello bye

# Say hello to the name.
hello:
	echo "Hello ${NAME}!"

bye:
	echo "Goodbye world!"

all: hello bye