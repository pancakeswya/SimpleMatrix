OS         := $(shell uname -s)

TEST_LIB   := gtest

CXX        := g++
CXXFLAGS   := -Wall -Wextra -Werror -std=c++17

AR         := ar
ARFLAGS    := -r -c -s

TEST_LDLIB := $(addprefix -l,$(TEST_LIB))

GCOV       := --coverage
LCOV       := lcov --no-external -c

VALGRIND   := valgrind --tool=memcheck --trace-children=yes --track-origins=yes --leak-check=full

OPEN       := $(if $(filter Linux,$(OS)),xdg-open,open)

CP         := cp
RM         := rm -rf
MAKEFLAGS  += --no-print-directory