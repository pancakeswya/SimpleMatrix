OS         := $(shell uname -s)

TEST_LIB   := gtest

NAME       := s21_matrix_oop.a

SRC_DIR    := src
INC_DIR    := includes
SRCS       := src/s21_matrix_oop.cc
OBJS       := $(SRCS:.cc=.o)

CXX        := g++
CXXFLAGS   := -Wall -Wextra -Werror -std=c++17 -I $(INC_DIR)

AR         := ar
ARFLAGS    := -r -c -s

TEST_LDLIB := $(addprefix -l,$(TEST_LIB))

GCOV       := --coverage
LCOV       := lcov --no-external -c

VALGRIND   := valgrind --tool=memcheck --trace-children=yes --track-origins=yes --leak-check=full

OPEN       := $(if $(filter Linux,$(OS)),xdg-open,open)

TEST_SRCS  := src/s21-matrix-test.cc
TEST_NAME  := s21-matrix-test
LCOV_NAME  := s21-matrix.info

REPORT_DIR := report

LDFLAGS    := $(addprefix -L,$(SRC_DIR))

CP         := cp
RM         := rm -rf
MAKEFLAGS  += --no-print-directory

all: $(NAME)

$(NAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

test: $(NAME)
	$(CXX) $(CXXFLAGS) $(TEST_SRCS) $(TEST_LDLIB) $(LDFLAGS) $< -o $(TEST_NAME)
	./$(TEST_NAME)

check-valgrind: test
	CK_FORK=NO $(VALGRIND) ./$(TEST_NAME)

check-style:
	clang-format -style=google -n $(SRC_DIR)/*.cc  $(INC_DIR)/*.h

gcov_report: $(NAME)
	$(CXX) $(CXXFLAGS) $(GCOV) -I $(INC_DIR) $(TEST_SRCS) $(SRCS) $(TEST_LDLIB) $(LDFLAGS) $< -o $(TEST_NAME)
	./$(TEST_NAME)
	$(LCOV) -t $(TEST_NAME) -d . -o $(LCOV_NAME)
	genhtml $(LCOV_NAME) -o $(REPORT_DIR)
	$(OPEN) $(REPORT_DIR)/index.html

clean:
	$(RM) $(NAME)
	$(RM) $(OBJS)
	$(RM) $(TEST_NAME)

fclean: clean
	$(RM) $(LCOV_NAME)
	$(RM) $(REPORT_DIR)
	$(RM) *.gcno *.gcda

rebuild:
	$(MAKE) clean
	$(MAKE) all