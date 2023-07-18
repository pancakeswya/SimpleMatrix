include Makefile.mk

NAME       := s21_matrix_oop.a

SRC_DIR    := src
SRCS       := src/s21_matrix_oop.cc
OBJS       := $(SRCS:.cc=.o)

TEST_SRCS  := src/s21-matrix-test.cc
TEST_NAME  := s21-matrix-test
LCOV_NAME  := s21-matrix.info

REPORT_DIR := report

LDFLAGS    := $(addprefix -L,$(SRC_DIR))

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
	clang-format -style=google -n $(SRC_DIR)/*.cc  $(SRC_DIR)/*.h

gcov_report: $(NAME)
	$(CXX) $(CXXFLAGS) $(GCOV) $(TEST_SRCS) $(SRCS) $(TEST_LDLIB) $(LDFLAGS) $< -o $(TEST_NAME)
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