#include "s21_matrix_oop.h"

#include <cmath>
#include <cstring>
#include <stdexcept>

namespace s21 {

Matrix::Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < kMinSize || cols < kMinSize) {
    throw std::out_of_range("Error: Invalid number of rows or columns");
  }
  matrix_ = new double*[rows];
  try {
    matrix_[0] = new double[rows * cols]();
  } catch (std::bad_alloc& exc) {
    delete matrix_;
    throw;
  }
  for (int i = 0; i < rows; i++) {
    matrix_[i] = matrix_[0] + i * cols;
  }
}

Matrix::Matrix() : Matrix(kMinSize, kMinSize) {}

Matrix::Matrix(const Matrix& other) : Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < other.rows_; i++) {
    memcpy(matrix_[i], other.matrix_[i], other.cols_ * sizeof(double));
  }
}

Matrix::Matrix(Matrix&& other) noexcept
    : matrix_(other.matrix_), rows_(other.rows_), cols_(other.cols_) {
  other.matrix_ = nullptr;
  other.rows_ = other.cols_ = 0;
}

Matrix::~Matrix() {
  if (matrix_) {
    delete[] matrix_[0];
    delete[] matrix_;
  }
}

bool Matrix::EqMatrix(const Matrix& other) const noexcept {
  if (!SizeIsEqual(other)) {
    return false;
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > kEpsilon) {
        return false;
      }
    }
  }
  return true;
}

void Matrix::SumMatrix(const Matrix& other) {
  if (!SizeIsEqual(other)) {
    throw std::invalid_argument("Error: Sum. Different sizes of matrices");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void Matrix::SubMatrix(const Matrix& other) {
  if (!SizeIsEqual(other)) {
    throw std::invalid_argument("Error: Sub. Different sizes of matrices");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void Matrix::MulNumber(double num) noexcept {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void Matrix::MulMatrix(const Matrix& other) {
  if (cols_ != other.rows_ || rows_ != other.cols_) {
    throw std::invalid_argument(
        "Error: First matrix cols != Second matrix rows or vice versa");
  }
  Matrix tmp(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int m = 0; m < other.rows_; m++) {
        tmp.matrix_[i][j] += matrix_[i][m] * other.matrix_[m][j];
      }
    }
  }
  Swap(tmp);
}

double Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::invalid_argument("Error: Matrix is not square");
  }
  return GetDet();
}

Matrix Matrix::Transpose() const {
  Matrix result(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }
  return result;
}

Matrix Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::invalid_argument("Error: Matrix is not square");
  } else if (rows_ == kMinSize) {
    throw std::invalid_argument("Error: Rows value should be greater than 1");
  }
  Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      result.matrix_[i][j] =
          std::pow(-1.0, i + j) * GetMinor(i, j).Determinant();
    }
  }
  return result;
}

Matrix Matrix::InverseMatrix() const {
  Matrix result(cols_, rows_);
  double det = Determinant();
  if (fabs(det) < kEpsilon) {
    throw std::invalid_argument("Error: Determinant of this matrix is equal 0");
  }
  if (cols_ == kMinSize && rows_ == kMinSize) {
    result.matrix_[0][0] = 1.0 / matrix_[0][0];
  } else {
    result = CalcComplements().Transpose() * (1.0 / det);
  }
  return result;
}

bool Matrix::operator==(const Matrix& other) const noexcept {
  return EqMatrix(other);
}

Matrix Matrix::operator-(const Matrix& other) const {
  Matrix tmp(*this);
  tmp -= other;
  return tmp;
}

Matrix Matrix::operator+(const Matrix& other) const {
  Matrix tmp(*this);
  tmp += other;
  return tmp;
}

Matrix Matrix::operator*(const Matrix& other) const {
  Matrix tmp(*this);
  tmp *= other;
  return tmp;
}

Matrix& Matrix::operator+=(const Matrix& other) {
  SumMatrix(other);
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& other) {
  SubMatrix(other);
  return *this;
}

Matrix& Matrix::operator=(const Matrix& other) {
  if (this != &other) {
    Matrix tmp(other);
    Swap(tmp);
  }
  return *this;
}

Matrix& Matrix::operator=(Matrix&& other) noexcept {
  Matrix tmp(std::move(other));
  if (this != &other) {
    Swap(tmp);
  }
  return *this;
}

Matrix& Matrix::operator*=(double rhs) {
  MulNumber(rhs);
  return *this;
}

Matrix& Matrix::operator*=(const Matrix& other) {
  MulMatrix(other);
  return *this;
}

const double& Matrix::operator()(int row, int col) const {
  if (!IsValidIndices(row, col)) {
    throw std::out_of_range("Error: Invalid value of rows or columns");
  }
  return matrix_[row][col];
}

double& Matrix::operator()(int row, int col) {
  if (!IsValidIndices(row, col)) {
    throw std::out_of_range("Error: Invalid value of rows or columns");
  }
  return matrix_[row][col];
}

double* Matrix::operator[](int row) {
  return matrix_[row];
}

Matrix::iterator Matrix::begin() noexcept {
  return iterator(&matrix_[0][0]);
}

Matrix::const_iterator Matrix::begin() const noexcept {
  return const_iterator(&matrix_[0][0]);
}

Matrix::iterator Matrix::end() noexcept {
  return iterator(&matrix_[0][rows_ * cols_]);
}

Matrix::const_iterator Matrix::end() const noexcept {
  return const_iterator(&matrix_[0][rows_ * cols_]);
}

int Matrix::GetRows() const noexcept { return rows_; }

int Matrix::GetCols() const noexcept { return cols_; }

void Matrix::SetCols(int cols) {
  if (cols < kMinSize) {
    throw std::out_of_range("Error: Invalid number of cols");
  }
  Matrix tmp(rows_, cols);
  int cols_to_copy = std::min(cols_, cols);
  for (int i = 0; i < rows_; i++) {
    memcpy(tmp.matrix_[i], matrix_[i], cols_to_copy * sizeof(double));
  }
  Swap(tmp);
}

void Matrix::SetRows(int rows) {
  if (rows < kMinSize) {
    throw std::out_of_range("Error: Invalid number of rows");
  }
  Matrix tmp(rows, cols_);
  int rows_to_copy = std::min(rows_, rows);
  for (int i = 0; i < rows_to_copy; i++) {
    memcpy(tmp.matrix_[i], matrix_[i], cols_ * sizeof(double));
  }
  Swap(tmp);
}

Matrix Matrix::GetMinor(int row, int col) const {
  Matrix minor(rows_ - 1, cols_ - 1);
  int sub_i = 0, sub_j;
  for (int i = 0; i < rows_; i++) {
    sub_j = 0;
    if (i != row) {
      for (int j = 0; j < cols_; j++) {
        if (j != col) {
          minor.matrix_[sub_i][sub_j++] = matrix_[i][j];
        }
      }
      sub_i++;
    }
  }
  return minor;
}

double Matrix::GetDet() const {
  double determinant = 0.0;
  if (rows_ == 1) {
    determinant = matrix_[0][0];
  } else if (rows_ == 2) {
    determinant = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    double sign = 1.0;
    for (int i = 0; i < rows_; i++) {
      determinant += sign * matrix_[0][i] * GetMinor(0, i).GetDet();
      sign = -sign;
    }
  }
  return determinant;
}

inline void Matrix::Swap(Matrix& other) noexcept {
  std::swap(matrix_, other.matrix_);
  std::swap(rows_, other.rows_);
  std::swap(cols_, other.cols_);
}

inline bool Matrix::IsValidIndices(int row, int col) const noexcept {
  return (row < rows_) && (col < cols_) && (row >= 0) && (col >= 0);
}

inline bool Matrix::SizeIsEqual(const Matrix& rhs) const noexcept {
  return rows_ == rhs.rows_ && cols_ == rhs.cols_;
}

Matrix operator*(double lhs, const Matrix& rhs) {
  Matrix tmp(rhs);
  tmp *= lhs;
  return tmp;
}

Matrix operator*(const Matrix& lhs, double rhs) {
  Matrix tmp(lhs);
  tmp *= rhs;
  return tmp;
}

}  // namespace s21