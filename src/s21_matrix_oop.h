#ifndef CPP_S21_MATRIXPLUS_1_SRC_S21MATRIX_OOP_H_
#define CPP_S21_MATRIXPLUS_1_SRC_S21MATRIX_OOP_H_

namespace s21 {

class Matrix {
 public:
  explicit Matrix(int, int);
  Matrix();
  Matrix(const Matrix &);
  Matrix(Matrix &&) noexcept;
  ~Matrix();
  bool EqMatrix(const Matrix &) const noexcept;
  void SumMatrix(const Matrix &);
  void SubMatrix(const Matrix &);
  void MulNumber(double) noexcept;
  void MulMatrix(const Matrix &);
  double Determinant() const;
  Matrix Transpose() const;
  Matrix CalcComplements() const;
  Matrix InverseMatrix() const;
  bool operator==(const Matrix &) const noexcept;
  Matrix operator-(const Matrix &) const;
  Matrix operator+(const Matrix &) const;
  Matrix operator*(const Matrix &) const;
  Matrix &operator=(const Matrix &);
  Matrix &operator=(Matrix &&) noexcept;
  Matrix &operator+=(const Matrix &);
  Matrix &operator-=(const Matrix &);
  Matrix &operator*=(const Matrix &);
  Matrix &operator*=(double);
  double operator()(int, int) const;
  double &operator()(int, int);
  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetCols(int);
  void SetRows(int);

 private:
  static constexpr int minSize = 1;
  static constexpr double epsilon = 1e-07;
  Matrix GetMinor(int, int) const;
  double GetDet() const;
  double &GetElement(int, int) const;
  void SwapMatrix(Matrix &) noexcept;
  bool ValidIndices(int, int) const noexcept;
  bool SizeIsEqual(const Matrix &) const noexcept;
  double **matrix_;
  int rows_, cols_;
};

Matrix operator*(const Matrix &, double);
Matrix operator*(double, const Matrix &);

}  // namespace s21

#endif  // CPP_S21_MATRIXPLUS_1_SRC_S21MATRIX_OOP_H_