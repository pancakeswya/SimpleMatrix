#ifndef STL_MATRIX_H_
#define STL_MATRIX_H_

#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <cassert>
#include <iostream>

// template matrix in stl style made for fun
namespace stl {

#if __cplusplus >= 202002L
template<typename Tp>
concept arithmetic = std::integral<Tp> or std::floating_point<Tp>;

template<typename Tp>
requires arithmetic<Tp>
#else
template<typename Tp, typename  = typename std::enable_if_t<std::is_arithmetic_v<Tp>, Tp>>
#endif
class matrix {
 public:
  using value_type = Tp;
  using reference = Tp&;
  using const_reference = const Tp&;
  using pointer = Tp*;
  using const_pointer = const Tp*;
  using size_type = size_t;

  using iterator = typename std::vector<Tp>::iterator;
  using const_iterator = typename std::vector<Tp>::const_iterator;

  // not ideal comparators
  template<typename T>
  std::enable_if_t<std::is_floating_point_v<T>, bool>
  is_equal(const T& rhs, const T& lhs) const noexcept {
    return std::abs(rhs - lhs) <= static_cast<T>(0.000001);
  }

  template<typename T>
  std::enable_if_t<!std::is_floating_point_v<T>, bool>
  is_equal(const T& rhs, const T& lhs) const noexcept {
    return rhs == lhs;
  }


  matrix() = default;
  matrix(std::initializer_list<std::initializer_list<Tp>> const& list)
    : rows_(list.size()), cols_(list.begin()->size()) {
    for (auto& i : list) {
      assert(i.size() == cols_);
      for (auto& j : i) {
        data_.push_back(j);
      }
    }
  }

  matrix(size_type rows, size_type cols)
    : rows_(rows), cols_(cols), data_(rows * cols) {}

  matrix(const matrix& other)
    : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {}

  matrix(matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), data_(std::move(other.data_)) {
    other.rows_ = other.cols_ = 0;
  }

  matrix& operator=(const matrix& other) {
    if (this != &other) {
      matrix tmp(other);
      swap(tmp);
    }
    return *this;
  }

  matrix& operator=(matrix&& other) noexcept {
    matrix tmp(std::move(other));
    if (this != &other) {
      swap(tmp);
    }
    return *this;
  }

  void swap(matrix& other) noexcept {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(data_, other.data_);
  }

  iterator begin() noexcept {
    return data_.begin();
  }

  const_iterator begin() const noexcept {
    return data_.begin();
  }

  iterator end() noexcept {
    return data_.end();
  }

  const_iterator end() const noexcept {
    return data_.end();
  }

  void clear() noexcept {
    data_.clear();
    cols_ = rows_ = 0;
  }

  pointer operator[](size_type row) noexcept {
    return &data_[row * cols_];
  }

  const_pointer operator[](size_type row) const noexcept {
    return &data_[row * cols_];
  }

  reference at(size_type row, size_type col) {
    if (row > rows_ || col > cols_) {
      throw std::out_of_range("operator[]. Row or col is out of range");
    }
    return data_[row * cols_ + col];
  }

  const_reference at(size_type row, size_type col) const {
    if (row > rows_ || col > cols_) {
      throw std::out_of_range("operator[]. Row or col is out of range");
    }
    return data_[row * cols_ + col];
  }

  size_type rows() const noexcept {
    return rows_;
  }

  size_type cols() const noexcept {
    return cols_;
  }

  bool equal(const matrix& other) const noexcept {
    if (!size_is_equal(other)) {
      return false;
    }
    for (size_type i = 0; i < data_.size(); ++i) {
      if (!is_equal(data_[i], other.data_[i])) {
        return false;
      }
    }
    return true;
  }

  void sum(const matrix& other) {
    if (!size_is_equal(other)) {
      throw std::invalid_argument("Sum. Different sizes of matrices");
    }
    for (size_type i = 0; i < data_.size(); ++i) {
        data_[i] += other.data_[i];
    }
  }

  void sub(const matrix& other) {
    if (!size_is_equal(other)) {
      throw std::invalid_argument("Sub. Different sizes of matrices");
    }
    for (int i = 0; i < data_.size(); ++i) {
        data_[i] -= other.data_[i];
    }
  }

  void multiply(const Tp& value) noexcept {
    for (int i = 0; i < data_.size(); ++i) {
        data_[i] *= value;
    }
  }

  void multiply(const matrix& other) {
    if (cols_ != other.rows_ || rows_ != other.cols_) {
      throw std::invalid_argument(
          "Multiply. First matrix cols != Second matrix rows or vice versa");
    }
    matrix tmp(rows_, other.cols_);
    for (size_type i = 0; i < rows_; ++i) {
      for (size_type j = 0; j < other.cols_; ++j) {
        for (size_type m = 0; m < other.rows_; ++m) {
          tmp[i][j] += (*this)[i][m] * other[m][j];
        }
      }
    }
    swap(tmp);
  }

  value_type determinant() const {
    if (rows_ != cols_) {
      throw std::invalid_argument("Determinant. Matrix is not square");
    }
    return det();
  }

  matrix transpose() const {
    matrix result(cols_, rows_);
    for (size_type i = 0; i < rows_; ++i) {
      for (size_type j = 0; j < cols_; ++j) {
        result[j][i] = (*this)[i][j];
      }
    }
    return result;
  }

  matrix complements() const {
    if (rows_ != cols_) {
      throw std::invalid_argument("Compliments: Matrix is not square");
    } else if (rows_ == 1) {
      throw std::invalid_argument("Compliment: Rows value should be greater than 1");
    }
    matrix result(rows_, cols_);
    for (size_type i = 0; i < rows_; ++i) {
      for (size_type j = 0; j < cols_; ++j) {
        result[i][j] = std::pow(-1.0, i + j) * minor(i, j).determinant();
      }
    }
    return result;
  }

  matrix inverse() const {
    matrix result(cols_, rows_);
    value_type det = determinant();
    if (is_equal(det, 0)) {
      throw std::invalid_argument("Inverse: Determinant of this matrix is equal 0");
    }
    if (cols_ == 1 && rows_ == 1) {
      result[0][0] = 1.0 / data_[0];
    } else {
      result = complements().transpose() * (1.0 / det);
    }
    return result;
  }

  bool operator==(const matrix& other) const noexcept {
    return equal(other);
  }

  bool operator!=(const matrix& other) const noexcept {
    return !equal(other);
  }

  matrix operator-(const matrix& other) const {
    matrix tmp(*this);
    tmp -= other;
    return tmp;
  }

  matrix operator+(const matrix& other) const {
    matrix tmp(*this);
    tmp += other;
    return tmp;
  }

  matrix operator*(const matrix& other) const {
    matrix tmp(*this);
    tmp *= other;
    return tmp;
  }

  matrix& operator+=(const matrix& other) {
    sum(other);
    return *this;
  }

  matrix& operator-=(const matrix& other) {
    sub(other);
    return *this;
  }

  matrix& operator*=(Tp rhs) {
    multiply(rhs);
    return *this;
  }

  matrix& operator*=(const matrix& other) {
    multiply(other);
    return *this;
  }

  friend matrix operator*(value_type lhs, const matrix& rhs) {
    matrix tmp(rhs);
    tmp *= lhs;
    return tmp;
  }

  friend matrix operator*(const matrix& lhs, value_type rhs) {
    matrix tmp(lhs);
    tmp *= rhs;
    return tmp;
  }

 protected:
  bool size_is_equal(const matrix& other) const noexcept {
    return rows_ == other.rows_ && cols_ == other.cols_;
  }

  matrix minor(size_type row, size_type col) const {
    matrix minor(rows_ - 1, cols_ - 1);
    size_type sub_i = 0, sub_j;
    for (size_type i = 0; i < rows_; ++i) {
      sub_j = 0;
      if (i != row) {
        for (size_type j = 0; j < cols_; ++j) {
          if (j != col) {
            minor[sub_i][sub_j++] = (*this)[i][j];
          }
        }
        sub_i++;
      }
    }
    return minor;
  }

  value_type det() const {
    value_type dtr;
    if (rows_ == 1) {
      dtr = data_[0];
    } else if (rows_ == 2) {
      dtr = data_[0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
    } else {
      int sign = 1;
      for (size_type i = 0; i < rows_;) {
        dtr = sign * data_[i++] * minor(0, i).det();
        sign = -sign;
      }
    }
    return dtr;
  }

 private:
  size_type rows_{};
  size_type cols_{};
  std::vector<Tp> data_;
};


} // namespace stl

#endif // STL_MATRIX_H_