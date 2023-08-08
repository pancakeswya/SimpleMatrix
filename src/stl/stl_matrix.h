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
    : m_rows(list.size()), m_cols(list.begin()->size()) {
    for (auto& i : list) {
      assert(i.size() == m_cols);
      for (auto& j : i) {
        m_data.push_back(j);
      }
    }
  }

  matrix(size_type rows, size_type cols)
    : m_rows(rows), m_cols(cols), m_data(rows * cols) {}

  matrix(const matrix& other)
    : m_rows(other.m_rows), m_cols(other.m_cols), m_data(other.m_data) {}

  matrix(matrix&& other) noexcept
    : m_rows(other.m_rows), m_cols(other.m_cols), m_data(std::move(other.m_data)) {
    other.m_rows = other.m_cols = 0;
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
    std::swap(m_rows, other.m_rows);
    std::swap(m_cols, other.m_cols);
    std::swap(m_data, other.m_data);
  }

  iterator begin() noexcept {
    return m_data.begin();
  }

  const_iterator begin() const noexcept {
    return m_data.begin();
  }

  iterator end() noexcept {
    return m_data.end();
  }

  const_iterator end() const noexcept {
    return m_data.end();
  }

  void clear() noexcept {
    m_data.clear();
    m_cols = m_rows = 0;
  }

  pointer operator[](size_type row) noexcept {
    return &m_data[row * m_cols];
  }

  const_pointer operator[](size_type row) const noexcept {
    return &m_data[row * m_cols];
  }

  reference at(size_type row, size_type col) {
    if (row > m_rows || col > m_cols) {
      throw std::out_of_range("operator[]. Row or col is out of range");
    }
    return m_data[row * m_cols + col];
  }

  const_reference at(size_type row, size_type col) const {
    if (row > m_rows || col > m_cols) {
      throw std::out_of_range("operator[]. Row or col is out of range");
    }
    return m_data[row * m_cols + col];
  }

  size_type rows() const noexcept {
    return m_rows;
  }

  size_type cols() const noexcept {
    return m_cols;
  }

  bool equal(const matrix& other) const noexcept {
    if (!size_is_equal(other)) {
      return false;
    }
    for (size_type i = 0; i < m_data.size(); ++i) {
      if (!is_equal(m_data[i], other.m_data[i])) {
        return false;
      }
    }
    return true;
  }

  void sum(const matrix& other) {
    if (!size_is_equal(other)) {
      throw std::invalid_argument("Sum. Different sizes of matrices");
    }
    for (size_type i = 0; i < m_data.size(); ++i) {
        m_data[i] += other.m_data[i];
    }
  }

  void sub(const matrix& other) {
    if (!size_is_equal(other)) {
      throw std::invalid_argument("Sub. Different sizes of matrices");
    }
    for (int i = 0; i < m_data.size(); ++i) {
        m_data[i] -= other.m_data[i];
    }
  }

  void multiply(const Tp& value) noexcept {
    for (int i = 0; i < m_data.size(); ++i) {
        m_data[i] *= value;
    }
  }

  void multiply(const matrix& other) {
    if (m_cols != other.m_rows || m_rows != other.m_cols) {
      throw std::invalid_argument(
          "Multiply. First matrix cols != Second matrix rows or vice versa");
    }
    matrix tmp(m_rows, other.m_cols);
    for (size_type i = 0; i < m_rows; ++i) {
      for (size_type j = 0; j < other.m_cols; ++j) {
        for (size_type m = 0; m < other.m_rows; ++m) {
          tmp[i][j] += (*this)[i][m] * other[m][j];
        }
      }
    }
    swap(tmp);
  }

  value_type determinant() const {
    if (m_rows != m_cols) {
      throw std::invalid_argument("Determinant. Matrix is not square");
    }
    return det();
  }

  matrix transpose() const {
    matrix result(m_cols, m_rows);
    for (size_type i = 0; i < m_rows; ++i) {
      for (size_type j = 0; j < m_cols; ++j) {
        result[j][i] = (*this)[i][j];
      }
    }
    return result;
  }

  matrix complements() const {
    if (m_rows != m_cols) {
      throw std::invalid_argument("Compliments: Matrix is not square");
    } else if (m_rows == 1) {
      throw std::invalid_argument("Compliment: Rows value should be greater than 1");
    }
    matrix result(m_rows, m_cols);
    for (size_type i = 0; i < m_rows; ++i) {
      for (size_type j = 0; j < m_cols; ++j) {
        result[i][j] = std::pow(-1.0, i + j) * minor(i, j).determinant();
      }
    }
    return result;
  }

  matrix inverse() const {
    matrix result(m_cols, m_rows);
    value_type det = determinant();
    if (is_equal(det, 0)) {
      throw std::invalid_argument("Inverse: Determinant of this matrix is equal 0");
    }
    if (m_cols == 1 && m_rows == 1) {
      result[0][0] = 1.0 / m_data[0];
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
    return m_rows == other.m_rows && m_cols == other.m_cols;
  }

  matrix minor(size_type row, size_type col) const {
    matrix minor(m_rows - 1, m_cols - 1);
    size_type sub_i = 0, sub_j;
    for (size_type i = 0; i < m_rows; ++i) {
      sub_j = 0;
      if (i != row) {
        for (size_type j = 0; j < m_cols; ++j) {
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
    if (m_rows == 1) {
      dtr = m_data[0];
    } else if (m_rows == 2) {
      dtr = m_data[0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
    } else {
      int sign = 1;
      for (size_type i = 0; i < m_rows;) {
        dtr = sign * m_data[i++] * minor(0, i).det();
        sign = -sign;
      }
    }
    return dtr;
  }

 private:
  size_type m_rows{};
  size_type m_cols{};
  std::vector<Tp> m_data;
};


} // namespace stl

#endif // STL_MATRIX_H_