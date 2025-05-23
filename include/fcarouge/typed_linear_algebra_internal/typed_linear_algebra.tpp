/* Typed Linear Algebra
Version 0.1.0
https://github.com/FrancoisCarouge/TypedLinearAlgebra

SPDX-License-Identifier: Unlicense

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org> */

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_TPP

namespace fcarouge {

namespace tla = typed_linear_algebra_internal;

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr bool
operator==(const typed_matrix<Matrix1, RowIndexes, ColumnIndexes> &lhs,
           const typed_matrix<Matrix2, RowIndexes, ColumnIndexes> &rhs) {
  return lhs.data == rhs.data;
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes, typename Indexes>
[[nodiscard]] inline constexpr auto
operator*(const typed_matrix<Matrix1, RowIndexes, Indexes> &lhs,
          const typed_matrix<Matrix2, Indexes, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<tla::product<Matrix1, Matrix2>>, RowIndexes,
                      ColumnIndexes>{lhs.data * rhs.data};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto operator*(Scalar lhs, const Matrix &rhs) {
  return tla::element<Matrix, 0, 0>{lhs * rhs.data(0)};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator*(Scalar lhs,
          const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<Matrix>, RowIndexes, ColumnIndexes>{
      lhs * rhs.data};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto
operator*(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return tla::element<Matrix, 0, 0>{lhs.data(0) * rhs};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator*(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return typed_matrix<tla::evaluate<Matrix>, RowIndexes, ColumnIndexes>{
      lhs.data * rhs};
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator+(const typed_matrix<Matrix1, RowIndexes, ColumnIndexes> &lhs,
          const typed_matrix<Matrix2, RowIndexes, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<Matrix1>, RowIndexes, ColumnIndexes>{
      lhs.data + rhs.data};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto operator+(const Matrix &lhs, Scalar rhs) {
  //! @todo Scalar will become Index with constraints.
  return tla::element<Matrix, 0, 0>{lhs.data(0) + rhs};
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator-(const typed_matrix<Matrix1, RowIndexes, ColumnIndexes> &lhs,
          const typed_matrix<Matrix2, RowIndexes, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<Matrix1>, RowIndexes, ColumnIndexes>{
      lhs.data - rhs.data};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto operator-(Scalar lhs, const Matrix &rhs) {
  return tla::element<Matrix, 0, 0>{lhs - rhs.data(0)};
}

template <typename Matrix1, typename Matrix2, typename RowIndexes1,
          typename RowIndexes2, typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator/(const typed_matrix<Matrix1, RowIndexes1, ColumnIndexes> &lhs,
          const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<tla::quotient<Matrix1, Matrix2>>,
                      RowIndexes1, RowIndexes2>{lhs.data / rhs.data};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator/(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return typed_matrix<tla::evaluate<Matrix>, RowIndexes, ColumnIndexes>{
      lhs.data / rhs};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto operator/(const Matrix &lhs, Scalar rhs) {
  return tla::element<Matrix, 0, 0>{lhs.data(0) / rhs};
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_TPP
