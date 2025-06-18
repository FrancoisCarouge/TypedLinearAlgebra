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

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <tla::algebraic OtherMatrix>
inline constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const typed_matrix<OtherMatrix, RowIndexes, ColumnIndexes> &other)
    : matrix{other.data()} {}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
inline constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const Matrix &other)
    : matrix{other} {}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
inline constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const element<0, 0> (&elements)[typed_matrix::rows * typed_matrix::columns])
  requires tla::uniform<typed_matrix> && tla::one_dimension<typed_matrix>
    : matrix{elements} {}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <tla::arithmetic Type>
inline constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const Type &value)
  requires tla::singleton<typed_matrix>
{
  data()(0, 0) = cast<underlying, Type>(value);
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename Type>
inline constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    std::initializer_list<std::initializer_list<Type>> row_list)
  requires tla::uniform<typed_matrix>
{
  for (std::size_t i{0}; const auto &row : row_list) {
    for (std::size_t j{0}; const auto &value : row) {
      data()(i, j) = cast<underlying, Type>(value);
      ++j;
    }
    ++i;
  }
}

//! @todo Combine the two constructors in ome?
//! @todo Verify if the types are the same, or assignable, for nicer error?
//! @todo Rewrite with a fold expression over the pack?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename... Types>
inline constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const Types &...values)
  requires tla::row<typed_matrix> && (not tla::column<typed_matrix>) &&
           tla::same_size<ColumnIndexes, std::tuple<Types...>>
{
  std::tuple value_pack{values...};
  tla::for_constexpr<0, typed_matrix::columns, 1>(
      [this, &value_pack](auto position) {
        auto value{std::get<position>(value_pack)};
        using type = std::remove_cvref_t<decltype(value)>;
        data()[position] = cast<underlying, type>(value);
      });
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename... Types>
inline constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const Types &...values)
  requires tla::column<typed_matrix> && (not tla::row<typed_matrix>) &&
           tla::same_size<RowIndexes, std::tuple<Types...>>
{
  std::tuple value_pack{values...};
  tla::for_constexpr<0, typed_matrix::rows, 1>(
      [this, &value_pack](auto position) {
        auto value{std::get<position>(value_pack)};
        using type = std::remove_cvref_t<decltype(value)>;
        data()[position] = cast<underlying, type>(value);
      });
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] inline constexpr typed_matrix<
    Matrix, RowIndexes, ColumnIndexes>::operator element<0, 0> &()
  requires tla::singleton<typed_matrix>
{
  return cast<element<0, 0> &, underlying &>(data()(0, 0));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] inline constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator[](this auto &&self,
                                                            std::size_t index)
  requires(tla::uniform<typed_matrix> && tla::one_dimension<typed_matrix>)
{
  return std::forward<decltype(self)>(self).data()(index);
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] inline constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator[](this auto &&self,
                                                            std::size_t row,
                                                            std::size_t column)
  requires tla::uniform<typed_matrix>
{
  return std::forward<decltype(self)>(self).data()(row, column);
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] inline constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator()(this auto &&self,
                                                            std::size_t index)
  requires tla::uniform<typed_matrix> && tla::one_dimension<typed_matrix>
{
  return std::forward<decltype(self)>(self).data()(index);
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] inline constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator()(this auto &&self,
                                                            std::size_t row,
                                                            std::size_t column)
  requires tla::uniform<typed_matrix>
{
  return std::forward<decltype(self)>(self).data()(row, column);
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr bool
operator==(const typed_matrix<Matrix1, RowIndexes, ColumnIndexes> &lhs,
           const typed_matrix<Matrix2, RowIndexes, ColumnIndexes> &rhs) {
  return lhs.data() == rhs.data();
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes, typename Indexes>
[[nodiscard]] inline constexpr auto
operator*(const typed_matrix<Matrix1, RowIndexes, Indexes> &lhs,
          const typed_matrix<Matrix2, Indexes, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<tla::product<Matrix1, Matrix2>>, RowIndexes,
                      ColumnIndexes>{lhs.data() * rhs.data()};
}

template <tla::arithmetic Scalar, typename Matrix>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto operator*(Scalar lhs, const Matrix &rhs) {
  return tla::element<Matrix, 0, 0>{lhs * rhs.data()(0)};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator*(Scalar lhs,
          const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<Matrix>, RowIndexes, ColumnIndexes>{
      lhs * rhs.data()};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto
operator*(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return tla::element<Matrix, 0, 0>{lhs.data()(0) * rhs};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator*(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return typed_matrix<tla::evaluate<Matrix>, RowIndexes, ColumnIndexes>{
      lhs.data() * rhs};
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator+(const typed_matrix<Matrix1, RowIndexes, ColumnIndexes> &lhs,
          const typed_matrix<Matrix2, RowIndexes, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<Matrix1>, RowIndexes, ColumnIndexes>{
      lhs.data() + rhs.data()};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto operator+(const Matrix &lhs, Scalar rhs) {
  //! @todo Scalar will become Index with constraints.
  return tla::element<Matrix, 0, 0>{lhs.data()(0) + rhs};
}

template <typename Matrix1, typename Matrix2, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator-(const typed_matrix<Matrix1, RowIndexes, ColumnIndexes> &lhs,
          const typed_matrix<Matrix2, RowIndexes, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<Matrix1>, RowIndexes, ColumnIndexes>{
      lhs.data() - rhs.data()};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto operator-(Scalar lhs, const Matrix &rhs) {
  return tla::element<Matrix, 0, 0>{lhs - rhs.data()(0)};
}

template <typename Matrix1, typename Matrix2, typename RowIndexes1,
          typename RowIndexes2, typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator/(const typed_matrix<Matrix1, RowIndexes1, ColumnIndexes> &lhs,
          const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes> &rhs) {
  return typed_matrix<tla::evaluate<tla::quotient<Matrix1, Matrix2>>,
                      RowIndexes1, RowIndexes2>{lhs.data() / rhs.data()};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
[[nodiscard]] inline constexpr auto
operator/(const typed_matrix<Matrix, RowIndexes, ColumnIndexes> &lhs,
          Scalar rhs) {
  return typed_matrix<tla::evaluate<Matrix>, RowIndexes, ColumnIndexes>{
      lhs.data() / rhs};
}

template <tla::arithmetic Scalar, typename Matrix, typename RowIndexes,
          typename ColumnIndexes>
  requires tla::singleton<Matrix>
[[nodiscard]] inline constexpr auto operator/(const Matrix &lhs, Scalar rhs) {
  return tla::element<Matrix, 0, 0>{lhs.data()(0) / rhs};
}

template <tla::arithmetic To, tla::arithmetic From>
struct element_caster<To, From> {
  [[nodiscard]] inline constexpr To operator()(const From &value) const {
    return value;
  }
};

template <tla::arithmetic To, tla::arithmetic From>
struct element_caster<To &, From &> {
  [[nodiscard]] inline constexpr To &operator()(From &value) const {
    return value;
  }
};
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_TPP
