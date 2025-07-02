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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TYPED_LINEAR_ALGEBRA_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TYPED_LINEAR_ALGEBRA_TPP

namespace fcarouge {
namespace tla = typed_linear_algebra_internal;

//! @todo Replace all calls to data() by direct matrix access?

//! @todo Verify types and storage (?) compatibility.
//! @todo Also add the equivalent operator=.
//! @todo Combine with the similar delcarations.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename Matrix2, typename RowIndexes2, typename ColumnIndexes2>
inline constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2> &other)
    : matrix{other.data()} {}

//! @todo Should the arithmetic constraint be dropped? The parameter renamed
//! to element systematically? A requirement of compatible conversion?
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

//! @todo Verify the list sizes at runtime? Deprecate?
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
[[nodiscard]] inline constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::
operator element<0, 0> &&(this auto &&self)
  requires tla::singleton<typed_matrix>
{
  // This is a form of `std::forward_like`, is there a simpler, or more compact
  // syntax?
  constexpr bool is_adding_const{
      std::is_const_v<std::remove_reference_t<decltype(self)>>};
  if constexpr (std::is_lvalue_reference_v<decltype(self) &&>) {
    if constexpr (is_adding_const)
      return cast<element<0, 0>, underlying>(
          std::forward<decltype(self)>(self).data()(0, 0));
    else
      return cast<element<0, 0> &, underlying &>(
          std::forward<decltype(self)>(self).data()(0, 0));
  } else {
    if constexpr (is_adding_const)
      return std::move(cast<element<0, 0>, underlying>(
          std::forward<decltype(self)>(self).data()(0, 0)));
    else
      return std::move(cast<element<0, 0> &&, underlying &&>(
          std::forward<decltype(self)>(self).data()(0, 0)));
  }
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

//! @todo Can we deduplicate with deducing this?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Row, std::size_t Column>
[[nodiscard]] inline constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at()
    -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, Row,
                    Column> &
  requires tla::in_range<
               Row, 0, typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows> &&
           tla::in_range<
               Column, 0,
               typed_matrix<Matrix, RowIndexes, ColumnIndexes>::columns>
{
  return cast<element<Row, Column> &, underlying &>(
      matrix(std::size_t{Row}, std::size_t{Column}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Row, std::size_t Column>
[[nodiscard]] inline constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at() const
    -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, Row,
                    Column>
  requires tla::in_range<
               Row, 0, typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows> &&
           tla::in_range<
               Column, 0,
               typed_matrix<Matrix, RowIndexes, ColumnIndexes>::columns>
{
  return cast<element<Row, Column>, underlying>(
      matrix(std::size_t{Row}, std::size_t{Column}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Index>
[[nodiscard]] inline constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at()
    -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, Index, 0> &
  requires tla::column<typed_matrix<Matrix, RowIndexes, ColumnIndexes>> &&
           tla::in_range<Index, 0,
                         typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows>
{
  return cast<element<Index, 0> &, underlying &>(matrix(std::size_t{Index}));
}

//! @todo Add row-vector overload.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Index>
[[nodiscard]] inline constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at() const
    -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, Index, 0>
  requires tla::column<typed_matrix<Matrix, RowIndexes, ColumnIndexes>> &&
           tla::in_range<Index, 0,
                         typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows>
{
  return cast<element<Index, 0>, underlying>(matrix(std::size_t{Index}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] inline constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::data(this auto &&self) {
  return std::forward<decltype(self)>(self).matrix;
}

} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TYPED_LINEAR_ALGEBRA_TPP
