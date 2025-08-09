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

//! @todo Verify types and storage (?) compatibility.
//! @todo Also add the equivalent operator=.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename Matrix2, typename RowIndexes2, typename ColumnIndexes2>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2> &other)
    : storage{other.data()} {}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const Matrix &other)
    : storage{other} {}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const element<0, 0> (&elements)[typed_matrix::rows * typed_matrix::columns])
  requires is_uniform_typed_matrix<typed_matrix> and
           is_one_dimension_typed_matrix<typed_matrix>
    : storage{elements} {}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const auto &value)
  requires is_singleton_typed_matrix<typed_matrix>
{
  using type = std::remove_cvref_t<decltype(value)>;
  storage(std::size_t{0}, std::size_t{0}) = cast<underlying, type>(value);
}

//! @todo Verify the list sizes at runtime? Deprecate?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename Type>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    std::initializer_list<std::initializer_list<Type>> row_list)
  requires is_uniform_typed_matrix<typed_matrix>
{
  for (std::size_t i{0}; const auto &row : row_list) {
    for (std::size_t j{0}; const auto &value : row) {
      storage(std::size_t{i}, std::size_t{j}) = cast<underlying, Type>(value);
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
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const Types &...values)
  requires is_row_typed_matrix<typed_matrix> and
           (not is_column_typed_matrix<typed_matrix>) and
           tla::same_size<ColumnIndexes, std::tuple<Types...>>
{
  std::tuple value_pack{values...};
  tla::for_constexpr<0, typed_matrix::columns, 1>(
      [this, &value_pack](auto position) {
        auto value{std::get<position>(value_pack)};
        using type = std::remove_cvref_t<decltype(value)>;
        storage(std::size_t{position}) = cast<underlying, type>(value);
      });
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename... Types>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const Types &...values)
  requires is_column_typed_matrix<typed_matrix> and
           (not is_row_typed_matrix<typed_matrix>) and
           tla::same_size<RowIndexes, std::tuple<Types...>>
{
  std::tuple value_pack{values...};
  tla::for_constexpr<0, typed_matrix::rows, 1>(
      [this, &value_pack](auto position) {
        auto value{std::get<position>(value_pack)};
        using type = std::remove_cvref_t<decltype(value)>;
        storage(std::size_t{position}) = cast<underlying, type>(value);
      });
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::
operator element<0, 0> &&(this auto &&self)
  requires is_singleton_typed_matrix<typed_matrix>
{
  // This is a form of `std::forward_like`, is there a simpler, or more compact
  // syntax?
  constexpr bool is_adding_const{
      std::is_const_v<std::remove_reference_t<decltype(self)>>};
  if constexpr (std::is_lvalue_reference_v<decltype(self) &&>) {
    if constexpr (is_adding_const)
      return cast<element<0, 0>, underlying>(
          std::forward<decltype(self)>(self).storage(std::size_t{0},
                                                     std::size_t{0}));
    else
      return cast<element<0, 0> &, underlying &>(
          std::forward<decltype(self)>(self).storage(std::size_t{0},
                                                     std::size_t{0}));
  } else {
    if constexpr (is_adding_const)
      return std::move(cast<element<0, 0>, underlying>(
          std::forward<decltype(self)>(self).storage(std::size_t{0},
                                                     std::size_t{0})));
    else
      return std::move(cast<element<0, 0> &&, underlying &&>(
          std::forward<decltype(self)>(self).storage(std::size_t{0},
                                                     std::size_t{0})));
  }
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator[](this auto &&self,
                                                            std::size_t index)
  requires(is_uniform_typed_matrix<typed_matrix> and
           is_one_dimension_typed_matrix<typed_matrix>)
{
  return std::forward<decltype(self)>(self).storage(std::size_t{index});
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator[](this auto &&self,
                                                            std::size_t row,
                                                            std::size_t column)
  requires is_uniform_typed_matrix<typed_matrix>
{
  return std::forward<decltype(self)>(self).storage(std::size_t{row},
                                                    std::size_t{column});
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator()(this auto &&self,
                                                            std::size_t index)
  requires is_uniform_typed_matrix<typed_matrix> and
           is_one_dimension_typed_matrix<typed_matrix>
{
  return std::forward<decltype(self)>(self).storage(std::size_t{index});
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator()(this auto &&self,
                                                            std::size_t row,
                                                            std::size_t column)
  requires is_uniform_typed_matrix<typed_matrix>
{
  return std::forward<decltype(self)>(self).storage(std::size_t{row},
                                                    std::size_t{column});
}

//! @todo Can we deduplicate with deducing this?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Row, std::size_t Column>
[[nodiscard]] constexpr auto
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
      storage(std::size_t{Row}, std::size_t{Column}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Row, std::size_t Column>
[[nodiscard]] constexpr auto
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
      storage(std::size_t{Row}, std::size_t{Column}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Index>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at()
    -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, Index, 0> &
  requires is_column_typed_matrix<
               typed_matrix<Matrix, RowIndexes, ColumnIndexes>> &&
           tla::in_range<Index, 0,
                         typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows>
{
  return cast<element<Index, 0> &, underlying &>(storage(std::size_t{Index}));
}

//! @todo Add row-vector overload.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Index>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at() const
    -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, Index, 0>
  requires is_column_typed_matrix<
               typed_matrix<Matrix, RowIndexes, ColumnIndexes>> &&
           tla::in_range<Index, 0,
                         typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows>
{
  return cast<element<Index, 0>, underlying>(storage(std::size_t{Index}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::data(this auto &&self) {
  return std::forward<decltype(self)>(self).storage;
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TYPED_LINEAR_ALGEBRA_TPP
