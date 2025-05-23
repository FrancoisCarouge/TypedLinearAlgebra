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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_HPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_HPP

//! @file
//! @brief Typed linear algebra implementation.
//!
//! @details Typed matrix, vectors, and operations.

#include "typed_linear_algebra_internal/format.hpp"
#include "typed_linear_algebra_internal/utility.hpp"

#include <concepts>
#include <cstddef>
#include <format>
#include <initializer_list>
#include <tuple>

namespace fcarouge {

namespace tla = typed_linear_algebra_internal;

//! @name Types
//! @{

//! @brief Strongly typed matrix.
//!
//! @details Compose a linear algebra backend matrix into a typed matrix. Row
//! and column indexes provide each element's index type.
//!
//! @tparam Matrix The underlying linear algebra matrix.
//! @tparam RowIndexes The packed types of the row indexes.
//! @tparam ColumnIndexes The packed types of the column indexes.
//!
//! @note Type safety cannot be guaranteed at compilation time without index
//! safety. The indexes can either be non-type template parameters or strong
//! types overloadings. Converting a runtime index to a dependent template type
//! is not possible in C++. A proxy reference could be used to allow traditional
//! assignment syntax but the runtime check and extra indirection are not
//! interesting tradeoffs. A template call operator can be used for getting a
//! type safe value but impractical syntax for setting. Without index safety,
//! the accepted tradeoff is a templated index `at<i, j>()` method.
//!
//! @note Deduction guides are tricky because a given element type comes from
//! a row and column index to be deduced.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
struct typed_matrix {
  static_assert(tla::algebraic<Matrix>);
  //! @todo Privatize this section.
public:
  //! @name Private Member Types
  //! @{

  //! @brief The type of the element's underlying storage.
  using underlying = tla::underlying_t<Matrix>;

  //! @}

  //! @name Private Member Functions
  //! @{

  //! @todo Can this be removed altogether?
  explicit inline constexpr typed_matrix(const Matrix &other) : data{other} {}

  //! @}

  //! @name Private Member Variables
  //! @{

  Matrix data;

  //! @}

public:
  //! @name Public Member Types
  //! @{

  //! @brief The tuple with the row components of the indexes.
  using row_indexes = RowIndexes;

  //! @brief The tuple with the column components of the indexes.
  using column_indexes = ColumnIndexes;

  //! @brief The type of the element at the given matrix indexes position.
  template <std::size_t RowIndex, std::size_t ColumnIndex>
  using element = tla::element<typed_matrix, RowIndex, ColumnIndex>;

  //! @}

  //! @name Public Member Variables
  //! @{

  //! @brief The count of rows.
  inline constexpr static std::size_t rows{tla::size<row_indexes>};

  //! @brief The count of rows.
  inline constexpr static std::size_t columns{tla::size<column_indexes>};

  //! @}

  //! @name Public Member Functions
  //! @{

  inline constexpr typed_matrix() = default;

  inline constexpr typed_matrix(const typed_matrix &other) = default;

  inline constexpr typed_matrix &operator=(const typed_matrix &other) = default;

  inline constexpr typed_matrix(typed_matrix &&other) = default;

  inline constexpr typed_matrix &operator=(typed_matrix &&other) = default;

  //! @todo Requires evaluated types of Matrix and OtherMatrix are identical?
  template <tla::algebraic OtherMatrix>
  inline constexpr typed_matrix(
      const typed_matrix<OtherMatrix, RowIndexes, ColumnIndexes> &other)
      : data{other.data} {}

  inline constexpr explicit typed_matrix(const element<0, 0> (
      &elements)[tla::size<RowIndexes> * tla::size<ColumnIndexes>])
    requires tla::uniform<typed_matrix> && tla::one_dimension<typed_matrix>
      : data{elements} {}

  template <tla::arithmetic Type>
    requires tla::singleton<typed_matrix>
  explicit inline constexpr typed_matrix(const Type &value) {
    data(0, 0) = tla::element_traits<underlying, Type>::to_underlying(value);
  }

  //! @todo Verify the list sizes at runtime?
  template <typename Type>
  inline constexpr explicit typed_matrix(
      std::initializer_list<std::initializer_list<Type>> row_list)
    requires tla::uniform<typed_matrix>
  {
    for (std::size_t i{0}; const auto &row : row_list) {
      for (std::size_t j{0}; const auto &value : row) {
        data(i, j) =
            tla::element_traits<underlying, Type>::to_underlying(value);
        ++j;
      }
      ++i;
    }
  }

  //! @todo Combine the two constructors in ome?
  //! @todo Verify if the types are the same, or assignable, for nicer error?
  //! @todo Rewrite with a fold expression over the pack?
  template <typename... Types>
    requires tla::row<typed_matrix> && (not tla::column<typed_matrix>) &&
             tla::same_size<ColumnIndexes, std::tuple<Types...>>
  explicit inline constexpr typed_matrix(const Types &...values) {
    std::tuple value_pack{values...};
    tla::for_constexpr<0, tla::size<ColumnIndexes>, 1>(
        [this, &value_pack](auto position) {
          auto value{std::get<position>(value_pack)};
          using type = std::remove_cvref_t<decltype(value)>;
          data[position] =
              tla::element_traits<underlying, type>::to_underlying(value);
        });
  }

  template <typename... Types>
    requires tla::column<typed_matrix> && (not tla::row<typed_matrix>) &&
             tla::same_size<RowIndexes, std::tuple<Types...>>
  inline constexpr typed_matrix(const Types &...values) {
    std::tuple value_pack{values...};
    tla::for_constexpr<0, tla::size<RowIndexes>, 1>(
        [this, &value_pack](auto position) {
          auto value{std::get<position>(value_pack)};
          using type = std::remove_cvref_t<decltype(value)>;
          data[position] =
              tla::element_traits<underlying, type>::to_underlying(value);
        });
  }

  [[nodiscard]] inline constexpr explicit(false) operator element<0, 0> &()
    requires tla::singleton<typed_matrix>
  {
    return tla::element_traits<underlying, element<0, 0>>::from_underlying(
        data(0, 0));
  }

  [[nodiscard]] inline constexpr auto &&operator[](this auto &&self,
                                                   std::size_t index)
    requires tla::uniform<typed_matrix> && tla::one_dimension<typed_matrix>
  {
    return std::forward<decltype(self)>(self).data(index);
  }

  [[nodiscard]] inline constexpr auto &&
  operator[](this auto &&self, std::size_t row, std::size_t column)
    requires tla::uniform<typed_matrix>
  {
    return std::forward<decltype(self)>(self).data(row, column);
  }

  [[nodiscard]] inline constexpr auto &&operator()(this auto &&self,
                                                   std::size_t index)
    requires tla::uniform<typed_matrix> && tla::one_dimension<typed_matrix>
  {
    return std::forward<decltype(self)>(self).data(index);
  }

  [[nodiscard]] inline constexpr auto &&
  operator()(this auto &&self, std::size_t row, std::size_t column)
    requires tla::uniform<typed_matrix>
  {
    return std::forward<decltype(self)>(self).data(row, column);
  }

  template <std::size_t Row, std::size_t Column>
    requires tla::in_range<Row, 0, tla::size<RowIndexes>> &&
             tla::in_range<Column, 0, tla::size<ColumnIndexes>>
  [[nodiscard]] inline constexpr element<Row, Column> &at() {
    return tla::element_traits<underlying, element<Row, Column>>::
        from_underlying(data(std::size_t{Row}, std::size_t{Column}));
  }

  template <std::size_t Index>
    requires tla::column<typed_matrix> &&
             tla::in_range<Index, 0, tla::size<RowIndexes>>
  [[nodiscard]] inline constexpr element<Index, 0> &at() {
    return tla::element_traits<underlying, element<Index, 0>>::from_underlying(
        data(std::size_t{Index}));
  }

  //! @}
};

//! @brief Strongly typed row vector.
template <typename Matrix, typename... ColumnIndexes>
using typed_row_vector =
    typed_matrix<Matrix, std::tuple<int>, std::tuple<ColumnIndexes...>>;

//! @brief Strongly typed column vector.
template <typename Matrix, typename... RowIndexes>
using typed_column_vector =
    typed_matrix<Matrix, std::tuple<RowIndexes...>, std::tuple<int>>;

//! @}

} // namespace fcarouge

#include "typed_linear_algebra_internal/typed_linear_algebra.tpp"

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_HPP
