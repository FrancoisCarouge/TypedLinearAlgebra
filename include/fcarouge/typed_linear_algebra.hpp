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

#include "typed_linear_algebra_forward.hpp"
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
private:
  static_assert(tla::algebraic<Matrix>);

  //! @todo Privatize this section.
public:
  //! @name Private Member Functions
  //! @{

  //! @brief
  explicit inline constexpr typed_matrix(const Matrix &other);

  //! @}

private:
  //! @name Private Member Variables
  //! @{

  //! @brief Underlying algebraic backend data storage.
  Matrix matrix;

  //! @}

public:
  //! @name Public Member Types
  //! @{

  //! @brief The type of the element's underlying storage.
  using underlying = tla::underlying_t<Matrix>;

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

  //! @brief
  inline constexpr typed_matrix() = default;

  //! @brief
  inline constexpr typed_matrix(const typed_matrix &other) = default;

  //! @brief
  inline constexpr typed_matrix &operator=(const typed_matrix &other) = default;

  //! @brief
  inline constexpr typed_matrix(typed_matrix &&other) = default;

  //! @brief
  inline constexpr typed_matrix &operator=(typed_matrix &&other) = default;

  //! @brief
  template <tla::algebraic OtherMatrix>
  inline constexpr typed_matrix(
      const typed_matrix<OtherMatrix, RowIndexes, ColumnIndexes> &other);

  //! @brief
  inline constexpr explicit typed_matrix(const element<0, 0> (
      &elements)[tla::size<RowIndexes> * tla::size<ColumnIndexes>])
    requires tla::uniform<typed_matrix> && tla::one_dimension<typed_matrix>;

  //! @brief
  template <tla::arithmetic Type>
  explicit inline constexpr typed_matrix(const Type &value)
    requires tla::singleton<typed_matrix>;

  //! @brief
  template <typename Type>
  inline constexpr explicit typed_matrix(
      std::initializer_list<std::initializer_list<Type>> row_list)
    requires tla::uniform<typed_matrix>;

  //! @brief
  template <typename... Types>
  explicit inline constexpr typed_matrix(const Types &...values)
    requires tla::row<typed_matrix> && (not tla::column<typed_matrix>) &&
             tla::same_size<ColumnIndexes, std::tuple<Types...>>;

  //! @brief
  template <typename... Types>
  inline constexpr typed_matrix(const Types &...values)
    requires tla::column<typed_matrix> && (not tla::row<typed_matrix>) &&
             tla::same_size<RowIndexes, std::tuple<Types...>>;

  //! @brief
  [[nodiscard]] inline constexpr explicit(false) operator element<0, 0> &()
    requires tla::singleton<typed_matrix>;

  //! @brief
  [[nodiscard]] inline constexpr auto &&operator[](this auto &&self,
                                                   std::size_t index)
    requires tla::uniform<typed_matrix> && tla::one_dimension<typed_matrix>;

  //! @brief
  [[nodiscard]] inline constexpr auto &&
  operator[](this auto &&self, std::size_t row, std::size_t column)
    requires tla::uniform<typed_matrix>;

  //! @brief
  [[nodiscard]] inline constexpr auto &&operator()(this auto &&self,
                                                   std::size_t index)
    requires tla::uniform<typed_matrix> && tla::one_dimension<typed_matrix>;

  //! @brief
  [[nodiscard]] inline constexpr auto &&
  operator()(this auto &&self, std::size_t row, std::size_t column)
    requires tla::uniform<typed_matrix>;

  template <std::size_t Row, std::size_t Column>
    requires tla::in_range<Row, 0, tla::size<RowIndexes>> &&
             tla::in_range<Column, 0, tla::size<ColumnIndexes>>
  [[nodiscard]] inline constexpr element<Row, Column> &at() {
    return tla::element_traits<underlying, element<Row, Column>>::
        from_underlying(matrix(std::size_t{Row}, std::size_t{Column}));
  }

  //! @todo Can we deduplicate with deducing this?
  template <std::size_t Row, std::size_t Column>
    requires tla::in_range<Row, 0, tla::size<RowIndexes>> &&
             tla::in_range<Column, 0, tla::size<ColumnIndexes>>
  [[nodiscard]] inline constexpr element<Row, Column> at() const {
    return tla::element_traits<underlying, element<Row, Column>>::
        from_underlying(matrix(std::size_t{Row}, std::size_t{Column}));
  }

  template <std::size_t Index>
    requires tla::column<typed_matrix> &&
             tla::in_range<Index, 0, tla::size<RowIndexes>>
  [[nodiscard]] inline constexpr element<Index, 0> &at() {
    return tla::element_traits<underlying, element<Index, 0>>::from_underlying(
        matrix(std::size_t{Index}));
  }

  //! @todo Can we deduplicate with deducing this?
  template <std::size_t Index>
    requires tla::column<typed_matrix> &&
             tla::in_range<Index, 0, tla::size<RowIndexes>>
  [[nodiscard]] inline constexpr element<Index, 0> at() const {
    return tla::element_traits<underlying, element<Index, 0>>::from_underlying(
        matrix(std::size_t{Index}));
  }

  [[nodiscard]] inline constexpr auto &&data(this auto &&self) {
    return std::forward<decltype(self)>(self).matrix;
  }

  //! @}
};

//! @brief Strongly typed row vector.
template <typename Matrix, typename... ColumnIndexes>
using typed_row_vector =
    typed_matrix<Matrix, tla::identity_index, std::tuple<ColumnIndexes...>>;

//! @brief Strongly typed column vector.
template <typename Matrix, typename... RowIndexes>
using typed_column_vector =
    typed_matrix<Matrix, std::tuple<RowIndexes...>, tla::identity_index>;

//! @}

} // namespace fcarouge

#include "typed_linear_algebra_internal/typed_linear_algebra.tpp"

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_HPP
