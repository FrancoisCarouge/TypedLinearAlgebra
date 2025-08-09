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
#include <utility>

namespace fcarouge {
namespace tla = typed_linear_algebra_internal;

//! @name Concepts
//! @{

//! @brief Concept of a typed matrix type.
template <typename Type>
concept is_typed_matrix = tla::is_typed_matrix<Type>;

//! @brief Concept of a singleton, one-element typed matrix type.
template <typename Type>
concept is_singleton_typed_matrix = tla::is_singleton_typed_matrix<Type>;

//! @brief Concept of a typed matrix in which all element types are the same.
//!
//! @details Matrices with uniform types are type safe even with the traditional
//! access operators.
//!
//! @note A matrix may be uniform with different row and column indexes.
template <typename Type>
concept is_uniform_typed_matrix = tla::is_uniform_typed_matrix<Type>;

//! @brief Concept of a typed matrix with only one dimension, vector, or column.
template <typename Type>
concept is_one_dimension_typed_matrix =
    tla::is_one_dimension_typed_matrix<Type>;

//! @brief Concept of a row typed matrix, vector.
template <typename Type>
concept is_row_typed_matrix = tla::is_row_typed_matrix<Type>;

//! @brief Concept of a column typed matrix, vector.
template <typename Type>
concept is_column_typed_matrix = tla::is_column_typed_matrix<Type>;

//! @}

//! @name Types
//! @{

//! @brief Strongly typed matrix.
//!
//! @details Compose a linear algebra backend matrix into a typed matrix. Row
//! and column indexes provide each element's index type.
//!
//! @tparam Matrix The underlying linear algebra matrix.
//! @tparam RowIndexes The tuple type of the row indexes.
//! @tparam ColumnIndexes The tuple type of the column indexes.
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
//!
//! @todo Don't limit the dimension to two? Use parameter pack of index tuples
//! for tensor types.
//! @todo Add complexity documentation where appropriate for the API?
//! @todo Document the lack of bound checking on the classical API or remove
//! altogether?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
class typed_matrix {
public:
  //! @name Public Member Types
  //! @{

  //! @brief The type of the composed, stored matrix.
  using matrix = Matrix;

  //! @brief The tuple with the row components of the indexes.
  using row_indexes = RowIndexes;

  //! @brief The tuple with the column components of the indexes.
  using column_indexes = ColumnIndexes;

  //! @brief The type of the element's underlying storage.
  using underlying = tla::underlying_t<Matrix>;

  //! @brief The type of the element at the given matrix indexes position.
  template <std::size_t RowIndex, std::size_t ColumnIndex>
  using element = tla::element<typed_matrix, RowIndex, ColumnIndex>;

  //! @}

  //! @name Public Member Variables
  //! @{

  //! @brief The count of rows.
  static inline constexpr std::size_t rows{std::tuple_size_v<row_indexes>};

  //! @brief The count of rows.
  static inline constexpr std::size_t columns{
      std::tuple_size_v<column_indexes>};

  //! @}

  //! @name Public Member Functions
  //! @{

  //! @brief Destruct a default typed matrix.
  constexpr ~typed_matrix() = default;

  //! @brief Construct a default typed matrix.
  //!
  //! @warning The initialization of the underlying matrix's storage follows the
  //! initialization behavior of the underlying matrix's type, which for some
  //! type means no initialization.
  constexpr typed_matrix() = default;

  //! @brief Copy construct a typed matrix.
  constexpr typed_matrix(const typed_matrix &other) = default;

  //! @brief Copy assign a typed matrix.
  constexpr typed_matrix &operator=(const typed_matrix &other) = default;

  //! @brief Move construct a typed matrix.
  constexpr typed_matrix(typed_matrix &&other) = default;

  //! @brief Move construct a typed matrix.
  constexpr typed_matrix &operator=(typed_matrix &&other) = default;

  //! @brief Copy construct the typed matrix.
  template <typename Matrix2, typename RowIndexes2, typename ColumnIndexes2>
  constexpr typed_matrix(
      const typed_matrix<Matrix2, RowIndexes2, ColumnIndexes2> &other);

  //! @brief Convert construct a typed matrix from an underlying matrix.
  //!
  //! @warning Useful for operations implementation where underlying data
  //! constrution is needed. Not recommended for convenience construction due to
  //! absence of type validation.
  constexpr explicit typed_matrix(const Matrix &other);

  //! @brief Convert construct a one-dimension uniformly typed matrix from
  //! array.
  //!
  //! @details Applicable to one-dimension matrix: column- or row-vector.
  //! Applicable to single-type matrix: uniform type of all elements.
  //! Single-argument constructors taking arrays of data should get implicit
  //! constructors.
  //!
  //! @param elements C-style array of elements of identical types.
  constexpr explicit typed_matrix(
      const element<0, 0> (&elements)[rows * columns])
    requires is_uniform_typed_matrix<typed_matrix> and
             is_one_dimension_typed_matrix<typed_matrix>;

  //! @brief Convert construct a singleton typed matrix from a single value.
  //!
  //! @details Applicable to singleton matrix: one element.
  //!
  //! @param value Element of compatible type.
  constexpr explicit typed_matrix(const auto &value)
    requires is_singleton_typed_matrix<typed_matrix>;

  //! @brief Convert construct a uniformly typed matrix from list-initializers.
  //!
  //! @details Applicable to matrix of uniform elements type. Single-argument
  //! constructors taking arrays of data should get implicit constructors.
  //!
  //! @param row_list List-initializers of list-initializer of elements.
  template <typename Type>
  constexpr explicit typed_matrix(
      std::initializer_list<std::initializer_list<Type>> row_list)
    requires is_uniform_typed_matrix<typed_matrix>;

  //! @brief Convert construct a row typed vector from elements.
  //!
  //! @details Applicable to one-dimension matrix: row-vector.
  //!
  //! @param values Parameter pack of elements.
  template <typename... Types>
  constexpr explicit typed_matrix(const Types &...values)
    requires is_row_typed_matrix<typed_matrix> and
             (not is_column_typed_matrix<typed_matrix>) and
             tla::same_size<ColumnIndexes, std::tuple<Types...>>;

  //! @brief Convert construct a column typed vector from elements.
  //!
  //! @details Applicable to one-dimension matrix: column-vector.
  //!
  //! @param values Parameter pack of elements.
  template <typename... Types>
  constexpr typed_matrix(const Types &...values)
    requires is_column_typed_matrix<typed_matrix> and
             (not is_row_typed_matrix<typed_matrix>) and
             tla::same_size<RowIndexes, std::tuple<Types...>>;

  //! @brief Access the singleton typed matrix element.
  //!
  //! @details Applicable to singleton matrix: one element. Returns a reference
  //! to the unique element of the typed matrix.
  [[nodiscard]] constexpr explicit(false)
  operator element<0, 0> &&(this auto &&self)
    requires is_singleton_typed_matrix<typed_matrix>;

  //! @brief Access the specified element.
  //!
  //! @details Applicable to one-dimension matrix: column- or row-vector.
  //! Applicable to single-type matrix: uniform type of all elements.
  //! Returns a reference to the element at the specified location.
  //!
  //! @param self Explicit object parameter deducing this: not user specified.
  //! @param index Position of the element to return.
  [[nodiscard]] constexpr auto &&operator[](this auto &&self, std::size_t index)
    requires(is_uniform_typed_matrix<typed_matrix> and
             is_one_dimension_typed_matrix<typed_matrix>);

  //! @brief Access the specified element.
  //!
  //! @details Applicable to single-type matrix: uniform type of all elements.
  //! Returns a reference to the element at the specified location.
  //!
  //! @param self Explicit object parameter deducing this: not user specified.
  //! @param row Row index of the element to return.
  //! @param column Column index of the element to return.
  [[nodiscard]] constexpr auto &&operator[](this auto &&self, std::size_t row,
                                            std::size_t column)
    requires is_uniform_typed_matrix<typed_matrix>;

  //! @brief Access the specified element.
  //!
  //! @details Applicable to one-dimension matrix: column- or row-vector.
  //! Applicable to single-type matrix: uniform type of all elements.
  //! Returns a reference to the element at the specified location.
  //!
  //! @param self Explicit object parameter deducing this: not user specified.
  //! @param index Position of the element to return.
  [[nodiscard]] constexpr auto &&operator()(this auto &&self, std::size_t index)
    requires is_uniform_typed_matrix<typed_matrix> and
             is_one_dimension_typed_matrix<typed_matrix>;

  //! @brief Access the specified element.
  //!
  //! @details Applicable to single-type matrix: uniform type of all elements.
  //! Returns a reference to the element at the specified location.
  //!
  //! @param self Explicit object parameter deducing this: not user specified.
  //! @param row Row index of the element to return.
  //! @param column Column index of the element to return.
  [[nodiscard]] constexpr auto &&operator()(this auto &&self, std::size_t row,
                                            std::size_t column)
    requires is_uniform_typed_matrix<typed_matrix>;

  //! @brief Access the specified element with compile-time bound checking.
  //!
  //! @details Returns a strongly typed reference to the element at the
  //! specified location.
  //!
  //! @tparam Row Row index of the element to return.
  //! @tparam Column Column index of the element to return.
  template <std::size_t Row, std::size_t Column>
  [[nodiscard]] constexpr auto
  at() -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, Row,
                       Column> &
    requires tla::in_range<
                 Row, 0,
                 typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows> &&
             tla::in_range<
                 Column, 0,
                 typed_matrix<Matrix, RowIndexes, ColumnIndexes>::columns>;

  //! @brief Access the specified element with compile-time bound checking.
  //!
  //! @details Returns a strongly typed element at the specified location.
  //!
  //! @tparam Row Row index of the element to return.
  //! @tparam Column Column index of the element to return.
  template <std::size_t Row, std::size_t Column>
  [[nodiscard]] constexpr auto
  at() const -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>,
                             Row, Column>
    requires tla::in_range<
                 Row, 0,
                 typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows> &&
             tla::in_range<
                 Column, 0,
                 typed_matrix<Matrix, RowIndexes, ColumnIndexes>::columns>;

  //! @brief Access the specified element with compile-time bound checking.
  //!
  //! @details Returns a strongly typed element at the specified location.
  //! Applicable to one-dimension matrix: column-vector.
  //!
  //! @tparam Index Position of the element to return.
  template <std::size_t Index>
  [[nodiscard]] constexpr auto at()
      -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, Index, 0>
          &
    requires is_column_typed_matrix<
                 typed_matrix<Matrix, RowIndexes, ColumnIndexes>> &&
             tla::in_range<
                 Index, 0,
                 typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows>;

  //! @brief Access the specified element with compile-time bound checking.
  //!
  //! @details Returns a strongly typed element at the specified location.
  //! Applicable to one-dimension matrix: column-vector.
  //!
  //! @tparam Index Position of the element to return.
  template <std::size_t Index>
  [[nodiscard]] constexpr auto at() const
      -> tla::element<typed_matrix<Matrix, RowIndexes, ColumnIndexes>, Index, 0>
    requires is_column_typed_matrix<
                 typed_matrix<Matrix, RowIndexes, ColumnIndexes>> &&
             tla::in_range<
                 Index, 0,
                 typed_matrix<Matrix, RowIndexes, ColumnIndexes>::rows>;

  //! @brief Direct access to the underlying storage.
  //!
  //! @details Reference to the underlying element storage.
  [[nodiscard]] constexpr auto &&data(this auto &&self);

  //! @}

private:
  //! @name Private Member Variables
  //! @{

  //! @brief Underlying algebraic backend data storage.
  Matrix storage;

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

//! @brief Typed matrix element conversions customization point.
//!
//! @details Specialize this template to allow conversion of element's type and
//! underlying type.
//!
//! @todo The call operator should be static once MSVC lands the support.
template <typename To, typename From> struct element_caster {
  [[nodiscard]] constexpr To operator()(From value) const;
};

//! @}

//! @brief Matrix element conversion customization point.
//!
//! @details Specialization of the element caster function objects allows the
//! end-user to permit underlying type conversions.
template <typename To, typename From>
static inline constexpr element_caster<To, From> cast{};
} // namespace fcarouge

#include "typed_linear_algebra_internal/cast.tpp"
#include "typed_linear_algebra_internal/operation.tpp"
#include "typed_linear_algebra_internal/typed_linear_algebra.tpp"

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_HPP
