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
concept same_as_typed_matrix = tla::same_as_typed_matrix<Type>;

//! @brief Concept of a typed matrix in which all element types are the same.
//!
//! @details Matrices with uniform types are type safe even with the traditional
//! access operators.
//!
//! @note A matrix may be uniform with different row and column indexes.
//!
//! @todo Identical types may be too conservative. Convertible may be
//! intended.
template <typename Type>
concept uniform_typed_matrix = tla::uniform_typed_matrix<Type>;

//! @brief Concept of a typed matrix with only one dimension, row, or column.
template <typename Type>
concept one_dimension_typed_matrix = tla::one_dimension_typed_matrix<Type>;

//! @brief Concept of a row, one-dimension typed matrix, vector.
template <typename Type>
concept row_typed_matrix = tla::row_typed_matrix<Type>;

//! @brief Concept of a column, one-dimension typed matrix, vector.
template <typename Type>
concept column_typed_matrix = tla::column_typed_matrix<Type>;

//! @brief Concept of a singleton, one-element, one-dimension typed matrix type.
template <typename Type>
concept singleton_typed_matrix = tla::singleton_typed_matrix<Type>;

//! @brief Concept of matrices of the same shape.
//!
//! @details The same shape of the matrices, that is they have the same number
//! of rows and columns.
template <typename Lhs, typename Rhs>
concept same_shape = tla::same_shape<Lhs, Rhs>;

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

  //! @brief The type of the composed matrix.
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

  //! @brief The number of dimensions in the matrix.
  static inline constexpr std::size_t rank{tla::rank<rows, columns>};

  //! @}

  //! @name Public Member Functions
  //! @{

  //! @brief Destruct a default typed matrix.
  constexpr ~typed_matrix() = default;

  //! @brief Construct a zero-initialized typed matrix.
  constexpr typed_matrix()
    requires std::default_initializable<Matrix>;

  //! @brief Copy construct a typed matrix.
  constexpr typed_matrix(const typed_matrix &other) = default;

  //! @brief Copy assign a typed matrix.
  constexpr typed_matrix &operator=(const typed_matrix &other) = default;

  //! @brief Move construct a typed matrix.
  constexpr typed_matrix(typed_matrix &&other) = default;

  //! @brief Move construct a typed matrix.
  constexpr typed_matrix &operator=(typed_matrix &&other) = default;

  //! @brief Copy construct generalization of a compatible typed matrix.
  //!
  //! @details Implicit conversions expected per default equivalency.
  constexpr explicit(false)
      typed_matrix(const same_as_typed_matrix auto &other);

  //! @brief Copy assign generalization of a compatible typed matrix.
  constexpr typed_matrix &operator=(const same_as_typed_matrix auto &other);

  //! @brief Move construct generalization of a compatible typed matrix.
  //!
  //! @details Implicit conversions expected per default equivalency.
  constexpr explicit(false) typed_matrix(same_as_typed_matrix auto &&other);

  //! @brief Move assign generalization of a compatible typed matrix.
  constexpr typed_matrix &operator=(same_as_typed_matrix auto &&other);

  //! @brief Convert construct a singleton typed matrix from a single value.
  //!
  //! @details Applicable to singleton matrix: one element.
  //!
  //! @param value Element of compatible type.
  constexpr explicit typed_matrix(const auto &value)
    requires singleton_typed_matrix<typed_matrix>;

  //! @brief Convert copy assign a singleton typed matrix from a single value.
  //!
  //! @details Applicable to singleton matrix: one element.
  //!
  //! @param value Element of compatible type.
  constexpr typed_matrix &operator=(const auto &value)
    requires singleton_typed_matrix<typed_matrix>;

  //! @brief Convert construct one-dimension uniformly typed matrix from array.
  //!
  //! @details Applicable to one-dimension matrix: column- or row-vector.
  //! Applicable to single-type matrix: uniform type of all elements.
  //! Single-argument constructors taking arrays of data should get implicit
  //! constructors.
  constexpr explicit(false)
      typed_matrix(const element<0, 0> (&elements)[rows * columns])
    requires uniform_typed_matrix<typed_matrix> and
             one_dimension_typed_matrix<typed_matrix>;

  //! @brief Copy assign one-dimension uniformly typed matrix from array.
  //!
  //! @details Applicable to one-dimension matrix: column- or row-vector.
  //! Applicable to single-type matrix: uniform type of all elements.
  //! Single-argument constructors taking arrays of data should get implicit
  //! constructors.
  constexpr typed_matrix &
  operator=(const element<0, 0> (&elements)[rows * columns])
    requires uniform_typed_matrix<typed_matrix> and
             one_dimension_typed_matrix<typed_matrix>;

  //! @brief Convert construct a uniformly typed matrix from list-initializers.
  //!
  //! @details Applicable to matrix of uniform elements type. Single-argument
  //! constructors taking arrays of data should get implicit constructors.
  //!
  //! @param row_list List-initializers of list-initializer of elements.
  template <typename Type>
  constexpr typed_matrix(
      std::initializer_list<std::initializer_list<Type>> row_list)
    requires uniform_typed_matrix<typed_matrix>;

  //! @brief Convert construct a row or column typed vector from elements.
  //!
  //! @details Applicable to one-dimension matrix. The first and second value
  //! parameter help in determining which constructor should the compiler call.
  //!
  //! @param first_value First element.
  //! @param second_value Second element.
  //! @param values Other elements.
  constexpr typed_matrix(const auto &first_value, const auto &second_value,
                         const auto &...values)
    requires one_dimension_typed_matrix<typed_matrix>;

  //! @brief Convert construct a typed matrix from an underlying matrix.
  //!
  //! @warning Useful for operations implementation where underlying data
  //! constrution is needed. Not recommended for convenience construction due to
  //! absence of type validation.
  //!
  //! @note Alternative design could evaluate feasability of private
  //! constructor, operator friendship, attorney-client, or key idioms.
  constexpr explicit typed_matrix(const Matrix &other);

  //! @brief Access the singleton typed matrix element.
  //!
  //! @details Applicable to singleton matrix: one element. Returns the unique
  //! element of the typed matrix.
  [[nodiscard]] constexpr explicit operator element<0, 0>(this auto &&self)
    requires singleton_typed_matrix<typed_matrix>;

  //! @brief Access the specified element.
  //!
  //! @details Applicable to matrices with identical element types. No bound
  //! checking. Returns a strongly typed element at the specified location. A
  //! reference is returned for non-const calls.
  //!
  //! @tparam Indexes Type(s) of the indexes. Use a template pack because some
  //! compilers have internal compiler errors with a placeholder type specifier.
  //! @param self Explicit object parameter deducing this: not user specified.
  //! @param indexes Position(s) of the element to return.
  template <typename... Indexes>
  [[nodiscard]] constexpr decltype(auto) operator[](this auto &&self,
                                                    Indexes... indexes)
    requires uniform_typed_matrix<typed_matrix> and (sizeof...(Indexes) >= rank)
  ;

  //! @brief Access the specified element.
  //!
  //! @details Applicable to matrices with identical element types. No bound
  //! checking. Returns a strongly typed element at the specified location. A
  //! reference is returned for non-const calls.
  //!
  //! @tparam Indexes Type(s) of the indexes. Use a template pack because some
  //! compilers have internal compiler errors with a placeholder type specifier.
  //! @param self Explicit object parameter deducing this: not user specified.
  //! @param indexes Position(s) of the element to return.
  template <typename... Indexes>
  [[nodiscard]] constexpr decltype(auto) operator()(this auto &&self,
                                                    Indexes... indexes)
    requires uniform_typed_matrix<typed_matrix> and (sizeof...(Indexes) >= rank)
  ;

  //! @brief Access the specified element with compile-time bound checking.
  //!
  //! @details Returns a strongly typed element at the specified location. A
  //! reference is returned for non-const calls.
  //!
  //! @tparam Indexes Position(s) of the element to return. Both row and column
  //! indexes for any matrix, one index for one-dimension matrices, no index for
  //! singleton matrices.
  template <auto... Indexes>
  [[nodiscard]] constexpr decltype(auto) at(this auto &&self)
    requires(sizeof...(Indexes) >= rank);

  //! @brief Direct access to the underlying storage.
  //!
  //! @details Reference to the underlying element storage.
  //!
  //! @warning Useful for operations implementation where underlying data
  //! access is needed. Not recommended for convenience access due to
  //! absence of type validation.
  [[nodiscard]] constexpr decltype(auto) data(this auto &&self);

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
template <typename To, typename From> struct element_caster {
  [[nodiscard]] static constexpr To operator()(From value);
};

//! @}

//! @name Algorithms
//! @{

[[nodiscard]] constexpr bool operator==(const same_as_typed_matrix auto &lhs,
                                        const same_as_typed_matrix auto &rhs);
[[nodiscard]] constexpr bool operator==(const singleton_typed_matrix auto &lhs,
                                        const auto &rhs);

[[nodiscard]] constexpr auto operator+(const same_as_typed_matrix auto &lhs,
                                       const same_as_typed_matrix auto &rhs);

[[nodiscard]] constexpr auto operator-(const same_as_typed_matrix auto &lhs,
                                       const same_as_typed_matrix auto &rhs);

[[nodiscard]] constexpr auto operator*(const same_as_typed_matrix auto &lhs,
                                       const same_as_typed_matrix auto &rhs);

[[nodiscard]] constexpr auto operator/(const same_as_typed_matrix auto &lhs,
                                       const same_as_typed_matrix auto &rhs);

[[nodiscard]] constexpr auto transposed(const same_as_typed_matrix auto &value);

//! @}

//! @name Adaptors
//! @{

//! @brief Matrix element conversion customization point.
//!
//! @details Specialization of the element caster function objects allows the
//! end-user to permit underlying type conversions.
template <typename To, typename From>
static inline constexpr element_caster<To, From> cast{};

//! @brief Factory function for partial template deduction.
//!
//! @warning Useful for operations implementation where underlying data
//! constrution is needed. Not recommended for convenience construction due to
//! absence of type validation.
//!
//! @note Alternative design could evaluate feasability of hiding this support.
template <typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto make_typed_matrix(auto &&value);

//! @brief Get function argument-dependent lookup overload.
//!
//! @details Also provides support for structured bindings.
template <int Index>
decltype(auto) get(one_dimension_typed_matrix auto &&value);

//! @}

} // namespace fcarouge

#include "typed_linear_algebra_internal/algorithm/add.tpp"
#include "typed_linear_algebra_internal/algorithm/divide.tpp"
#include "typed_linear_algebra_internal/algorithm/equal_to.tpp"
#include "typed_linear_algebra_internal/algorithm/matrix_product.tpp"
#include "typed_linear_algebra_internal/algorithm/product.tpp"
#include "typed_linear_algebra_internal/algorithm/scale.tpp"
#include "typed_linear_algebra_internal/algorithm/substract.tpp"
#include "typed_linear_algebra_internal/algorithm/transposed.tpp"
#include "typed_linear_algebra_internal/cast.tpp"
#include "typed_linear_algebra_internal/tuple.tpp"
#include "typed_linear_algebra_internal/typed_linear_algebra.tpp"

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_HPP
