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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_UTILITY_HPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_UTILITY_HPP

#include <concepts>
#include <cstddef>
#include <tuple>
#include <type_traits>

namespace fcarouge::typed_linear_algebra_internal {

//! @brief Arithmetic concept.
//!
//! @details Any integer or floating point type.
template <typename Type>
concept arithmetic = std::integral<Type> || std::floating_point<Type>;

template <typename Type> struct repacker {
  using type = Type;
};

template <template <typename...> typename Pack, typename... Types>
struct repacker<Pack<Types...>> {
  using type = std::tuple<Types...>;

  static inline constexpr std::size_t size{sizeof...(Types)};
};

template <auto Value, auto... Values> struct first_value {
  static constexpr auto value{Value};
};

template <typename Type> struct not_implemented {
  template <auto Size>
  inline constexpr explicit not_implemented(
      [[maybe_unused]] const char (&message)[Size]) {
    // The argument message is printed in the compiler error output.
  }

  static constexpr auto type_dependent_false{sizeof(Type) != sizeof(Type)};
  static constexpr auto missing{type_dependent_false};

  static_assert(missing, "This type is not implemented. See compiler message.");
};

template <typename Pack> using repack = repacker<Pack>::type;

//! @brief Unpack the first type of the type template parameter pack.
//!
//! @details Shorthand for `std::tuple_element_t<0, std::tuple<Types...>>`.
template <typename... Types>
using first = std::tuple_element_t<0, std::tuple<Types...>>;

//! @brief Compile-time for loop.
//!
//! @details Help compilers with non-type template parameters on members.
template <std::size_t Begin, std::size_t End, std::size_t Increment,
          typename Function>
inline constexpr void for_constexpr(Function &&function) {
  if constexpr (Begin < End) {
    function(std::integral_constant<std::size_t, Begin>());
    typed_linear_algebra_internal::for_constexpr<Begin + Increment, End,
                                                 Increment>(
        std::forward<Function>(function));
  }
}

//! @brief Size of tuple-like types.
//!
//! @details Convenient short form. In place of `std::tuple_size_v`.
template <typename Pack> inline constexpr auto size{repacker<Pack>::size};

template <auto... Values>
inline constexpr auto first_v{first_value<Values...>::value};

template <typename Type, std::size_t Size> struct tupler {
  template <typename = std::make_index_sequence<Size>> struct helper;

  template <std::size_t... Indexes>
  struct helper<std::index_sequence<Indexes...>> {
    template <std::size_t> using wrap = Type;

    using type = std::tuple<wrap<Indexes>...>;
  };

  using type = typename helper<>::type;
};

//! @brief An alias for making a tuple of the same type.
template <typename Type, std::size_t Size>
using tuple_n_type = typename tupler<Type, Size>::type;

//! @name Types
//! @{

//! @brief Linear algebra divides expression type specialization point.
//!
//! @details Matrix division is a mathematical abuse of terminology. Informally
//! defined as multiplication by the inverse. Similarly to division by zero in
//! real numbers, there exists matrices that are not invertible. Remember the
//! division operation is not commutative. Matrix inversion can be avoided by
//! solving `X * rhs = lhs` for `rhs` through a decomposer. There exists several
//! ways to decompose and solve the equation. Implementations trade off
//! numerical stability, triangularity, symmetry, space, time, etc. Dividing an
//! `R1 x C` matrix by an `R2 x C` matrix results in an `R1 x R2` matrix.
template <typename Lhs, typename Rhs> struct divides {
  [[nodiscard]] inline constexpr auto
  operator()(const Lhs &lhs, const Rhs &rhs) const -> decltype(lhs / rhs);
};

//! @brief Divider helper type.
template <typename Lhs, typename Rhs>
using quotient =
    std::invoke_result_t<divides<Lhs, Rhs>, const Lhs &, const Rhs &>;

//! @brief Linear algebra evaluates override expression lazy evaluation
//! specialization point.
template <typename Type> struct evaluates {
  [[nodiscard]] inline constexpr auto operator()() const -> Type;
};

//! @brief Evaluater helper type.
template <typename Type> using evaluate = std::invoke_result_t<evaluates<Type>>;

//! @brief Linear algebra transposes specialization point.
template <typename Type> struct transposes {
  [[nodiscard]] inline constexpr auto operator()(const Type &value) const {
    return value;
  }
};

template <typename Type>
  requires requires(Type value) { value.transpose(); }
struct transposes<Type> {
  [[nodiscard]] inline constexpr auto operator()(const Type &value) const {
    return value.transpose();
  }
};

//! @brief Transposer helper type.
template <typename Type>
using transpose = std::invoke_result_t<transposes<Type>, const Type &>;

//! @brief Transpose helper function.
//!
//! @details Enable readable linear algebra notation.
template <typename Type> auto t(const Type &value) {
  return transposes<Type>{}(value);
}

//! @brief Type multiplies expression type specialization point.
template <typename Lhs, typename Rhs> struct multiplies {
  [[nodiscard]] inline constexpr auto
  operator()(const Lhs &lhs, const Rhs &rhs) const -> decltype(lhs * rhs);
};

//! @brief Helper type to deduce the result type of the product.
template <typename Lhs, typename Rhs>
using product =
    std::invoke_result_t<multiplies<Lhs, Rhs>, const Lhs &, const Rhs &>;

//! @brief Type minus, subtraction expression type specialization point.
template <typename Lhs, typename Rhs> struct minus {
  [[nodiscard]] inline constexpr auto
  operator()(const Lhs &lhs, const Rhs &rhs) const -> decltype(lhs - rhs);
};

//! @}

//! @brief The underlying storage type of the matrix's elements.
template <typename Matrix>
using underlying_t =
    std::remove_cvref_t<decltype(std::declval<Matrix>()(0, 0))>;

//! @brief The type of the element at the given matrix indexes position.
template <typename Matrix, std::size_t RowIndex, std::size_t ColumnIndex>
using element =
    product<std::tuple_element_t<RowIndex, typename Matrix::row_indexes>,
            std::tuple_element_t<ColumnIndex, typename Matrix::column_indexes>>;

//! @brief Every element types of the matrix are the same.
//!
//! @details Matrices with uniform types are type safe even with the traditional
//! operators.
//!
//! @note A matrix may be uniform with different row and column indexes.
//!
//! @todo There may be a way to write this concepts via two fold expressions.
template <typename Matrix>
concept uniform = []() {
  bool result{true};

  typed_linear_algebra_internal::for_constexpr<0, Matrix::rows, 1>([&result](auto i) {
    typed_linear_algebra_internal::for_constexpr<0, Matrix::columns, 1>([&result, &i](auto j) {
      result &= std::is_same_v<element<Matrix, i, j>, element<Matrix, 0, 0>>;
    });
  });

  return result;
}();

//! @brief The index is within the range, inclusive.
template <std::size_t Index, std::size_t Begin, std::size_t End>
concept in_range = Begin <= Index && Index <= End;

//! @brief The given matrix is a single column.
template <typename Matrix>
concept column = Matrix::columns == 1;

//! @brief The matrix is a single row.
template <typename Matrix>
concept row = Matrix::rows == 1;

//! @brief The given matrix is a single dimension, that is a row or a column.
template <typename Matrix>
concept one_dimension = column<Matrix> || row<Matrix>;

//! @brief The given row and column indexes form a singleton matrix.
template <typename Matrix>
concept singleton = column<Matrix> && row<Matrix>;

//! @brief The packs have the same count of types.
template <typename Pack1, typename Pack2>
concept same_size = size<Pack1> == size<Pack2>;

//! @brief Element traits for conversions.
template <typename Underlying, typename Type> struct element_traits {
  [[nodiscard]] static inline constexpr Underlying to_underlying(Type value) {
    return value;
  }

  [[nodiscard]] static inline constexpr Type &
  from_underlying(Underlying &value) {
    return value;
  }
};

} // namespace fcarouge::typed_linear_algebra_internal

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_UTILITY_HPP
