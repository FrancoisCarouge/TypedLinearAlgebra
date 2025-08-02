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
#include <functional>
#include <tuple>
#include <type_traits>

namespace fcarouge::typed_linear_algebra_internal {
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
  [[nodiscard]] constexpr auto
  operator()(const Lhs &lhs, const Rhs &rhs) const -> decltype(lhs / rhs);
};

template <typename Lhs, typename Rhs>
using quotient =
    std::invoke_result_t<divides<Lhs, Rhs>, const Lhs &, const Rhs &>;

template <> struct divides<std::identity, std::identity> {
  [[nodiscard]] constexpr auto
  operator()(const std::identity &lhs,
             const std::identity &rhs) const -> std::identity;
};

template <typename Lhs> struct divides<Lhs, std::identity> {
  [[nodiscard]] constexpr auto
  operator()(const Lhs &lhs, const std::identity &rhs) const -> Lhs;
};

template <typename Rhs> struct divides<std::identity, Rhs> {
  [[nodiscard]] constexpr auto
  operator()(const std::identity &lhs,
             const Rhs &rhs) const -> quotient<quotient<Rhs, Rhs>, Rhs>;
};

template <typename Lhs, typename... Types>
struct divides<Lhs, std::tuple<Types...>> {
  [[nodiscard]] constexpr auto operator()(const Lhs &lhs,
                                          const std::tuple<Types...> &rhs) const
      -> std::tuple<quotient<Lhs, Types>...>;
};

template <typename... Types>
struct divides<std::identity, std::tuple<Types...>> {
  [[nodiscard]] constexpr auto operator()(const std::identity &lhs,
                                          const std::tuple<Types...> &rhs) const
      -> std::tuple<quotient<std::identity, Types>...>;
};

template <typename Rhs, typename... Types>
struct divides<std::tuple<Types...>, Rhs> {
  [[nodiscard]] constexpr auto
  operator()(const std::tuple<Types...> &lhs,
             const Rhs &rhs) const -> std::tuple<quotient<Types, Rhs>...>;
};

template <typename... Types>
struct divides<std::tuple<Types...>, std::identity> {
  [[nodiscard]] constexpr auto
  operator()(const std::tuple<Types...> &lhs,
             const std::identity &rhs) const -> std::tuple<Types...>;
};

template <typename... Types1, typename... Types2>
struct divides<std::tuple<Types1...>, std::tuple<Types2...>> {
  [[nodiscard]] constexpr auto operator()(const std::tuple<Types1...> &lhs,
                                          const std::tuple<Types2...> &rhs)
      const -> std::tuple<quotient<Types1, Types2>...>;
};

template <typename Lhs, typename Rhs> struct multiplies {
  [[nodiscard]] constexpr auto
  operator()(const Lhs &lhs, const Rhs &rhs) const -> decltype(lhs * rhs);
};

template <typename Lhs, typename Rhs>
using product = std::invoke_result_t<multiplies<Lhs, Rhs>, Lhs, Rhs>;

template <typename Lhs> struct multiplies<Lhs, std::identity> {
  [[nodiscard]] constexpr auto
  operator()(const Lhs &lhs, const std::identity &rhs) const -> Lhs;
};

template <typename Rhs> struct multiplies<std::identity, Rhs> {
  [[nodiscard]] constexpr auto operator()(const std::identity &lhs,
                                          const Rhs &rhs) const -> Rhs;
};

template <typename Lhs, typename... Types>
struct multiplies<Lhs, std::tuple<Types...>> {
  [[nodiscard]] constexpr auto operator()(const Lhs &lhs,
                                          const std::tuple<Types...> &rhs) const
      -> std::tuple<product<Lhs, Types>...>;
};

template <typename Rhs, typename... Types>
struct multiplies<std::tuple<Types...>, Rhs> {
  [[nodiscard]] constexpr auto
  operator()(const std::tuple<Types...> &lhs,
             const Rhs &rhs) const -> std::tuple<product<Types, Rhs>...>;
};

template <typename... Types>
struct multiplies<std::tuple<Types...>, std::identity> {
  [[nodiscard]] constexpr auto
  operator()(const std::tuple<Types...> &lhs,
             const std::identity &rhs) const -> std::tuple<Types...>;
};

template <typename... Types1, typename... Types2>
struct multiplies<std::tuple<Types1...>, std::tuple<Types2...>> {
  [[nodiscard]] constexpr auto operator()(const std::tuple<Types1...> &lhs,
                                          const std::tuple<Types2...> &rhs)
      const -> std::tuple<product<Types1, Types2>...>;
};

template <std::size_t Begin, std::size_t End, std::size_t Increment,
          typename Function>
constexpr void for_constexpr(Function &&function) {
  if constexpr (Begin < End) {
    function(std::integral_constant<std::size_t, Begin>());
    for_constexpr<Begin + Increment, End, Increment>(
        std::forward<Function>(function));
  }
}

template <typename Type>
using underlying_t =
    std::remove_cvref_t<decltype(std::declval<std::remove_cvref_t<Type>>()(0,
                                                                           0))>;

template <typename Type>
concept is_typed_matrix = std::same_as<
    std::remove_cvref_t<Type>,
    typed_matrix<typename std::remove_cvref_t<Type>::matrix,
                 typename std::remove_cvref_t<Type>::row_indexes,
                 typename std::remove_cvref_t<Type>::column_indexes>>;

template <typename Type, std::size_t RowIndex, std::size_t ColumnIndex>
using element = std::remove_cvref_t<product<
    std::tuple_element_t<RowIndex,
                         typename std::remove_cvref_t<Type>::row_indexes>,
    std::tuple_element_t<ColumnIndex,
                         typename std::remove_cvref_t<Type>::column_indexes>>>;

//! @todo There may be a way to write this concepts via two fold expressions.
template <typename Type>
concept is_uniform_typed_matrix =
    is_typed_matrix<Type> and ([]() {
      bool result{true};

      for_constexpr<0, std::remove_cvref_t<Type>::rows, 1>([&result](auto i) {
        for_constexpr<0, std::remove_cvref_t<Type>::columns, 1>([&result,
                                                                 &i](auto j) {
          result &= std::is_same_v<element<Type, i, j>, element<Type, 0, 0>>;
        });
      });

      return result;
    }());

template <typename Type>
concept is_column_typed_matrix =
    is_typed_matrix<Type> and (std::remove_cvref_t<Type>::columns == 1);

template <typename Type>
concept is_row_typed_matrix =
    is_typed_matrix<Type> and (std::remove_cvref_t<Type>::rows == 1);

template <typename Type>
concept is_one_dimension_typed_matrix =
    is_typed_matrix<Type> and
    (is_column_typed_matrix<Type> or is_row_typed_matrix<Type>);

template <typename Type>
concept is_singleton_typed_matrix =
    is_typed_matrix<Type> and is_column_typed_matrix<Type> and
    is_row_typed_matrix<Type>;

template <typename Type, std::size_t Size> struct tupler {
  template <typename = std::make_index_sequence<Size>> struct helper;

  template <std::size_t... Indexes>
  struct helper<std::index_sequence<Indexes...>> {
    template <std::size_t> using wrap = Type;

    using type = std::tuple<wrap<Indexes>...>;
  };

  using type = typename helper<>::type;
};

template <typename Type, std::size_t Size>
using tuple_n_type = typename tupler<Type, Size>::type;

using identity_index = std::tuple<std::identity>;
} // namespace fcarouge::typed_linear_algebra_internal

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_UTILITY_HPP
