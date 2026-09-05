/* Typed Linear Algebra
Version 0.3.0
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

#include "fcarouge/typed_linear_algebra_forward.hpp"

#include <concepts>
#include <cstddef>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

namespace fcarouge::typed_linear_algebra_internal {
//! @brief Linear algebra divides expression type specialization point.
//!
//! @details Matrix division is a mathematical abuse of terminology. Informally
//! defined as multiplication by the inverse. Similarly to division by zero in
//! real numbers, there exist matrices that are not invertible. Remember the
//! division operation is not commutative. Matrix inversion can be avoided by
//! solving `X * rhs = lhs` for `rhs` through a decomposer. There exist several
//! ways to decompose and solve the equation. Implementations trade off
//! numerical stability, triangularity, symmetry, space, time, etc. Dividing an
//! `R1 x C` matrix by an `R2 x C` matrix results in an `R1 x R2` matrix.
template <typename Lhs, typename Rhs> struct divides {
  [[nodiscard]] static constexpr auto operator()(const Lhs &lhs, const Rhs &rhs)
      -> decltype(lhs / rhs);
};

template <typename Lhs, typename Rhs>
using quotient =
    std::invoke_result_t<divides<Lhs, Rhs>, const Lhs &, const Rhs &>;

template <> struct divides<std::identity, std::identity> {
  [[nodiscard]] static constexpr auto operator()(const std::identity &lhs,
                                                 const std::identity &rhs)
      -> std::identity;
};

template <typename Lhs> struct divides<Lhs, std::identity> {
  [[nodiscard]] static constexpr auto operator()(const Lhs &lhs,
                                                 const std::identity &rhs)
      -> Lhs;
};

template <typename Rhs> struct divides<std::identity, Rhs> {
  [[nodiscard]] static constexpr auto operator()(const std::identity &lhs,
                                                 const Rhs &rhs)
      -> quotient<quotient<Rhs, Rhs>, Rhs>;
};

template <typename Lhs, typename... Types>
struct divides<Lhs, std::tuple<Types...>> {
  [[nodiscard]] static constexpr auto
  operator()(const Lhs &lhs, const std::tuple<Types...> &rhs)
      -> std::tuple<quotient<Lhs, Types>...>;
};

template <typename... Types>
struct divides<std::identity, std::tuple<Types...>> {
  [[nodiscard]] static constexpr auto
  operator()(const std::identity &lhs, const std::tuple<Types...> &rhs)
      -> std::tuple<quotient<std::identity, Types>...>;
};

template <typename Rhs, typename... Types>
struct divides<std::tuple<Types...>, Rhs> {
  [[nodiscard]] static constexpr auto
  operator()(const std::tuple<Types...> &lhs, const Rhs &rhs)
      -> std::tuple<quotient<Types, Rhs>...>;
};

template <typename... Types>
struct divides<std::tuple<Types...>, std::identity> {
  [[nodiscard]] static constexpr auto
  operator()(const std::tuple<Types...> &lhs, const std::identity &rhs)
      -> std::tuple<Types...>;
};

template <typename... Types1, typename... Types2>
struct divides<std::tuple<Types1...>, std::tuple<Types2...>> {
  [[nodiscard]] static constexpr auto
  operator()(const std::tuple<Types1...> &lhs, const std::tuple<Types2...> &rhs)
      -> std::tuple<quotient<Types1, Types2>...>;
};

template <typename Lhs, typename Rhs> struct multiplies {
  [[nodiscard]] static constexpr auto operator()(const Lhs &lhs, const Rhs &rhs)
      -> decltype(lhs * rhs);
};

template <typename Lhs, typename Rhs>
using product = std::invoke_result_t<multiplies<Lhs, Rhs>, Lhs, Rhs>;

template <typename Lhs>
  requires(not std::same_as<Lhs, std::identity>)
struct multiplies<Lhs, std::identity> {
  [[nodiscard]] static constexpr auto operator()(const Lhs &lhs,
                                                 const std::identity &rhs)
      -> Lhs;
};

template <typename Rhs> struct multiplies<std::identity, Rhs> {
  [[nodiscard]] static constexpr auto operator()(const std::identity &lhs,
                                                 const Rhs &rhs) -> Rhs;
};

template <typename Lhs, typename... Types>
  requires(not std::same_as<Lhs, std::identity>)
struct multiplies<Lhs, std::tuple<Types...>> {
  [[nodiscard]] static constexpr auto
  operator()(const Lhs &lhs, const std::tuple<Types...> &rhs)
      -> std::tuple<product<Lhs, Types>...>;
};

template <typename Rhs, typename... Types>
struct multiplies<std::tuple<Types...>, Rhs> {
  [[nodiscard]] static constexpr auto
  operator()(const std::tuple<Types...> &lhs, const Rhs &rhs)
      -> std::tuple<product<Types, Rhs>...>;
};

//! @todo Some specialization with identity may be redundant, to remove.
template <typename... Types>
struct multiplies<std::tuple<Types...>, std::identity> {
  [[nodiscard]] static constexpr auto
  operator()(const std::tuple<Types...> &lhs, const std::identity &rhs)
      -> std::tuple<Types...>;
};

template <typename... Types1, typename... Types2>
struct multiplies<std::tuple<Types1...>, std::tuple<Types2...>> {
  [[nodiscard]] static constexpr auto
  operator()(const std::tuple<Types1...> &lhs, const std::tuple<Types2...> &rhs)
      -> std::tuple<product<Types1, Types2>...>;
};

template <std::size_t... Is, typename F>
constexpr void
for_constexpr_detail([[maybe_unused]] std::index_sequence<Is...> indexes,
                     F &&f) {
  (std::forward<F>(f)(std::integral_constant<std::size_t, Is>{}), ...);
}

//! @todo Remove for C++26 P1789 Library Support for Expansion Statements.
template <std::size_t Size, typename Function>
constexpr void for_constexpr(Function &&function) {
  for_constexpr_detail(std::make_index_sequence<Size>{},
                       std::forward<Function>(function));
}

template <typename Type> struct underlying {
  [[nodiscard]] static constexpr auto operator()()
    requires requires { typename Type::underlying; }
  {
    return typename Type::underlying{};
  }
  [[nodiscard]] static constexpr auto operator()()
    requires requires { typename Type::Scalar; }
  {
    return typename Type::Scalar{};
  }
  [[nodiscard]] static constexpr auto operator()()
    requires requires { typename Type::element_type; }
  {
    return std::remove_cvref_t<typename Type::element_type>{};
  }
};

template <typename Type>
using underlying_t = std::invoke_result_t<underlying<Type>>;

template <typename Type>
concept same_as_typed_matrix = std::same_as<
    std::remove_cvref_t<Type>,
    typed_matrix<typename std::remove_cvref_t<Type>::matrix,
                 typename std::remove_cvref_t<Type>::row_indexes,
                 typename std::remove_cvref_t<Type>::column_indexes>>;

template <std::size_t Rows, std::size_t Columns>
constexpr std::size_t rank{[] {
  if constexpr (Rows > 1 && Columns > 1) {
    return 2;
  } else if constexpr (Rows == 1 && Columns == 1) {
    return 0;
  } else {
    return 1;
  }
}()};

//! @brief The type of the element at the given row and column matrix
//! indexes.
//!
//! @details Internal helper, oblivious to the matrix's rank. Reused by the
//! rank-checked `element_t` specializations below and by algorithms that
//! must visit every row and column position of a matrix regardless of its
//! rank, for example homogeneity checks.
template <typename Type, std::size_t RowIndex, std::size_t ColumnIndex>
using element_at = std::remove_cvref_t<product<
    std::tuple_element_t<RowIndex,
                         typename std::remove_cvref_t<Type>::row_indexes>,
    std::tuple_element_t<ColumnIndex,
                         typename std::remove_cvref_t<Type>::column_indexes>>>;

//! @brief Linear algebra element type specialization point.
//!
//! @details The count of indexes must match the rank of the matrix: no
//! index for a singleton matrix, one index for a row or column vector, two
//! indexes for any other matrix. Naming `element_t<Type, Indexes...>::type`
//! is always well-formed regardless of the count of indexes matching the
//! rank of the matrix: standard library facilities may need to complete
//! this placeholder `type` while probing members excluded by their own
//! `requires` clause, so completion must stay harmless. Only an actual
//! construction attempt, the signature of a genuine mismatch, is
//! diagnosed, by deleting every constructor.
template <typename Type, std::size_t... Indexes> struct element_t {
  struct type {
    type() = delete;
    type(const type &) = delete;
    type(type &&) = delete;
    ~type() = delete;
    type &operator=(const type &) = delete;
    type &operator=(type &&) = delete;
  };
};

template <typename Type, std::size_t RowIndex, std::size_t ColumnIndex>
  requires(rank<std::remove_cvref_t<Type>::rows,
                std::remove_cvref_t<Type>::columns> == 2)
struct element_t<Type, RowIndex, ColumnIndex> {
  using type = element_at<Type, RowIndex, ColumnIndex>;
};

template <typename Type, std::size_t Index>
  requires(rank<std::remove_cvref_t<Type>::rows,
                std::remove_cvref_t<Type>::columns> == 1)
struct element_t<Type, Index> {
  using type = element_at<Type, Index / std::remove_cvref_t<Type>::columns,
                          Index % std::remove_cvref_t<Type>::columns>;
};

template <typename Type>
  requires(rank<std::remove_cvref_t<Type>::rows,
                std::remove_cvref_t<Type>::columns> == 0)
struct element_t<Type> {
  using type = element_at<Type, 0, 0>;
};

template <typename Type, std::size_t... Indexes>
using element = element_t<Type, Indexes...>::type;

template <typename Type> constexpr auto is_uniform_typed_matrix() -> bool {
  using matrix = std::remove_cvref_t<Type>;
  bool result{true};

  for_constexpr<matrix::rows>([&result](auto i) {
    for_constexpr<matrix::columns>([&result, &i](auto j) {
      result &=
          std::is_same_v<element_at<matrix, i, j>, element_at<matrix, 0, 0>>;
    });
  });

  return result;
}

template <typename Type>
concept uniform_typed_matrix =
    same_as_typed_matrix<Type> and is_uniform_typed_matrix<Type>();

template <typename Type>
concept column_typed_matrix =
    same_as_typed_matrix<Type> and (std::remove_cvref_t<Type>::columns == 1);

template <typename Type>
concept row_typed_matrix =
    same_as_typed_matrix<Type> and (std::remove_cvref_t<Type>::rows == 1);

template <typename Type, auto Rank>
concept rank_typed_matrix = same_as_typed_matrix<Type> and
                            (rank<std::remove_cvref_t<Type>::rows,
                                  std::remove_cvref_t<Type>::columns> == Rank);

template <typename Lhs, typename Rhs>
concept same_shape =
    same_as_typed_matrix<Lhs> and same_as_typed_matrix<Rhs> and
    (std::remove_cvref_t<Lhs>::rows == std::remove_cvref_t<Rhs>::rows) &&
    (std::remove_cvref_t<Lhs>::columns == std::remove_cvref_t<Rhs>::columns);

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

template <typename> struct is_integral_constant_t;

template <std::size_t Index>
struct is_integral_constant_t<std::integral_constant<std::size_t, Index>>
    : std::bool_constant<true> {};

template <typename>
struct is_integral_constant_t : std::bool_constant<false> {};

template <typename Type>
concept index = is_integral_constant_t<Type>::value;

//! @brief A type presenting the standard tuple-like protocol with exactly
//! `Size` gettable elements.
//!
//! @details Matches any type `std::tuple_size` reports as having `Size`
//! elements, each retrievable through an unqualified, index-templated `get`
//! call resolved by argument-dependent lookup, e.g. `mp_units::quantity`
//! (vector-represented quantities specialize `std::tuple_size` and
//! `std::tuple_element`, and expose `get<Index>`), `std::tuple`,
//! `std::array`, or `std::pair`. This is the same protocol structured
//! bindings rely on. Unlike duck-typing a single member name, this does not
//! collide with unrelated types: `std::chrono::duration`, for example, does
//! not specialize `std::tuple_size`.
template <typename Type, std::size_t Size>
concept tuple_like =
    requires {
      {
        std::tuple_size<std::remove_cvref_t<Type>>::value
      } -> std::convertible_to<std::size_t>;
    } and
    (std::tuple_size_v<std::remove_cvref_t<Type>> ==
     Size) and[]<std::size_t... Indexes>(std::index_sequence<Indexes...>) {
  return (requires(const Type &value) { get<Indexes>(value); } and ...);
}
(std::make_index_sequence<Size>());

//! @brief A tuple-like value convertible to the given one-dimension typed
//! matrix, but neither that typed matrix nor already convertible to its
//! underlying matrix.
//!
//! @details Practical to disambiguate the typed vector's tuple-like
//! convert-constructor from the copy/move and raw-matrix constructors: both
//! `Value` being another typed matrix and `Value` being convertible to
//! `TypedMatrix`'s underlying matrix are excluded, else those constructors
//! would compete with this one. See `tuple_like`.
template <typename Value, typename TypedMatrix>
concept other_tuple_like_vector =
    rank_typed_matrix<TypedMatrix, 1> and (not same_as_typed_matrix<Value>) and
    (not std::convertible_to<
        Value, typename std::remove_cvref_t<TypedMatrix>::matrix>) and
    tuple_like<Value, std::remove_cvref_t<TypedMatrix>::rows *
                          std::remove_cvref_t<TypedMatrix>::columns>;

template <char... Digits> constexpr auto parse_digits() -> std::size_t {
  static_assert((('0' <= Digits && Digits <= '9') && ...),
                "Characters must only be digits.");

  std::size_t number{0};
  ((number = number * 10 + (Digits - '0')), ...);

  return number;
}

//! @brief Concept of two types with an implicit conversion in at least one
//! direction.
//!
//! @details Satisfied when `Lhs` is implicitly convertible to `Rhs` or `Rhs`
//! is implicitly convertible to `Lhs`. The relation is symmetric and, because
//! every type is convertible to itself, reflexive. "Inter" here means "either
//! way", not "both ways".
//!
//! @tparam Lhs The first type of the pair.
//! @tparam Rhs The second type of the pair.
template <typename Lhs, typename Rhs>
concept are_interconvertible =
    std::is_convertible_v<Lhs, Rhs> or std::is_convertible_v<Rhs, Lhs>;

//! @brief Concept of two types with no implicit conversion in either direction.
//!
//! @details The negation of `are_interconvertible`. Never satisfied for two
//! identical types.
//!
//! @tparam Lhs The first type of the pair.
//! @tparam Rhs The second type of the pair.
template <typename Lhs, typename Rhs>
concept are_not_interconvertible = not are_interconvertible<Lhs, Rhs>;

//! @brief Whether every pair of distinct element positions of the typed matrix
//! holds types with no implicit conversion between them.
//!
//! @details Visits every ordered pair of row/column positions, skipping a
//! position paired with itself, and requires `are_not_interconvertible` for
//! each. Rank-oblivious: uses `element_at` so it applies to singleton, vector,
//! and two-dimension typed matrices alike.
//!
//! @tparam Type The typed matrix type to inspect.
//!
//! @todo Evaluate compile-time performance. Instantiation, evaluation, and
//! compilation of this function may be expensive. We may want to use a pure
//! concept approach. The ordered-pair traversal also tests each unordered
//! pair twice.
template <typename Type> constexpr auto is_distinct_typed_matrix() -> bool {
  using matrix = std::remove_cvref_t<Type>;
  bool result{true};

  for_constexpr<matrix::rows>([&result](auto i) {
    for_constexpr<matrix::columns>([&result, &i](auto j) {
      for_constexpr<matrix::rows>([&result, &i, &j](auto k) {
        for_constexpr<matrix::columns>([&result, &i, &j, &k](auto l) {
          if constexpr (i != k || j != l) {
            result &= are_not_interconvertible<element_at<matrix, i, j>,
                                               element_at<matrix, k, l>>;
          }
        });
      });
    });
  });

  return result;
}

//! @brief Concept of a typed matrix whose element types have no pairwise
//! implicit conversion.
//!
//! @details Satisfied when `Type` is a typed matrix and no element type is
//! implicitly convertible to another element type at a different position.
//! See the public `fcarouge::distinct_typed_matrix` for the rationale and the
//! relationship with `uniform_typed_matrix`.
template <typename Type>
concept distinct_typed_matrix =
    same_as_typed_matrix<Type> and is_distinct_typed_matrix<Type>();
} // namespace fcarouge::typed_linear_algebra_internal

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_UTILITY_HPP
