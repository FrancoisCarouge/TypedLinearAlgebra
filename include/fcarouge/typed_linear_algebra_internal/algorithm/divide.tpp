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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_DIVIDE_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_DIVIDE_TPP

namespace fcarouge {
namespace tla = typed_linear_algebra_internal;

//! @details Matrix division is a mathematical abuse of terminology. Informally
//! defined as multiplication by the inverse. Similarly to division by zero in
//! real numbers, there exists matrices that are not invertible. Remember the
//! division operation is not commutative. Matrix inversion can be avoided by
//! solving `X * rhs = lhs` for `rhs` through a decomposer. There exists several
//! ways to decompose and solve the equation. Implementations trade off
//! numerical stability, triangularity, symmetry, space, time, etc. Dividing an
//! `R1 x C` matrix by an `R2 x C` matrix results in an `R1 x R2` matrix.
//!
//! @todo Combine? Generalize?
[[nodiscard]] constexpr auto operator/(const same_as_typed_matrix auto &lhs,
                                       const same_as_typed_matrix auto &rhs) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;

  //! @todo Convert to a requires clause?
  static_assert(lhs_matrix::columns == rhs_matrix::columns,
                "Matrix division requires compatible sizes.");

  using lhs_row_indexes = typename lhs_matrix::row_indexes;
  using lhs_column_indexes = typename lhs_matrix::column_indexes;
  using rhs_row_indexes = typename rhs_matrix::row_indexes;
  using rhs_column_indexes = typename rhs_matrix::column_indexes;
  using row_indexes =
      tla::quotient<lhs_row_indexes,
                    std::tuple_element_t<0, lhs_column_indexes>>;
  using column_indexes =
      tla::quotient<std::tuple_element_t<0, rhs_column_indexes>,
                    rhs_row_indexes>;

  //! @todo Add type verification, perhaps with a generalization of the
  //! multiplication verification?

  return make_typed_matrix<row_indexes, column_indexes>(lhs.data() /
                                                        rhs.data());
}

[[nodiscard]] constexpr auto operator/(const singleton_typed_matrix auto &lhs,
                                       const auto &rhs)
  requires(not same_as_typed_matrix<decltype(rhs)>)
{
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using element = typename matrix::template element<0, 0>;

  return element{lhs} / rhs;
}

[[nodiscard]] constexpr auto operator/(const auto &lhs,
                                       const singleton_typed_matrix auto &rhs)
  requires(not same_as_typed_matrix<decltype(lhs)>)
{
  //! @todo Should there be constraints on the type?
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using element = typename matrix::template element<0, 0>;

  return lhs / element{rhs};
}

[[nodiscard]] constexpr auto operator/(const singleton_typed_matrix auto &lhs,
                                       const singleton_typed_matrix auto &rhs) {
  //! @todo Should there be constraints on the type?
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;
  using lhs_element = typename lhs_matrix::template element<0, 0>;
  using rhs_element = typename rhs_matrix::template element<0, 0>;

  return lhs_element{lhs} / rhs_element{rhs};
}

[[nodiscard]] constexpr auto operator/(const auto &lhs,
                                       const column_typed_matrix auto &rhs)
  requires(not same_as_typed_matrix<decltype(lhs)>)
{
  //! @todo Should there be constraints on the type?
  using type = std::remove_cvref_t<decltype(lhs)>;
  using matrix = std::remove_cvref_t<decltype(rhs)>;
  using underlying = typename matrix::underlying;
  using row_indexes = tla::quotient<type, typename matrix::column_indexes>;
  using column_indexes =
      tla::quotient<std::identity, typename matrix::row_indexes>;

  //! @todo Add type verification, perhaps with a generalization of the
  //! multiplication verification?

  return make_typed_matrix<row_indexes, column_indexes>(
      cast<underlying, type>(lhs) / rhs.data());
}

[[nodiscard]] constexpr auto operator/(const same_as_typed_matrix auto &lhs,
                                       const auto &rhs)
  requires(not same_as_typed_matrix<decltype(rhs)>)
{
  //! @todo Should there be constraints on the type?
  using type = std::remove_cvref_t<decltype(rhs)>;
  using matrix = std::remove_cvref_t<decltype(lhs)>;
  using underlying = typename matrix::underlying;
  using row_indexes = tla::quotient<typename matrix::row_indexes, type>;
  using column_indexes = typename matrix::column_indexes;

  return make_typed_matrix<row_indexes, column_indexes>(
      lhs.data() / cast<underlying, type>(rhs));
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_DIVIDE_TPP
