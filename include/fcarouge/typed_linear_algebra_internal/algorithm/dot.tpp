/* Typed Linear Algebra
Version 0.4.0
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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_DOT_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_DOT_TPP

#if __has_include(<linalg>)

#include <linalg>

#endif

namespace fcarouge {
[[nodiscard]] constexpr auto dot(const rank_typed_matrix<1> auto &lhs,
                                 const rank_typed_matrix<1> auto &rhs);

namespace internal {
[[nodiscard]] constexpr auto dot(const rank_typed_matrix<1> auto &lhs,
                                 const rank_typed_matrix<1> auto &rhs) {
#if __has_include(<linalg>)
  using std::linalg::dot;
#endif

  if constexpr (requires { lhs.data().dot(rhs.data()); }) {
    return lhs.data().dot(rhs.data());
  } else if constexpr (requires { fcarouge::dot(lhs.data(), rhs.data()); }) {
    return fcarouge::dot(lhs.data(), rhs.data());
  } else if constexpr (requires { dot(lhs.data(), rhs.data()); }) {
    return dot(lhs.data(), rhs.data());
  } else if constexpr (requires {
                         dot(tla::as_vector_span(lhs),
                             tla::as_vector_span(rhs));
                       }) {
    return dot(tla::as_vector_span(lhs), tla::as_vector_span(rhs));
  } else {
    static_assert(sizeof(lhs) == 0,
                  "Dot is not supported for this linear algebra backend.");
  }
}
} // namespace internal

//! @brief Sum of the products of the corresponding elements of two vectors.
//!
//! @details Row or column orientation is not required to match between the
//! two vectors, only their element count. Delegated to the linear algebra
//! backend via `internal::dot`.
//!
//! @param lhs The first vector.
//! @param rhs The second vector.
//!
//! @return The sum of the products of the corresponding elements of `lhs`
//! and `rhs`.
[[nodiscard]] constexpr auto dot(const rank_typed_matrix<1> auto &lhs,
                                 const rank_typed_matrix<1> auto &rhs) {
  using lhs_matrix = std::remove_cvref_t<decltype(lhs)>;
  using rhs_matrix = std::remove_cvref_t<decltype(rhs)>;

  static_assert(lhs_matrix::rows * lhs_matrix::columns ==
                    rhs_matrix::rows * rhs_matrix::columns,
                "Dot product requires vectors of the same size.");

  using first_term = tla::product<typename lhs_matrix::template element<0>,
                                  typename rhs_matrix::template element<0>>;

  tla::for_constexpr<lhs_matrix::rows * lhs_matrix::columns>([](auto i) {
    using term = tla::product<typename lhs_matrix::template element<i>,
                              typename rhs_matrix::template element<i>>;

    static_assert(std::is_convertible_v<term, first_term>,
                  "Dot product requires compatible element types.");
  });

  using underlying = typename lhs_matrix::underlying;

  return cast<first_term, underlying>(internal::dot(lhs, rhs));
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_DOT_TPP
