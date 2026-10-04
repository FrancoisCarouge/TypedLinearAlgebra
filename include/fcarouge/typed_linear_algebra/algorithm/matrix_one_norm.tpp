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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_ONE_NORM_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_ONE_NORM_TPP

#if __has_include(<linalg>)

#include <linalg>

#endif

#include <concepts>
#include <type_traits>

namespace fcarouge {
[[nodiscard]] constexpr auto
matrix_one_norm(const uniform_typed_matrix auto &value);

namespace typed_linear_algebra::internal {
//! @brief One norm of the backend storage, in the underlying type.
//!
//! @details Tries each backend's own spelling, most to least specific: an
//! Eigen `.cwiseAbs()` column sum maximum, a qualified recursion into
//! `nested_typed_eigen`'s composed typed matrix storage, an ADL-found
//! Armadillo `sum(abs(X), 0)` maximum, then `std::linalg::matrix_one_norm`.
//!
//! @note Neither Eigen's `lpNorm<1>()` nor Armadillo's `norm(X, 1)` is the
//! matrix one norm: the former is the sum of all absolute elements, and the
//! latter turns into it for any one-row or one-column input, for example a row
//! vector, whose matrix one norm is its largest absolute element instead.
[[nodiscard]] constexpr auto
matrix_one_norm(const uniform_typed_matrix auto &value) {
  if constexpr (requires {
                  value.data().cwiseAbs().colwise().sum().maxCoeff();
                }) {
    return value.data().cwiseAbs().colwise().sum().maxCoeff();
  } else if constexpr (requires { fcarouge::matrix_one_norm(value.data()); }) {
    return fcarouge::matrix_one_norm(value.data());
  } else if constexpr (requires { sum(abs(value.data()), 0).max(); }) {
    return sum(abs(value.data()), 0).max();
  }
#if __has_include(<linalg>)
  else if constexpr (requires { std::linalg::matrix_one_norm(value.data()); }) {
    return std::linalg::matrix_one_norm(value.data());
  }
#endif
  else {
    static_assert(sizeof(value) == 0,
                  "One norm is not supported for this linear algebra backend.");
  }
}
} // namespace typed_linear_algebra::internal

//! @brief One norm of a matrix.
//!
//! @details The maximum, over the columns, of the sum of the absolute values of
//! the column's elements. For a column vector, the sum of the absolute values
//! of its elements. For a row vector, the largest absolute value of its
//! elements. For a singleton, the absolute value. The elements must be of a
//! uniform type supporting addition, so that the norm shares the element's
//! type: a matrix of lengths has a length norm. Affine types, for example
//! `std::chrono::time_point`, are rejected.
//!
//! @param value The typed matrix.
//!
//! @return The one norm of `value`, of the element type.
//!
//! @see std::linalg::matrix_one_norm
[[nodiscard]] constexpr auto
matrix_one_norm(const uniform_typed_matrix auto &value) {
  using matrix = std::remove_cvref_t<decltype(value)>;
  using element = tla::element_at<matrix, 0, 0>;
  using underlying = typename matrix::underlying;

  static_assert(
      requires(element lhs, element rhs) {
        { lhs + rhs } -> std::convertible_to<element>;
      }, "One norm requires addable element types.");

  return cast<element, underlying>(tla::matrix_one_norm(value));
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_ONE_NORM_TPP
