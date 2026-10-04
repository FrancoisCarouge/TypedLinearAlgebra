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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_FROB_NORM_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_FROB_NORM_TPP

#if __has_include(<linalg>)

#include <linalg>

#endif

#include <concepts>
#include <type_traits>

namespace fcarouge {
[[nodiscard]] constexpr auto
matrix_frob_norm(const uniform_typed_matrix auto &value);

namespace typed_linear_algebra::internal {
//! @brief Frobenius norm of the backend storage, in the underlying type.
//!
//! @details Tries each backend's own spelling, most to least specific: an
//! Eigen `.norm()` member, a qualified recursion into `nested_typed_eigen`'s
//! composed typed matrix storage, an ADL-found Armadillo `norm(X, "fro")`,
//! then `std::linalg::matrix_frob_norm`.
[[nodiscard]] constexpr auto
matrix_frob_norm(const uniform_typed_matrix auto &value) {
  if constexpr (requires { value.data().norm(); }) {
    return value.data().norm();
  } else if constexpr (requires { fcarouge::matrix_frob_norm(value.data()); }) {
    return fcarouge::matrix_frob_norm(value.data());
  } else if constexpr (requires { norm(value.data(), "fro"); }) {
    return norm(value.data(), "fro");
  }
#if __has_include(<linalg>)
  else if constexpr (requires {
                       std::linalg::matrix_frob_norm(value.data());
                     }) {
    return std::linalg::matrix_frob_norm(value.data());
  }
#endif
  else {
    static_assert(
        sizeof(value) == 0,
        "Frobenius norm is not supported for this linear algebra backend.");
  }
}
} // namespace typed_linear_algebra::internal

//! @brief Frobenius norm of a matrix.
//!
//! @details The square root of the sum of the squares of the absolute values of
//! all the elements. For a row or column vector, the Euclidean L2 norm. For a
//! singleton, the absolute value. The elements must be of a uniform type
//! supporting addition, so that the norm shares the element's type: a matrix of
//! lengths has a length norm. Affine types, for example
//! `std::chrono::time_point`, are rejected.
//!
//! @param value The typed matrix.
//!
//! @return The Frobenius norm of `value`, of the element type.
//!
//! @see std::linalg::matrix_frob_norm
[[nodiscard]] constexpr auto
matrix_frob_norm(const uniform_typed_matrix auto &value) {
  using matrix = std::remove_cvref_t<decltype(value)>;
  using element = tla::element_at<matrix, 0, 0>;
  using underlying = typename matrix::underlying;

  static_assert(
      requires(element lhs, element rhs) {
        { lhs + rhs } -> std::convertible_to<element>;
      }, "Frobenius norm requires addable element types.");

  return cast<element, underlying>(tla::matrix_frob_norm(value));
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MATRIX_FROB_NORM_TPP
