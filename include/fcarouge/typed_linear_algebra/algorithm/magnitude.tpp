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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MAGNITUDE_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MAGNITUDE_TPP

namespace fcarouge {
[[nodiscard]] constexpr auto
matrix_frob_norm(const uniform_typed_matrix auto &value);

//! @brief Euclidean L2 norm of a row or column vector.
//!
//! @details The Frobenius norm of a vector, under its geometric name.
//! Restricted to vectors, where the norm has this conventional meaning.
[[nodiscard]] constexpr auto magnitude(const uniform_typed_matrix auto &value) {
  static_assert(
      rank_typed_matrix<decltype(value), 1>,
      "The magnitude operation only supports vector types at this time.");

  return fcarouge::matrix_frob_norm(value);
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_ALGORITHM_MAGNITUDE_TPP
