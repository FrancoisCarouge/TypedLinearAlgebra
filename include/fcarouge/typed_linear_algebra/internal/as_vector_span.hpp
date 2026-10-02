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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_AS_VECTOR_SPAN_HPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_AS_VECTOR_SPAN_HPP

//! @file
//! @brief Rank one view of a typed vector's storage, for `std::linalg`.
//!
//! @details Kept apart from `utility.hpp` so that only the backends providing
//! `std::linalg` pull in `<mdspan>`.

#if __has_include(<linalg>)

#include <cstddef>
#include <mdspan>
#include <type_traits>

namespace fcarouge::typed_linear_algebra::internal {
//! @brief Reinterprets a row or column typed vector's contiguous, rank two,
//! n-by-one or one-by-n storage as the rank one span required by
//! `std::linalg`'s vector concepts.
//!
//! @details A typed row or column vector, `rank_typed_matrix<1>`, is stored
//! as a rank two, n-by-one or one-by-n, underlying matrix, unlike the rank
//! one `in-vector`, `out-vector` shapes expected by `std::linalg`.
template <typename Type> constexpr auto as_vector_span(Type &value) {
  using matrix = std::remove_cvref_t<Type>;
  using underlying = typename matrix::underlying;

  return std::mdspan<underlying,
                     std::extents<std::size_t, matrix::rows * matrix::columns>>(
      value.data().data_handle());
}
} // namespace fcarouge::typed_linear_algebra::internal

#endif

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_AS_VECTOR_SPAN_HPP
