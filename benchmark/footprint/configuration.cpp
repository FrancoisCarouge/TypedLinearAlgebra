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

#include "fcarouge/linalg.hpp"

#include <cstddef>

namespace fcarouge::benchmark {
//! @brief The number of rows of the matrices.
inline constexpr std::size_t rows{${ROWS}};

//! @brief The number of columns of the matrices.
inline constexpr std::size_t columns{${COLUMNS}};

//! @brief The typed matrix and its factor of the transposed shape.
using typed = matrix<float, rows, columns>;
using typed_factor = matrix<float, columns, rows>;

//! @brief The freestanding storage of the typed matrix and of its factor.
using storage = freestanding::matrix<float, rows, columns>;
using storage_factor = freestanding::matrix<float, columns, rows>;
} // namespace fcarouge::benchmark

//! @brief The entry points of the operations, each linked alone, of the typed
//! matrix and of its freestanding storage.
//!
//! @details The operands are the caller's. A matrix result is statically
//! allocated, a scalar result returned. The equality has no storage operation:
//! the library implements it over the typed elements.
extern "C" {
using fcarouge::benchmark::storage;
using fcarouge::benchmark::storage_factor;
using fcarouge::benchmark::typed;
using fcarouge::benchmark::typed_factor;

[[gnu::flatten]] auto typed_addition_step(const typed &lhs, const typed &rhs)
    -> const decltype(lhs + rhs) * {
  static decltype(lhs + rhs) result;
  result = lhs + rhs;
  return &result;
}

[[gnu::flatten]] auto storage_addition_step(const storage &lhs,
                                            const storage &rhs) -> const
    decltype(lhs + rhs) * {
  static decltype(lhs + rhs) result;
  result = lhs + rhs;
  return &result;
}

[[gnu::flatten]] auto typed_subtraction_step(const typed &lhs, const typed &rhs)
    -> const decltype(lhs - rhs) * {
  static decltype(lhs - rhs) result;
  result = lhs - rhs;
  return &result;
}

[[gnu::flatten]] auto storage_subtraction_step(const storage &lhs,
                                               const storage &rhs) -> const
    decltype(lhs - rhs) * {
  static decltype(lhs - rhs) result;
  result = lhs - rhs;
  return &result;
}

[[gnu::flatten]] auto typed_negation_step(const typed &value) -> const
    decltype(-value) * {
  static decltype(-value) result;
  result = -value;
  return &result;
}

[[gnu::flatten]] auto storage_negation_step(const storage &value) -> const
    decltype(-value) * {
  static decltype(-value) result;
  result = -value;
  return &result;
}

[[gnu::flatten]] auto typed_product_step(const typed &lhs,
                                         const typed_factor &rhs) -> const
    decltype(lhs * rhs) * {
  static decltype(lhs * rhs) result;
  result = lhs * rhs;
  return &result;
}

[[gnu::flatten]] auto storage_product_step(const storage &lhs,
                                           const storage_factor &rhs) -> const
    decltype(lhs * rhs) * {
  static decltype(lhs * rhs) result;
  result = lhs * rhs;
  return &result;
}

[[gnu::flatten]] auto typed_scale_step(const typed &lhs, float rhs) -> const
    decltype(lhs * rhs) * {
  static decltype(lhs * rhs) result;
  result = lhs * rhs;
  return &result;
}

[[gnu::flatten]] auto storage_scale_step(const storage &lhs, float rhs) -> const
    decltype(lhs * rhs) * {
  static decltype(lhs * rhs) result;
  result = lhs * rhs;
  return &result;
}

[[gnu::flatten]] auto typed_left_scale_step(float lhs, const typed &rhs)
    -> const decltype(lhs * rhs) * {
  static decltype(lhs * rhs) result;
  result = lhs * rhs;
  return &result;
}

[[gnu::flatten]] auto storage_left_scale_step(float lhs, const storage &rhs)
    -> const decltype(lhs * rhs) * {
  static decltype(lhs * rhs) result;
  result = lhs * rhs;
  return &result;
}

[[gnu::flatten]] auto typed_equal_to_step(const typed &lhs, const typed &rhs)
    -> bool {
  return lhs == rhs;
}
}
