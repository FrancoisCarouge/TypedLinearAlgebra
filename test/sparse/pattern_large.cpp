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

#include <concepts>
#include <cstddef>

namespace fcarouge::test {
namespace {
//! @test Verifies the compile-time cost is practical at the supported scale, a
//! few hundred stored elements, for each pattern operation.
inline constexpr std::size_t size{100};
inline constexpr auto wide{sparse::tridiagonal<size>};

static_assert(wide.count == 298);
static_assert(sparse::pattern_transpose<wide> == wide);
static_assert(sparse::pattern_union<wide, sparse::diagonal<size>> == wide);
static_assert(sparse::pattern_intersection<wide, sparse::diagonal<size>> ==
              sparse::diagonal<size>);
static_assert(sparse::product_term_count<wide, wide> == 890);
static_assert(sparse::pattern_product<wide, wide>.count == 494);

//! @test Verifies a matrix of that scale only stores its elements.
using matrix_type = sparse_matrix<double, wide>;

static_assert(sizeof(matrix_type) == 298 * sizeof(double));

//! @test Verifies the operations on a matrix of that scale in constant
//! evaluation: a tridiagonal squared is pentadiagonal, and the identity is
//! neutral.
inline constexpr auto square{[] -> auto {
  sparse::matrix<double, wide> value;

  for (auto &element : value.values) {
    element = 2.;
  }

  return value * value;
}()};

static_assert(decltype(square)::structure == sparse::banded<size, size, 2, 2>);
static_assert(square[0, 0] == 8.);
static_assert(square[1, 1] == 12.);
static_assert(square[0, 2] == 4.);
static_assert(square[0, 3] == 0.);
static_assert(square[99, 99] == 8.);

inline constexpr auto identity{[] -> auto {
  sparse::matrix<double, sparse::diagonal<size>> value;

  for (auto &element : value.values) {
    element = 1.;
  }

  return value;
}()};

static_assert(
    std::same_as<decltype(identity * square),
                 sparse::matrix<double, decltype(square)::structure>>);
static_assert(identity * square == square);
static_assert(square * identity == square);
} // namespace
} // namespace fcarouge::test
