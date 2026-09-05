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

#include "fcarouge/linalg.hpp"

#include <concepts>

namespace fcarouge::test {
namespace {
//! @test Verifies the same_as_typed_matrix concept holds at every level of a
//! typed matrix nested inside another typed matrix's `Matrix` parameter, not
//! only at the outermost level.
using outer = matrix<double, 2, 2>;      // typed_matrix<matrix1<...>, ...>
using middle = matrix1<double, 2, 2>;    // typed_matrix<matrix0<...>, ...>
using innermost = matrix0<double, 2, 2>; // typed_matrix<Eigen::Matrix, ...>

static_assert(same_as_typed_matrix<outer>);
static_assert(same_as_typed_matrix<middle>);
static_assert(same_as_typed_matrix<innermost>);

// Each level's `matrix` member type names the next level in, confirming the
// nesting is real rather than the aliases collapsing to one type.
static_assert(std::same_as<outer::matrix, middle>);
static_assert(std::same_as<middle::matrix, innermost>);
static_assert(same_as_typed_matrix<typename outer::matrix>);
static_assert(same_as_typed_matrix<typename middle::matrix>);

// The innermost `Matrix` parameter is the plain Eigen storage, not a typed
// matrix: nesting bottoms out.
static_assert(not same_as_typed_matrix<typename innermost::matrix>);
} // namespace
} // namespace fcarouge::test
