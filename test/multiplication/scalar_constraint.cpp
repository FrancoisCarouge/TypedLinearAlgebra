/* Typed Linear Algebra
Version 0.5.0
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

namespace fcarouge::test {
namespace {
//! @brief A type without multiplication.
struct unrelated {};

using matrix1x1 = matrix<double, 1, 1>;
using matrix2x2 = matrix<double, 2, 2>;

//! @brief A requires-expression only yields `false` on substitution failure
//! within a template, hence the concept.
template <typename Lhs, typename Rhs>
concept multipliable_by = requires(Lhs lhs, Rhs rhs) { lhs * rhs; };

//! @test Verifies the scalar multiplication constraint is SFINAE-friendly, for
//! a singleton and a matrix: a scale multiplying the row indexes is accepted,
//! on either side, and any other type is excluded from overload resolution
//! rather than hard-erroring.
static_assert(multipliable_by<matrix1x1, double>);
static_assert(multipliable_by<double, matrix1x1>);
static_assert(not multipliable_by<matrix1x1, unrelated>);
static_assert(not multipliable_by<unrelated, matrix1x1>);
static_assert(multipliable_by<matrix2x2, double>);
static_assert(multipliable_by<double, matrix2x2>);
static_assert(not multipliable_by<matrix2x2, unrelated>);
static_assert(not multipliable_by<unrelated, matrix2x2>);
} // namespace
} // namespace fcarouge::test
