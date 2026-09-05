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

#include "fcarouge/eigen.hpp"
#include "fcarouge/linalg.hpp"

#include <string>
#include <tuple>

namespace fcarouge::test {
namespace {
//! @test Verifies `distinct_typed_matrix` over scalar elements, Eigen backend:
//! shape, cv-qualification, and non-typed-matrix rejection. A scalar matrix
//! carries one element type everywhere, so only its singleton is distinct.

// Singleton: no position pair to compare, distinct (and uniform) vacuously.
static_assert(distinct_typed_matrix<matrix<double, 1, 1>>);
static_assert(uniform_typed_matrix<matrix<double, 1, 1>>);

// More than one element: the sole type repeats, so never distinct.
static_assert(not distinct_typed_matrix<matrix<double, 1, 3>>);
static_assert(not distinct_typed_matrix<matrix<double, 3, 1>>);
static_assert(not distinct_typed_matrix<matrix<double, 2, 2>>);
static_assert(not distinct_typed_matrix<matrix<float, 4, 2>>);
static_assert(not distinct_typed_matrix<row_vector<double, 5>>);
static_assert(not distinct_typed_matrix<column_vector<double, 3>>);

// cv-qualifiers and references are stripped first.
using singleton = matrix<double, 1, 1>;
static_assert(distinct_typed_matrix<const singleton>);
static_assert(distinct_typed_matrix<volatile singleton>);
static_assert(distinct_typed_matrix<const singleton &>);
static_assert(distinct_typed_matrix<singleton &&>);

// Not a typed matrix: fails the leading conjunct, no hard error.
static_assert(not distinct_typed_matrix<void>);
static_assert(not distinct_typed_matrix<int>);
static_assert(not distinct_typed_matrix<double>);
static_assert(not distinct_typed_matrix<std::string>);
static_assert(not distinct_typed_matrix<std::tuple<double>>);
struct incomplete;
static_assert(not distinct_typed_matrix<incomplete>);
static_assert(not distinct_typed_matrix<eigen::matrix<double, 1, 1>>);
static_assert(not distinct_typed_matrix<matrix<double, 1, 1> *>);

// Independent from `uniform_typed_matrix`: they meet only on the singleton.
static_assert(distinct_typed_matrix<matrix<double, 1, 1>> ==
              uniform_typed_matrix<matrix<double, 1, 1>>);
static_assert(not distinct_typed_matrix<matrix<double, 2, 2>> and
              uniform_typed_matrix<matrix<double, 2, 2>>);
} // namespace
} // namespace fcarouge::test
