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

#include <tuple>

namespace fcarouge::test {
namespace {
//! @test Verifies the same_shape concept over a built-in element type, Eigen
//! backend.

using singleton = matrix<double, 1, 1>;
using row = matrix<double, 1, 3>;
using column = matrix<double, 3, 1>;
using plane = matrix<double, 2, 3>;

// Positive: a typed matrix always has the same shape as itself, the relation
// is reflexive over every rank.
static_assert(same_shape<singleton, singleton>);
static_assert(same_shape<row, row>);
static_assert(same_shape<column, column>);
static_assert(same_shape<plane, plane>);

// Positive: the shape is the row and column counts only. The element types,
// and the backend element storage, play no part; two matrices of different
// element types but identical dimensions have the same shape.
static_assert(same_shape<matrix<double, 2, 3>, matrix<float, 2, 3>>);
static_assert(same_shape<matrix<double, 2, 3>, matrix<int, 2, 3>>);
static_assert(same_shape<row_vector<double, 4>, row_vector<float, 4>>);

// Positive: the relation is symmetric in its two arguments.
static_assert(same_shape<matrix<double, 2, 3>, matrix<float, 2, 3>> ==
              same_shape<matrix<float, 2, 3>, matrix<double, 2, 3>>);
static_assert(same_shape<row, column> == same_shape<column, row>);

// Positive, special case: cv-qualification and references are stripped with
// `std::remove_cvref_t` before the dimensions are compared.
static_assert(same_shape<const plane, plane>);
static_assert(same_shape<plane &, const plane &&>);
static_assert(same_shape<const plane &, volatile plane>);
static_assert(same_shape<const volatile plane &, plane>);

// Negative: a differing row count is a different shape.
static_assert(not same_shape<plane, matrix<double, 3, 3>>);
static_assert(not same_shape<singleton, column>);

// Negative: a differing column count is a different shape.
static_assert(not same_shape<plane, matrix<double, 2, 4>>);
static_assert(not same_shape<singleton, row>);

// Negative: a matrix and its transpose shape are not the same shape, even
// though they hold the same number of elements.
static_assert(not same_shape<plane, matrix<double, 3, 2>>);
static_assert(not same_shape<row, column>);

// Negative: a row vector and the column vector of the same length are both
// rank one with the same element count, but not the same shape.
static_assert(not same_shape<row_vector<double, 4>, column_vector<double, 4>>);

// Negative, special case: both arguments must be typed matrices. Any other
// type makes the concept unsatisfied through a substitution failure in the
// immediate context, never a hard error on a missing `::rows`/`::columns`.
static_assert(not same_shape<plane, double>);
static_assert(not same_shape<double, plane>);
static_assert(not same_shape<plane, void>);
static_assert(not same_shape<plane, std::tuple<double>>);
static_assert(not same_shape<int, float>);

// Negative, special case: the underlying backend matrix, and a pointer to a
// typed matrix, are not typed matrices themselves.
static_assert(not same_shape<plane, eigen::matrix<double, 2, 3>>);
static_assert(not same_shape<eigen::matrix<double, 2, 3>, plane>);
static_assert(not same_shape<plane, plane *>);

// Negative, special case: an incomplete type is probed safely, the concept is
// simply unsatisfied.
struct incomplete;
static_assert(not same_shape<plane, incomplete>);

// Agrees with the shape-carrying member counts and the shape concepts.
static_assert(same_shape<row, matrix<double, row::rows, row::columns>>);
static_assert(same_shape<row, row> and row_typed_matrix<row> and
              not column_typed_matrix<row>);
static_assert(same_shape<singleton, singleton> and
              rank_typed_matrix<singleton, 0>);
static_assert(not same_shape<row, column> and rank_typed_matrix<row, 1> and
              rank_typed_matrix<column, 1>);
} // namespace
} // namespace fcarouge::test
