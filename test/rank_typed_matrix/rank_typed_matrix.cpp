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

#include <cstddef>
#include <tuple>

namespace fcarouge::test {
namespace {
//! @test Verifies the rank_typed_matrix concept over a built-in element type,
//! Eigen backend.

using singleton = matrix<double, 1, 1>;
using row = matrix<double, 1, 6>;
using column = matrix<double, 6, 1>;
using plane = matrix<double, 3, 4>;

// Singleton is rank 0.
static_assert(rank_typed_matrix<singleton, 0>);
static_assert(not rank_typed_matrix<singleton, 1>);
static_assert(not rank_typed_matrix<singleton, 2>);

// Row and column vectors are rank 1.
static_assert(rank_typed_matrix<row, 1>);
static_assert(rank_typed_matrix<column, 1>);
static_assert(not rank_typed_matrix<row, 0>);
static_assert(not rank_typed_matrix<row, 2>);
static_assert(not rank_typed_matrix<column, 0>);
static_assert(not rank_typed_matrix<column, 2>);

// Two-dimension matrices are rank 2.
static_assert(rank_typed_matrix<plane, 2>);
static_assert(not rank_typed_matrix<plane, 0>);
static_assert(not rank_typed_matrix<plane, 1>);

// A one-by-one is a singleton, not a rank 1 vector.
static_assert(rank_typed_matrix<matrix<double, 1, 1>, 0>);
static_assert(rank_typed_matrix<row_vector<double, 4>, 1>);
static_assert(rank_typed_matrix<column_vector<double, 4>, 1>);
static_assert(rank_typed_matrix<row_vector<double, 1>, 0>);
static_assert(rank_typed_matrix<column_vector<double, 1>, 0>);

// Ranks outside zero, one, two never match.
static_assert(not rank_typed_matrix<singleton, 3>);
static_assert(not rank_typed_matrix<plane, 3>);
static_assert(not rank_typed_matrix<plane, -1>);

// The rank argument accepts any integral spelling.
static_assert(rank_typed_matrix<singleton, 0U>);
static_assert(rank_typed_matrix<plane, 2L>);
static_assert(rank_typed_matrix<plane, std::size_t{2}>);

// Reference and cv-qualification are transparent.
static_assert(rank_typed_matrix<const plane, 2>);
static_assert(rank_typed_matrix<plane &, 2>);
static_assert(rank_typed_matrix<const plane &&, 2>);

// Non typed matrix types never match, no hard error.
static_assert(not rank_typed_matrix<double, 0>);
static_assert(not rank_typed_matrix<eigen::matrix<double, 1, 1>, 0>);
static_assert(not rank_typed_matrix<std::tuple<double>, 1>);
static_assert(not rank_typed_matrix<int[3], 1>);

// Agrees with the member rank and with the shape concepts.
static_assert(rank_typed_matrix<singleton, singleton::rank>);
static_assert(rank_typed_matrix<row, row::rank>);
static_assert(rank_typed_matrix<plane, plane::rank>);
static_assert(rank_typed_matrix<singleton, 0> and
              row_typed_matrix<singleton> and column_typed_matrix<singleton>);
static_assert(rank_typed_matrix<row, 1> and row_typed_matrix<row> and
              not column_typed_matrix<row>);
static_assert(rank_typed_matrix<column, 1> and column_typed_matrix<column> and
              not row_typed_matrix<column>);
} // namespace
} // namespace fcarouge::test
