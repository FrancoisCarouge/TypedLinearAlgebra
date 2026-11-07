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

namespace fcarouge::test {
namespace {
using sparse::diagonal;
using sparse::lower_triangular;
using sparse::make_pattern;
using sparse::pattern_intersection;
using sparse::pattern_product;
using sparse::pattern_transpose;
using sparse::pattern_union;
using sparse::tridiagonal;
using sparse::upper_triangular;

inline constexpr auto rectangular{make_pattern<2, 3>({{0, 2}, {1, 0}})};
inline constexpr auto first_row{make_pattern<3, 3>({{0, 0}, {0, 1}, {0, 2}})};
inline constexpr auto last_column{make_pattern<3, 3>({{0, 2}, {1, 2}, {2, 2}})};

//! @test Verifies the union, and that it is commutative and idempotent.
static_assert(pattern_union<first_row, last_column>.count == 5);
static_assert(pattern_union<first_row, last_column> ==
              pattern_union<last_column, first_row>);
static_assert(pattern_union<first_row, first_row> == first_row);
static_assert(pattern_union<first_row, sparse::zero<3, 3>> == first_row);
static_assert(pattern_union<diagonal<3>, tridiagonal<3>> == tridiagonal<3>);
static_assert(pattern_union<lower_triangular<3>, upper_triangular<3>> ==
              sparse::dense<3, 3>);

//! @test Verifies the intersection.
static_assert(pattern_intersection<first_row, last_column> ==
              make_pattern<3, 3>({{0, 2}}));
static_assert(pattern_intersection<first_row, sparse::zero<3, 3>>.count == 0);
static_assert(pattern_intersection<lower_triangular<3>, upper_triangular<3>> ==
              diagonal<3>);

//! @test Verifies the transposition, its shape, and that it is an involution
//! that keeps the row-major order.
static_assert(pattern_transpose<rectangular> ==
              make_pattern<3, 2>({{0, 1}, {2, 0}}));
static_assert(pattern_transpose<pattern_transpose<rectangular>> == rectangular);
static_assert(pattern_transpose<lower_triangular<3>> == upper_triangular<3>);
static_assert(pattern_transpose<diagonal<3>> == diagonal<3>);
static_assert(pattern_transpose<first_row> ==
              make_pattern<3, 3>({{0, 0}, {1, 0}, {2, 0}}));

//! @test Verifies the structural product of well-known families.
static_assert(pattern_product<diagonal<3>, diagonal<3>> == diagonal<3>);
static_assert(pattern_product<diagonal<3>, tridiagonal<3>> == tridiagonal<3>);
static_assert(pattern_product<tridiagonal<3>, diagonal<3>> == tridiagonal<3>);
static_assert(pattern_product<tridiagonal<4>, tridiagonal<4>> ==
              sparse::banded<4, 4, 2, 2>);
static_assert(pattern_product<lower_triangular<3>, lower_triangular<3>> ==
              lower_triangular<3>);
static_assert(pattern_product<upper_triangular<3>, upper_triangular<3>> ==
              upper_triangular<3>);
static_assert(pattern_product<lower_triangular<3>, upper_triangular<3>> ==
              sparse::dense<3, 3>);

//! @test Verifies a structural zero annihilates the product.
static_assert(pattern_product<sparse::zero<3, 3>, sparse::dense<3, 3>>.count ==
              0);
static_assert(pattern_product<sparse::dense<3, 3>, sparse::zero<3, 3>>.count ==
              0);

//! @test Verifies the product of two non-empty patterns that never meet.
inline constexpr auto first_column{make_pattern<2, 2>({{0, 0}, {1, 0}})};
inline constexpr auto second_row{make_pattern<2, 2>({{1, 0}, {1, 1}})};

static_assert(pattern_product<first_column, second_row>.count == 0);

//! @test Verifies the product's shape.
static_assert(pattern_product<rectangular, sparse::dense<3, 4>>.rows == 2);
static_assert(pattern_product<rectangular, sparse::dense<3, 4>>.columns == 4);

//! @test Verifies transposition reverses the product.
static_assert(pattern_transpose<pattern_product<first_row, last_column>> ==
              pattern_product<pattern_transpose<last_column>,
                              pattern_transpose<first_row>>);

//! @test Verifies the count of elementary products, not of result positions.
static_assert(sparse::product_term_count<diagonal<3>, diagonal<3>> == 3);
static_assert(sparse::product_term_count<tridiagonal<3>, tridiagonal<3>> == 17);
static_assert(
    sparse::product_term_count<sparse::dense<3, 3>, sparse::dense<3, 3>> == 27);

//! @test Verifies the subset relation.
static_assert(sparse::is_subset_of<diagonal<3>, tridiagonal<3>>);
static_assert(not sparse::is_subset_of<tridiagonal<3>, diagonal<3>>);
static_assert(sparse::is_subset_of<sparse::zero<3, 3>, diagonal<3>>);
static_assert(sparse::is_subset_of<diagonal<3>, diagonal<3>>);
static_assert(not sparse::is_subset_of<diagonal<3>, diagonal<2>>);
} // namespace
} // namespace fcarouge::test
