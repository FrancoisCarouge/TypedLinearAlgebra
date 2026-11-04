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
//! @test Verifies the storage compiles for a freestanding implementation, where
//! the compiler supports one.
#ifndef _MSC_VER
static_assert(__STDC_HOSTED__ == 0);
#endif

using band_matrix = sparse_matrix<double, sparse::tridiagonal<3>>;

inline constexpr band_matrix band{
    sparse::make_matrix<band_matrix>(1., 2., sparse::structural_zero, 3., 4.,
                                     5., sparse::structural_zero, 6., 7.)};

//! @test Verifies the operations without a hosted implementation, the
//! square root's included.
static_assert((band + band).at<0, 1>() == 4.);
static_assert((band - (band * 2.)).at<0, 1>() == -2.);
static_assert((-band).at<0, 1>() == -2.);
static_assert((band * band).at<1, 1>() == 52.);
static_assert((band * 2.).at<2, 2>() == 14.);
static_assert((band / 2.).at<0, 1>() == 1.);
static_assert(transposed(band).at<1, 0>() == 2.);
static_assert(matrix_one_norm(band) == 12.);
static_assert(matrix_inf_norm(band) == 13.);
static_assert(
    matrix_frob_norm(sparse::make_matrix<sparse_matrix<
                         double, sparse::make_pattern<2, 2>({{0, 1}, {1, 0}})>>(
        sparse::structural_zero, 3., 4., sparse::structural_zero)) == 5.);
} // namespace
} // namespace fcarouge::test
