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
//! @test Verifies the counts of the named patterns.
static_assert(sparse::dense<2, 3>.count == 6);
static_assert(sparse::dense<1, 1>.count == 1);
static_assert(sparse::zero<2, 3>.count == 0);
static_assert(sparse::diagonal<4>.count == 4);
static_assert(sparse::tridiagonal<4>.count == 10);
static_assert(sparse::tridiagonal<1>.count == 1);
static_assert(sparse::lower_triangular<4>.count == 10);
static_assert(sparse::upper_triangular<4>.count == 10);

//! @test Verifies the shape of a rectangular band.
static_assert(sparse::banded<2, 4, 0, 1>.count == 4);
static_assert(sparse::banded<4, 2, 1, 0>.count == 4);

//! @test Verifies the positions of the named patterns.
static_assert(sparse::diagonal<3>.stores(1, 1));
static_assert(not sparse::diagonal<3>.stores(0, 1));
static_assert(sparse::tridiagonal<3>.stores(0, 1));
static_assert(sparse::tridiagonal<3>.stores(1, 0));
static_assert(not sparse::tridiagonal<3>.stores(0, 2));
static_assert(sparse::lower_triangular<3>.stores(2, 0));
static_assert(not sparse::lower_triangular<3>.stores(0, 2));
static_assert(sparse::upper_triangular<3>.stores(0, 2));
static_assert(not sparse::upper_triangular<3>.stores(2, 0));

//! @test Verifies the named patterns agree where they coincide.
static_assert(sparse::banded<3, 3, 2, 2> == sparse::dense<3, 3>);
static_assert(sparse::banded<3, 3, 0, 0> == sparse::diagonal<3>);
static_assert(sparse::diagonal<1> == sparse::dense<1, 1>);
} // namespace
} // namespace fcarouge::test
