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

#include <array>
#include <concepts>
#include <type_traits>

namespace fcarouge::test {
namespace {
using sparse::diagonal;
using sparse::tridiagonal;

//! @test Verifies a matrix stores only its pattern's elements: there is no
//! index nor pointer of any kind in the object, the pattern is in the type.
static_assert(sizeof(sparse_matrix<double, diagonal<3>>) == 3 * sizeof(double));
static_assert(sizeof(sparse_matrix<double, tridiagonal<3>>) ==
              7 * sizeof(double));
static_assert(sizeof(sparse_matrix<float, tridiagonal<3>>) ==
              7 * sizeof(float));
static_assert(sizeof(sparse_matrix<double, sparse::dense<4, 4>>) ==
              16 * sizeof(double));
static_assert(sizeof(sparse_matrix<double, diagonal<16>>) ==
              16 * sizeof(double));

//! @test Verifies the memory saved over a dense matrix of the same shape.
static_assert(sizeof(sparse_matrix<double, diagonal<16>>) * 16 ==
              sizeof(sparse_matrix<double, sparse::dense<16, 16>>));

//! @test Verifies the storage has the properties of an array of values.
using storage = sparse::matrix<double, tridiagonal<3>>;

static_assert(std::is_trivially_copyable_v<storage>);
static_assert(std::is_standard_layout_v<storage>);
static_assert(std::is_trivially_destructible_v<storage>);
static_assert(std::is_nothrow_move_constructible_v<storage>);
static_assert(std::same_as<decltype(storage::values), std::array<double, 7>>);

//! @test Verifies the typed matrix adds nothing to its storage.
static_assert(
    std::is_trivially_copyable_v<sparse_matrix<double, tridiagonal<3>>>);

//! @test Verifies a matrix without structural zero is laid out as a dense
//! matrix, in row-major order.
inline constexpr sparse::matrix<double, sparse::dense<2, 2>> dense_values{
    1., 2., 3., 4.};

static_assert(dense_values.values[1] == 2.);
static_assert(dense_values[0, 1] == 2.);
static_assert(dense_values[1, 0] == 3.);
} // namespace
} // namespace fcarouge::test
