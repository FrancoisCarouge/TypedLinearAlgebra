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
#include <tuple>

namespace fcarouge::test {
namespace {
using sparse::banded;
using sparse::dense;
using sparse::diagonal;
using sparse::lower_triangular;
using sparse::make_pattern;
using sparse::tridiagonal;
using sparse::upper_triangular;
using sparse::zero;

template <auto Pattern> using typed = sparse_matrix<double, Pattern>;
template <auto Pattern> using storage = sparse::matrix<double, Pattern>;

//! @test Verifies the sum and the difference store the union of their terms.
static_assert(
    std::same_as<decltype(typed<diagonal<3>>{} + typed<tridiagonal<3>>{}),
                 typed<tridiagonal<3>>>);
static_assert(
    std::same_as<decltype(typed<tridiagonal<3>>{} - typed<diagonal<3>>{}),
                 typed<tridiagonal<3>>>);
static_assert(std::same_as<decltype(typed<lower_triangular<3>>{} +
                                    typed<upper_triangular<3>>{}),
                           typed<dense<3, 3>>>);
static_assert(std::same_as<decltype(typed<zero<3, 3>>{} + typed<diagonal<3>>{}),
                           typed<diagonal<3>>>);

//! @test Verifies the product stores the structural product of its factors.
static_assert(
    std::same_as<decltype(typed<diagonal<3>>{} * typed<diagonal<3>>{}),
                 typed<diagonal<3>>>);
static_assert(
    std::same_as<decltype(typed<diagonal<3>>{} * typed<tridiagonal<3>>{}),
                 typed<tridiagonal<3>>>);
static_assert(
    std::same_as<decltype(typed<tridiagonal<4>>{} * typed<tridiagonal<4>>{}),
                 typed<banded<4, 4, 2, 2>>>);
static_assert(std::same_as<decltype(typed<lower_triangular<3>>{} *
                                    typed<lower_triangular<3>>{}),
                           typed<lower_triangular<3>>>);
static_assert(std::same_as<decltype(typed<zero<3, 3>>{} * typed<dense<3, 3>>{}),
                           typed<zero<3, 3>>>);

//! @test Verifies the transposition transposes the pattern, and the shape.
static_assert(std::same_as<decltype(transposed(typed<lower_triangular<3>>{})),
                           typed<upper_triangular<3>>>);

inline constexpr auto rectangular{make_pattern<2, 3>({{0, 2}, {1, 0}})};
inline constexpr auto rectangular_transposed{
    make_pattern<3, 2>({{0, 1}, {2, 0}})};

static_assert(std::same_as<decltype(transposed(typed<rectangular>{})),
                           typed<rectangular_transposed>>);

//! @test Verifies the negation and the scaling keep the pattern.
static_assert(
    std::same_as<decltype(-typed<tridiagonal<3>>{}), typed<tridiagonal<3>>>);
static_assert(std::same_as<decltype(typed<tridiagonal<3>>{} * 2.),
                           typed<tridiagonal<3>>>);
static_assert(std::same_as<decltype(2. * typed<tridiagonal<3>>{}),
                           typed<tridiagonal<3>>>);

//! @test Verifies a matrix-vector product is sparse when the factors are.
inline constexpr auto first_element{make_pattern<3, 1>({{0, 0}})};

static_assert(
    std::same_as<decltype(storage<diagonal<3>>{} * storage<first_element>{}),
                 storage<first_element>>);
static_assert(
    std::same_as<decltype(storage<tridiagonal<3>>{} * storage<first_element>{}),
                 storage<make_pattern<3, 1>({{0, 0}, {1, 0}})>>);

//! @test Verifies the element types are the ones of the indexes whatever the
//! pattern.
static_assert(std::same_as<typed<diagonal<3>>::element<0, 1>, double>);
static_assert(std::same_as<typed<diagonal<3>>::underlying, double>);
static_assert(
    std::same_as<sparse_matrix<float, diagonal<3>>::underlying, float>);

//! @test Verifies the indexes default to the representation type.
static_assert(
    std::same_as<typed<diagonal<2>>::row_indexes, std::tuple<double, double>>);
} // namespace
} // namespace fcarouge::test
