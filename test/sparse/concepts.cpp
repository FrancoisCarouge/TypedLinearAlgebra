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
#include <cstddef>
#include <type_traits>
#include <utility>

namespace fcarouge::test {
namespace {
using sparse::diagonal;
using sparse::make_pattern;

using diagonal_matrix = sparse_matrix<double, diagonal<2>>;
using column = sparse_matrix<double, make_pattern<2, 1>({{1, 0}})>;
using row = sparse_matrix<double, make_pattern<1, 2>({{0, 1}})>;
using singleton = sparse_matrix<double, sparse::dense<1, 1>>;
using empty_singleton = sparse_matrix<double, sparse::zero<1, 1>>;

template <typename Matrix, std::size_t Row, std::size_t Column>
concept writable_at =
    std::same_as<decltype(std::declval<typename Matrix::matrix &>()
                              [std::integral_constant<std::size_t, Row>{},
                               std::integral_constant<std::size_t, Column>{}]),
                 double &>;

template <typename Matrix, auto... Indexes>
concept readable_at =
    requires(const Matrix value) { value.template at<Indexes...>(); };

//! @test Verifies a typed sparse matrix satisfies the library's concepts,
//! according to its shape and its element types, whatever its pattern.
static_assert(same_as_typed_matrix<diagonal_matrix>);
static_assert(uniform_typed_matrix<diagonal_matrix>);
static_assert(rank_typed_matrix<diagonal_matrix, 2>);
static_assert(rank_typed_matrix<column, 1>);
static_assert(rank_typed_matrix<row, 1>);
static_assert(column_typed_matrix<column>);
static_assert(row_typed_matrix<row>);
static_assert(rank_typed_matrix<singleton, 0>);
static_assert(rank_typed_matrix<empty_singleton, 0>);
static_assert(
    same_shape<diagonal_matrix, sparse_matrix<double, sparse::dense<2, 2>>>);

//! @test Verifies every position can be read, a structural zero included, but
//! only a stored position has a reference for a write, at compile time.
static_assert(readable_at<diagonal_matrix, 0, 0>);
static_assert(readable_at<diagonal_matrix, 0, 1>);
static_assert(writable_at<diagonal_matrix, 0, 0>);
static_assert(writable_at<diagonal_matrix, 1, 1>);
static_assert(not writable_at<diagonal_matrix, 0, 1>);
static_assert(not writable_at<diagonal_matrix, 1, 0>);

//! @test Verifies the rule for a vector.
static_assert(readable_at<column, 0>);
static_assert(not writable_at<column, 0, 0>);
static_assert(writable_at<column, 1, 0>);
static_assert(readable_at<row, 0>);
static_assert(not writable_at<row, 0, 0>);
static_assert(writable_at<row, 0, 1>);

//! @test Verifies the rule for the singleton.
static_assert(writable_at<singleton, 0, 0>);
static_assert(readable_at<empty_singleton>);
static_assert(not writable_at<empty_singleton, 0, 0>);

//! @test Verifies a list construction, whose positions are only known at run
//! time, is only available without a structural zero, by the reference it
//! needs.
static_assert(
    std::same_as<
        decltype(std::declval<sparse::matrix<double, diagonal<2>> &>()[0, 1]),
        double>);
static_assert(
    std::same_as<decltype(std::declval<sparse::matrix<
                              double, sparse::dense<2, 2>> &>()[0, 1]),
                 double &>);

//! @test Verifies a by-value read of a non-const matrix does not yield an
//! assignable element, a structural zero included.
template <typename Matrix>
concept element_assignable =
    requires(Matrix value) { value.template at<0, 0>() = 1.; };

static_assert(not element_assignable<diagonal_matrix>);

//! @test Verifies a matrix converts when it cannot lose an element.
template <auto Pattern> using typed = sparse_matrix<double, Pattern>;

static_assert(std::convertible_to<sparse::matrix<double, diagonal<2>>,
                                  sparse::matrix<double, sparse::dense<2, 2>>>);
static_assert(
    not std::convertible_to<sparse::matrix<double, sparse::dense<2, 2>>,
                            sparse::matrix<double, diagonal<2>>>);
static_assert(
    not std::constructible_from<sparse::matrix<double, sparse::dense<2, 2>>,
                                sparse::matrix<float, diagonal<2>>>);
static_assert(
    std::convertible_to<typed<diagonal<2>>, typed<sparse::dense<2, 2>>>);
} // namespace
} // namespace fcarouge::test
