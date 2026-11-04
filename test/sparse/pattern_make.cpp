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

namespace fcarouge::test {
namespace {
using sparse::make_pattern;
using sparse::position;

inline constexpr auto unordered{
    make_pattern<3, 4>({{2, 3}, {0, 1}, {1, 2}, {0, 0}})};

//! @test Verifies the positions are sorted in row-major order, whatever the
//! order they are given in.
static_assert(unordered.positions[0] == position{0, 0});
static_assert(unordered.positions[1] == position{0, 1});
static_assert(unordered.positions[2] == position{1, 2});
static_assert(unordered.positions[3] == position{2, 3});

//! @test Verifies the pattern's shape and count.
static_assert(unordered.rows == 3);
static_assert(unordered.columns == 4);
static_assert(unordered.count == 4);

//! @test Verifies a position's index is its row-major rank, and a structural
//! zero has none.
static_assert(unordered.index_of(0, 0) == 0);
static_assert(unordered.index_of(0, 1) == 1);
static_assert(unordered.index_of(1, 2) == 2);
static_assert(unordered.index_of(2, 3) == 3);
static_assert(unordered.index_of(1, 0) == sparse::absent);
static_assert(unordered.index_of(2, 2) == sparse::absent);
static_assert(unordered.stores(2, 3));
static_assert(not unordered.stores(0, 3));

//! @test Verifies a lookup past the last position, or out of the matrix, finds
//! no element.
static_assert(unordered.index_of(2, 4) == sparse::absent);
static_assert(unordered.index_of(3, 0) == sparse::absent);

//! @test Verifies the order of the given positions does not change the type or
//! the value: two equal patterns are the same matrix type.
inline constexpr auto ascending{make_pattern<2, 2>({{0, 0}, {1, 1}})};
inline constexpr auto descending{make_pattern<2, 2>({{1, 1}, {0, 0}})};

static_assert(ascending == descending);
static_assert(std::same_as<sparse::matrix<double, ascending>,
                           sparse::matrix<double, descending>>);

//! @test Verifies an empty pattern stores nothing.
inline constexpr auto empty{make_pattern<2, 3>()};

static_assert(empty.count == 0);
static_assert(not empty.stores(0, 0));
static_assert(empty.index_of(0, 0) == sparse::absent);

//! @test Verifies a pattern generated from a predicate.
inline constexpr auto is_even{[](std::size_t row, std::size_t column) -> bool {
  return (row + column) % 2 == 0;
}};
inline constexpr auto checkerboard{sparse::generate<3, 3, is_even>()};

static_assert(checkerboard.count == 5);
static_assert(checkerboard.stores(1, 1));
static_assert(not checkerboard.stores(0, 1));

//! @test Verifies the generated and the listed patterns of the same positions
//! are the same matrix type.
static_assert(checkerboard ==
              make_pattern<3, 3>({{0, 0}, {0, 2}, {1, 1}, {2, 0}, {2, 2}}));
} // namespace
} // namespace fcarouge::test
