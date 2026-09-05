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

#include "fcarouge/linalg.hpp"

#include <au/units/meters.hh>
#include <au/units/seconds.hh>

#include <functional>
#include <tuple>

namespace fcarouge::test {
namespace {
using position = au::QuantityD<au::Meters>;
using kilo_position = au::QuantityD<au::Kilo<au::Meters>>;
using velocity = au::QuantityD<au::UnitQuotientT<au::Meters, au::Seconds>>;

template <typename Rows, typename Columns>
using shaped = matrix<double, Rows, Columns>;
using id = std::tuple<std::identity>;

//! @test Verifies `distinct_typed_matrix` over Au quantity types, Eigen
//! backend.

// Unrelated dimensions: distinct.
static_assert(
    distinct_typed_matrix<shaped<std::tuple<position, velocity>, id>>);

// Off-diagonal collision: {position, velocity} x {velocity, position} repeats
// area/time at (0,0) and (1,1). Not distinct.
static_assert(
    not distinct_typed_matrix<shaped<std::tuple<position, velocity>,
                                     std::tuple<velocity, position>>>);

// Same type repeated: uniform, not distinct.
static_assert(
    not distinct_typed_matrix<shaped<std::tuple<position, position>, id>>);
static_assert(uniform_typed_matrix<shaped<std::tuple<position, position>, id>>);

// Same dimension, different unit, mutually convertible (m, km): not distinct,
// not uniform.
static_assert(
    not distinct_typed_matrix<shaped<std::tuple<position, kilo_position>, id>>);
static_assert(
    not uniform_typed_matrix<shaped<std::tuple<position, kilo_position>, id>>);

// A raw representation next to a quantity stays distinct: Au requires an
// explicit conversion either way.
static_assert(distinct_typed_matrix<shaped<std::tuple<double, position>, id>>);
} // namespace
} // namespace fcarouge::test
