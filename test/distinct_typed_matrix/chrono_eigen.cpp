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

#include <chrono>
#include <functional>
#include <ratio>
#include <tuple>

namespace fcarouge::test {
namespace {
using seconds = std::chrono::duration<double>;
using milliseconds = std::chrono::duration<double, std::milli>;
using minutes = std::chrono::duration<double, std::ratio<60>>;
using timestamp = std::chrono::time_point<std::chrono::steady_clock, seconds>;

template <typename Rows, typename Columns>
using shaped = matrix<double, Rows, Columns>;
using id = std::tuple<std::identity>;

//! @test Verifies `distinct_typed_matrix` over std::chrono index types, Eigen
//! backend.

// Singleton: one position, distinct (and uniform) vacuously.
static_assert(distinct_typed_matrix<shaped<std::tuple<seconds>, id>>);
static_assert(uniform_typed_matrix<shaped<std::tuple<seconds>, id>>);

// A duration and a time point neither convert to one another nor have a common
// type: distinct.
static_assert(
    distinct_typed_matrix<shaped<std::tuple<seconds, timestamp>, id>>);
static_assert(
    distinct_typed_matrix<shaped<id, std::tuple<timestamp, minutes>>>);

// Same duration type repeated: uniform, not distinct.
static_assert(
    not distinct_typed_matrix<shaped<std::tuple<seconds, seconds>, id>>);
static_assert(uniform_typed_matrix<shaped<std::tuple<seconds, seconds>, id>>);

// Any two duration types share a std::common_type (a duration of their common
// period), so two duration positions are never distinct - near periods (s, ms)
// or far (s, min) alike.
static_assert(
    not distinct_typed_matrix<shaped<std::tuple<seconds, milliseconds>, id>>);
static_assert(
    not uniform_typed_matrix<shaped<std::tuple<seconds, milliseconds>, id>>);
static_assert(
    not distinct_typed_matrix<shaped<id, std::tuple<seconds, minutes>>>);
} // namespace
} // namespace fcarouge::test
