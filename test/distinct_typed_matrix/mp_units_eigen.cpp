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

#include <mp-units/systems/isq.h>
#include <mp-units/systems/si.h>

#include <functional>
#include <tuple>

namespace fcarouge::test {
namespace {
using mp_units::si::unit_symbols::km;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;

template <auto Reference>
using quantity = mp_units::quantity<Reference, double>;
using position = quantity<mp_units::isq::length[m]>;
using kilo_position = quantity<mp_units::isq::length[km]>;
using velocity = quantity<mp_units::isq::velocity[m / s]>;
using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;

template <typename Rows, typename Columns>
using shaped = matrix<double, Rows, Columns>;
using id = std::tuple<std::identity>;

//! @test Verifies `distinct_typed_matrix` over mp-units quantities: distinct
//! when no two element positions share an implicit conversion target.

// Singleton: distinct and uniform, vacuously.
static_assert(distinct_typed_matrix<shaped<std::tuple<position>, id>>);
static_assert(uniform_typed_matrix<shaped<std::tuple<position>, id>>);

// Column and row vectors of unrelated quantity kinds: distinct.
static_assert(
    distinct_typed_matrix<shaped<std::tuple<position, velocity>, id>>);
static_assert(
    distinct_typed_matrix<shaped<id, std::tuple<position, velocity>>>);

// Repeated element type: uniform, never distinct.
static_assert(
    not distinct_typed_matrix<shaped<std::tuple<velocity, velocity>, id>>);
static_assert(uniform_typed_matrix<shaped<std::tuple<velocity, velocity>, id>>);

// Two dimensions, four unrelated types, one of them the raw representation.
static_assert(
    distinct_typed_matrix<shaped<std::tuple<double, position>,
                                 std::tuple<velocity, acceleration>>>);

// Different but mutually convertible units (m, km): the case distinctness adds
// over non-uniformity. Neither distinct nor uniform.
static_assert(
    not distinct_typed_matrix<shaped<std::tuple<position, kilo_position>, id>>);
static_assert(
    not uniform_typed_matrix<shaped<std::tuple<position, kilo_position>, id>>);

// Wider than two: the non-adjacent convertible pair is still caught.
static_assert(distinct_typed_matrix<
              shaped<id, std::tuple<position, velocity, acceleration>>>);
static_assert(not distinct_typed_matrix<
              shaped<id, std::tuple<position, velocity, kilo_position>>>);

// Collision on positions sharing neither row nor column: {position, velocity} x
// {velocity, position} repeats length*speed at (0,0) and (1,1).
static_assert(
    not distinct_typed_matrix<shaped<std::tuple<position, velocity>,
                                     std::tuple<velocity, position>>>);

// Non-square, all six element types distinct.
static_assert(distinct_typed_matrix<
              shaped<std::tuple<double, position>,
                     std::tuple<velocity, acceleration, position>>>);

// cv-qualifiers and references are stripped first.
using mx = shaped<std::tuple<position, velocity>, id>;
static_assert(distinct_typed_matrix<const mx>);
static_assert(distinct_typed_matrix<mx &>);
static_assert(distinct_typed_matrix<const mx &&>);

// Not a typed matrix.
static_assert(not distinct_typed_matrix<void>);
static_assert(not distinct_typed_matrix<position>);
static_assert(not distinct_typed_matrix<std::tuple<position, velocity>>);

// Known over-approximation: a shared conversion target reachable only through a
// user-defined conversion, or a common base, is invisible; still reported
// distinct.
namespace shared_base {
struct base {};
struct first : base {};
struct second : base {};
static_assert(distinct_typed_matrix<shaped<std::tuple<first, second>, id>>);
} // namespace shared_base
} // namespace
} // namespace fcarouge::test
