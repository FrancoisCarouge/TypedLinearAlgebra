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

#include <functional>
#include <tuple>

namespace fcarouge::test {
namespace {
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::s;
template <auto Reference>
using quantity = mp_units::quantity<Reference, double>;
using position = quantity<mp_units::isq::length[m]>;
using velocity = quantity<mp_units::isq::velocity[m / s]>;

using identity = std::tuple<std::identity>;
using positions = matrix<double, std::tuple<position, position>, identity>;
using state = matrix<double, std::tuple<position, velocity>, identity>;

//! @brief A requires-expression only yields `false` on substitution failure
//! within a template, hence the concept.
template <typename Type>
concept one_normable = requires(Type value) { matrix_one_norm(value); };

//! @test Verifies the one norm constraint is SFINAE-friendly: a
//! uniform matrix is accepted, a non-uniform one is excluded from overload
//! resolution rather than hard-erroring.
static_assert(one_normable<positions>);
static_assert(not one_normable<state>);
} // namespace
} // namespace fcarouge::test
