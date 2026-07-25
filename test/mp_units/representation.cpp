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

#include <cassert>
#include <format>
#include <print>
#include <utility>

namespace fcarouge::test {
namespace {
// Inject the magnitude customization-point object (CPO, a variable with a call
// operator) so it competes with the library's magnitude free function.
using mp_units::magnitude;
using namespace mp_units;

// Compile-time guarantee tests from:
// mp-unit's test/runtime/linear_algebra_test.cpp
using vec3 =
    typed_column_vector<Eigen::Vector<double, 3>, double, double, double>;

// The library vector type is accepted as a vector representation.
static_assert(RepresentationOf<vec3, quantity_tensor_order::vector>);

// A library scalar multiplication may yield a lazy expression template; the
// representation machinery must canonicalize it back to the concrete vector
// type rather than store the proxy (which would dangle). `decltype(vec3{}
// * 2.0)` is exactly such a proxy for Eigen/Blaze.
static_assert(
    std::same_as<
        representation_canonical_type_t<decltype(std::declval<vec3>() * 2.0)>,
        vec3>);

// Consequently, arithmetic on vector quantities stores the concrete vector
// type, not a proxy.
static_assert(std::same_as<decltype(vec3(1, 2, 3) * isq::velocity[m / s] *
                                    (2. * isq::duration[s]))::rep,
                           vec3>);

// A vector quantity is NOT a representation type: a quantity can never be
// nested as another quantity's representation (`value_type_t<quantity>` is
// the quantity itself, which `disable_representation` rejects).
static_assert(!detail::VectorRepresentation<decltype(vec3(0, 0, 0) *
                                                     isq::velocity[m / s])>);
static_assert(!RepresentationOf<decltype(vec3(0, 0, 0) * isq::velocity[m / s]),
                                quantity_tensor_order::vector>);

// `magnitude()` of a vector quantity is a scalar quantity in the same unit.
static_assert(
    QuantityOf<decltype(magnitude(vec3(3, 4, 0) * isq::velocity[m / s])),
               isq::speed>);
} // namespace
} // namespace fcarouge::test
