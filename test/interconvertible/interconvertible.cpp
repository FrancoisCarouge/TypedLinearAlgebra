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

#include "fcarouge/typed_linear_algebra.hpp"

#include <type_traits>

namespace fcarouge::test {
namespace {
// A type is interconvertible with itself.
static_assert(tla::are_interconvertible<double, double>);
static_assert(not tla::are_not_interconvertible<double, double>);

// Unrelated types have no conversion in either direction.
static_assert(tla::are_not_interconvertible<double, double *>);
static_assert(tla::are_not_interconvertible<int, int *>);

// A conversion in a single direction is sufficient, and the relation is
// symmetric in the argument order.
static_assert(std::is_convertible_v<int *, void *>);
static_assert(not std::is_convertible_v<void *, int *>);
static_assert(tla::are_interconvertible<int *, void *>);
static_assert(tla::are_interconvertible<void *, int *>);

// Implicit numeric conversions count.
static_assert(tla::are_interconvertible<float, double>);
static_assert(tla::are_interconvertible<int, double>);
} // namespace
} // namespace fcarouge::test
