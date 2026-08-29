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
#include <cstddef>
#include <mdspan>

namespace fcarouge::test {
namespace {
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;

//! @test Verifies the equal to algorithm for a column vector of
//! heterogeneous types with the mdspan backend.
[[maybe_unused]] const auto test{[] -> int {
  double storage_a[]{9., 10.};
  double storage_b[]{9., 10.};
  double storage_c[]{8., 10.};
  std::mdspan span_a{&storage_a[0], std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_b{&storage_b[0], std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_c{&storage_c[0], std::extents<std::size_t, 2, 1>{}};
  column_vector<double, decltype(1. * s), decltype(1. * s2)> a{span_a};
  column_vector<double, decltype(1. * s), decltype(1. * s2)> b{span_b};
  column_vector<double, decltype(1. * s), decltype(1. * s2)> c{span_c};

  assert(a == b);
  assert(a != c);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
