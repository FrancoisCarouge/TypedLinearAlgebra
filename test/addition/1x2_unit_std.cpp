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
using literals::operator""_i;
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

namespace {
//! @test Verifies the addition operator with non-trivial types.
[[maybe_unused]] const auto test{[] {
  using position = quantity<mp_units::isq::length[m]>;
  using velocity = quantity<mp_units::isq::velocity[m / s]>;

  double storage_a[]{0., 0.};
  double storage_b[]{0., 0.};
  double storage_r[]{0., 0.};

  std::mdspan span_a{&storage_a[0], std::extents<std::size_t, 1, 2>{}};
  std::mdspan span_b{&storage_b[0], std::extents<std::size_t, 1, 2>{}};
  std::mdspan span_r{&storage_r[0], std::extents<std::size_t, 1, 2>{}};

  row_vector<representation, position, velocity> a{span_a};
  row_vector<representation, position, velocity> b{span_b};
  row_vector<representation, position, velocity> r{span_r};

  a.at<0_i>(1. * m);
  a.at<1_i>(2. * m / s);
  b.at<0_i>(3. * m);
  b.at<1_i>(4. * m / s);

  add(a, b, r);

  assert(4. * m == r.at<0>());
  assert(4. * m == r.at<0_i>());
  assert(4. * m == r[0_i]);
  assert(4. * m == r(0_i));

  assert(6. * m / s == r.at<1>());
  assert(6. * m / s == r.at<1_i>());
  assert(6. * m / s == r[1_i]);
  assert(6. * m / s == r(1_i));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
