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

namespace fcarouge::test {
using literals::operator""_i;
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

namespace {
//! @test Verifies the singleton by singleton matrix multiplication operator.
[[maybe_unused]] const auto test{[] {
  using length = quantity<mp_units::isq::length[m]>;
  using area = quantity<mp_units::isq::area[m2]>;

  double storage_a{0.};
  double storage_b{0.};
  double storage_r{0.};

  std::mdspan span_a{&storage_a, std::extents<std::size_t, 1, 1>{}};
  std::mdspan span_b{&storage_b, std::extents<std::size_t, 1, 1>{}};
  std::mdspan span_r{&storage_r, std::extents<std::size_t, 1, 1>{}};

  row_vector<representation, length> a{span_a};
  row_vector<representation, length> b{span_b};
  row_vector<representation, area> r{span_r};

  a = 2. * m;
  b = 3. * m;
  r = a * b;

  assert(6. * m2 == r.at());
  assert(6. * m2 == r[]);
  assert(6. * m2 == r());
  assert(6. * m2 == r);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
