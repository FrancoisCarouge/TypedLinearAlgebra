/* Typed Linear Algebra
Version 0.2.0
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

namespace fcarouge::test {
using literals::operator""_i;
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

namespace {
//! @test Verifies the scalar multiplication operator.
[[maybe_unused]] const auto test{[] {
  using length = quantity<mp_units::isq::length[m]>;
  using area = quantity<mp_units::isq::area[m2]>;

  const matrix<representation, std::tuple<length, length>,
               std::tuple<length, length>>
      a{{1. * m2, 2. * m2}, {3. * m2, 4. * m2}};
  const matrix<representation, std::tuple<length, length>,
               std::tuple<area, area>>
      r{a * 2. * m};

  assert(r(0_i, 0_i) == 2. * m3);
  assert(r(0_i, 1_i) == 4. * m3);
  assert(r(1_i, 0_i) == 6. * m3);
  assert(r(1_i, 1_i) == 8. * m3);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
