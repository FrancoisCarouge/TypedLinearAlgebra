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

#include <cassert>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the singleton by singleton matrix addition operator.
[[maybe_unused]] const auto test{[] {
  using au::symbols::m;

  constexpr auto m2{au::squared(m)};

  using length = au::QuantityD<au::Meters>;

  // Intended:
  // const row_vector<representation, length> a{3. * m};

  using area = au::QuantityD<au::UnitPowerT<au::Meters, 2>>;

  const row_vector<representation, area> a{3. * m2};

  const row_vector<representation, length> b{2. * m};

  const row_vector<representation, length> r{a + b};

  assert(5. * m == r.at());
  assert(5. * m == r[]);
  assert(5. * m == r());
  assert(5. * m == r);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
