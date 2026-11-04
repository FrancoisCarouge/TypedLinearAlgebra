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
of the public at large_unit and to the detriment of our heirs and
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

#include <au/prefix.hh>
#include <au/units/meters.hh>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the dot product rejects terms of different scales, meters and
//! kilometers.
[[maybe_unused]] const auto test{[] -> int {
  using small_unit = au::QuantityD<au::Meters>;
  using large_unit = au::QuantityD<au::Kilo<au::Meters>>;
  using mixed_column = column_vector<representation, small_unit, large_unit>;
  using small_column = column_vector<representation, small_unit, small_unit>;
  using scalars = column_vector<representation, representation, representation>;

  const small_unit first{au::meters(500.)};
  const large_unit second{au::kilo(au::meters)(3.)};
  const mixed_column values{first, second};
  const scalars ones{1., 1.};

  // Intended: [[maybe_unused]] const auto result{dot(small_column{values},
  // ones)};
  [[maybe_unused]] const auto result{dot(values, ones)};

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
