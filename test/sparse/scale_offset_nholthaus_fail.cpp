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

#include <units/temperature.h>

#include <cassert>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies a conversion into a structural zero does not compile when it
//! would change its value: zero degrees Celsius is 273.15 kelvin.
[[maybe_unused]] const auto test{[] {
  using celsius = units::temperature::celsius<representation>;
  using kelvin = units::temperature::kelvin<representation>;

  constexpr auto top{sparse::make_pattern<2, 1>({{0, 0}})};

  using dense_celsius =
      sparse_column_vector<representation, sparse::dense<2, 1>, celsius,
                           celsius>;
  using dense_kelvin =
      sparse_column_vector<representation, sparse::dense<2, 1>, kelvin, kelvin>;
  using sparse_celsius =
      sparse_column_vector<representation, top, celsius, celsius>;
  using sparse_kelvin =
      sparse_column_vector<representation, top, kelvin, kelvin>;

  sparse_celsius value;

  value.at<0>(celsius{1.});

  // Intended:
  // const dense_kelvin converted{dense_celsius{value}};
  const sparse_kelvin converted{value};
  assert(converted.at<0>() == kelvin{274.15});

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
