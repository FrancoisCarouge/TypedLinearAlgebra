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

#include <units/length.h>

#include <cassert>
#include <functional>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the nholthaus/units quantities operations this library
//! already computes correctly on meters and kilometers: the elements, equality,
//! singletons, and uniform units. Contrasts the other tests in this directory.
[[maybe_unused]] const auto test{[] -> int {
  using small_unit = units::length::meters<representation>;
  using large_unit = units::length::kilometers<representation>;
  using small_column = column_vector<representation, small_unit, small_unit>;
  using large_column = column_vector<representation, large_unit, large_unit>;
  using small_singleton = column_vector<representation, small_unit>;
  using large_singleton = column_vector<representation, large_unit>;
  using gain =
      matrix<representation, std::tuple<representation, representation>,
             std::tuple<std::identity, std::identity>>;

  const small_unit first{500.};
  const large_unit second{3.};

  // Each element keeps the value it was constructed from.
  const column_vector<representation, small_unit, large_unit> mixed{first,
                                                                    second};

  assert(mixed.at<0>() == first);
  assert(mixed.at<1>() == second);

  // Equality compares the values, not the numbers.
  const small_column small_ones{small_unit{1000.}, small_unit{1000.}};
  const large_column large_ones{large_unit{1.}, large_unit{1.}};

  assert(small_ones == large_ones);

  // A singleton addition and subtraction compute on the values.
  const small_singleton small_singleton_value{first};
  const large_singleton large_singleton_value{second};

  assert(small_singleton_value + large_singleton_value == first + second);
  assert(small_singleton_value - large_singleton_value == first - second);

  // Uniform units need no scaling.
  const gain ones{{1., 1.}, {1., 1.}};
  const small_column uniform_product{ones * small_column{first, first}};

  assert(uniform_product.at<0>() == first + first);

  const large_column uniform_sum{large_column{second, second} +
                                 large_column{second, second}};

  assert(uniform_sum.at<0>() == second + second);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
