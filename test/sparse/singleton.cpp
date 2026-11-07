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

#include <cassert>

namespace fcarouge::test {
namespace {
//! @test Verifies the singleton matrix, stored or a structural zero.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  sparse_matrix<double, sparse::dense<1, 1>> value{2.};

  assert(value.at() == 2.);
  assert(value[] == 2.);
  assert(value() == 2.);
  assert(value == 2.);

  value.at(3.);

  assert(value.at() == 3.);

  value = 4.;

  assert(value.at() == 4.);

  [[maybe_unused]] const double converted{value};

  assert(converted == 4.);

  // Arithmetic is on the element.
  const sparse_matrix<double, sparse::dense<1, 1>> other{5.};

  assert(value + other == 9.);
  assert(value - other == -1.);
  assert(value * other == 20.);
  assert(value / other == 0.8);
  assert(value * 2. == 8.);
  assert(2. * value == 8.);
  assert(value + 1. == 5.);
  assert(1. + value == 5.);

  // A structural zero singleton reads, and operates, as a zero.
  const sparse_matrix<double, sparse::zero<1, 1>> none;

  assert(none.at() == 0.);
  assert(none == 0.);
  assert(none + other == 5.);
  assert(none * other == 0.);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
