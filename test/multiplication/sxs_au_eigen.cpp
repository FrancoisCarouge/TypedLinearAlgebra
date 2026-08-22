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
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the square by square matrix multiplication operator.
[[maybe_unused]] const auto test{[] -> int {
  using length = au::QuantityD<au::Meters>;
  using area = au::QuantityD<au::UnitPowerT<au::Meters, 2>>;
  using indexes = std::tuple<length, length>;
  using result_indexes = std::tuple<area, area>;

  const matrix<representation, indexes, indexes> a{
      {au::squared(au::meters)(1.), au::squared(au::meters)(2.)},
      {au::squared(au::meters)(3.), au::squared(au::meters)(4.)}};
  const matrix<representation, indexes, indexes> b{
      {au::squared(au::meters)(5.), au::squared(au::meters)(6.)},
      {au::squared(au::meters)(7.), au::squared(au::meters)(8.)}};
  const matrix<representation, result_indexes, result_indexes> r{a * b};

  assert((r.at<0, 0>() == au::pow<4>(au::meters)(19.)));
  assert((r.at<0, 1>() == au::pow<4>(au::meters)(22.)));
  assert((r.at<1, 0>() == au::pow<4>(au::meters)(43.)));
  assert((r.at<1, 1>() == au::pow<4>(au::meters)(50.)));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
