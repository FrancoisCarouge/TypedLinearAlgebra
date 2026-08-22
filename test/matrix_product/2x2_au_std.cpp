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
#include <cstddef>
#include <mdspan>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the matrix product algorithm for a two-by-two matrix
//! shape.
[[maybe_unused]] const auto test{[] -> int {
  using length = au::QuantityD<au::Meters>;
  using indexes = std::tuple<length, length>;

  double storage_a[4]{};
  double storage_b[4]{};
  double storage_r[4]{};

  std::mdspan span_a{&storage_a[0], std::extents<std::size_t, 2, 2>{}};
  std::mdspan span_b{&storage_b[0], std::extents<std::size_t, 2, 2>{}};
  std::mdspan span_r{&storage_r[0], std::extents<std::size_t, 2, 2>{}};

  matrix<representation, indexes, indexes> a{span_a};
  matrix<representation, indexes, indexes> b{span_b};
  matrix<representation, indexes, indexes> r{span_r};

  a.at<0, 0>(au::squared(au::meters)(1.));
  a.at<0, 1>(au::squared(au::meters)(2.));
  a.at<1, 0>(au::squared(au::meters)(3.));
  a.at<1, 1>(au::squared(au::meters)(4.));

  b.at<0, 0>(au::squared(au::meters)(5.));
  b.at<0, 1>(au::squared(au::meters)(6.));
  b.at<1, 0>(au::squared(au::meters)(7.));
  b.at<1, 1>(au::squared(au::meters)(8.));

  matrix_product(a, b, r);

  assert((r.at<0, 0>() == au::squared(au::meters)(19.)));
  assert((r.at<0, 1>() == au::squared(au::meters)(22.)));
  assert((r.at<1, 0>() == au::squared(au::meters)(43.)));
  assert((r.at<1, 1>() == au::squared(au::meters)(50.)));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
