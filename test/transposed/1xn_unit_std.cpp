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
#include <cstddef>
#include <mdspan>
#include <tuple>

namespace fcarouge::test {
namespace {
//! @test Verifies the transposed algorithm.
[[maybe_unused]] const auto test{[] {
  double storage[]{9., 10.};
  std::mdspan span{&storage[0], std::extents<std::size_t, 2, 1>{}};
  column_vector<double, decltype(1. * s), decltype(1. * s2)> n{span};

  assert(n.at<0>() == 9. * s);
  assert(n.at<1>() == 10. * s2);
  static_assert(std::same_as<decltype(n)::row_indexes,
                             std::tuple<decltype(1. * s), decltype(1. * s2)>>);

  n.at<0>() = 42. * s;
  n.at<1>() = 43. * s2;

  assert(n.at<0>() == 42. * s);
  assert(n.at<1>() == 43. * s2);

  auto nᵀ{transposed(n)};

  assert(nᵀ.at<0>() == 42. * s);
  assert(nᵀ.at<1>() == 43. * s2);

  static_assert(std::same_as<decltype(nᵀ)::column_indexes,
                             std::tuple<decltype(1. * s), decltype(1. * s2)>>);

  n.at<0>() = 24. * s;
  n.at<1>() = 25. * s2;

  assert(nᵀ.at<0>() == 24. * s);
  assert(nᵀ.at<1>() == 25. * s2);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
