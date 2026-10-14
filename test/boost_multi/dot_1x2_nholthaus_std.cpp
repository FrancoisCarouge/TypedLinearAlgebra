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

#include <boost/multi/array.hpp>

#include <units/area.h>
#include <units/length.h>

#include <cassert>
#include <cstddef>
#include <mdspan>
#include <memory>

namespace fcarouge::test {
using units::m;
using units::m2;
using representation = double;

namespace {
//! @test Verifies the dot product of two row vectors of quantities viewing
//! memory owned by Boost.Multi arrays.
[[maybe_unused]] const auto test{[] -> int {
  using length = units::length::meters<representation>;
  using storage = boost::multi::inplace_array<representation, 2>;

  storage storage_a({1, 2}, 0.);
  storage storage_b({1, 2}, 0.);

  storage_a.elements() = {2., -3.};
  storage_b.elements() = {4., 5.};

  std::mdspan span_a{std::to_address(storage_a.data_elements()),
                     std::extents<std::size_t, 1, 2>{}};
  std::mdspan span_b{std::to_address(storage_b.data_elements()),
                     std::extents<std::size_t, 1, 2>{}};

  row_vector<representation, length, length> a{span_a};
  row_vector<representation, length, length> b{span_b};

  assert(a.at<1>() == -3. * m);
  assert(dot(a, b) == -7. * m2);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
