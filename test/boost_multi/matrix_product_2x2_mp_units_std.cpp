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

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <mdspan>
#include <memory>
#include <tuple>

namespace fcarouge::test {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::m2;

namespace {
//! @test Verifies the matrix product algorithm for two-by-two matrices viewing
//! memory owned by Boost.Multi arrays, read and written in row-major order.
[[maybe_unused]] const auto test{[] -> int {
  using length = quantity<mp_units::isq::length[m]>;
  using indexes = std::tuple<length, length>;
  using storage = boost::multi::inplace_array<representation, 2>;

  storage storage_a({2, 2}, 0.);
  storage storage_b({2, 2}, 0.);
  storage storage_r({2, 2}, 0.);

  storage_a.elements() = {1., 2., 3., 4.};
  storage_b.elements() = {5., 6., 7., 8.};

  std::mdspan span_a{std::to_address(storage_a.data_elements()),
                     std::extents<std::size_t, 2, 2>{}};
  std::mdspan span_b{std::to_address(storage_b.data_elements()),
                     std::extents<std::size_t, 2, 2>{}};
  std::mdspan span_r{std::to_address(storage_r.data_elements()),
                     std::extents<std::size_t, 2, 2>{}};

  matrix<representation, indexes, indexes> a{span_a};
  matrix<representation, indexes, indexes> b{span_b};
  matrix<representation, indexes, indexes> r{span_r};

  matrix_product(a, b, r);

  assert((a.at<0, 1>() == 2. * m2));
  assert((a.at<1, 0>() == 3. * m2));
  assert((r.at<0, 0>() == 19. * m2));
  assert((r.at<0, 1>() == 22. * m2));
  assert((r.at<1, 0>() == 43. * m2));
  assert((r.at<1, 1>() == 50. * m2));
  assert(
      std::ranges::equal(storage_r.elements(), std::array{19., 22., 43., 50.}));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
