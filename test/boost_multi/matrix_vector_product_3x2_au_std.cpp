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

#include <au/units/meters.hh>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <mdspan>
#include <memory>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the matrix-vector product algorithm for a non-square,
//! three-by-two matrix and a two element column vector viewing memory owned by
//! Boost.Multi arrays, producing a three element column vector.
[[maybe_unused]] const auto test{[] -> int {
  using au::symbols::m;

  constexpr auto m3{au::cubed(m)};

  using length = au::QuantityD<au::Meters>;
  using volume = au::QuantityD<au::UnitPowerT<au::Meters, 3>>;
  using row_indexes = std::tuple<length, length, length>;
  using column_indexes = std::tuple<length, length>;
  using storage = boost::multi::inplace_array<representation, 2>;

  storage storage_a({3, 2}, 0.);
  storage storage_x({2, 1}, 0.);
  storage storage_y({3, 1}, 0.);

  storage_a.elements() = {1., 2., 3., 4., 5., 6.};
  storage_x.elements() = {1., -1.};

  std::mdspan span_a{std::to_address(storage_a.data_elements()),
                     std::extents<std::size_t, 3, 2>{}};
  std::mdspan span_x{std::to_address(storage_x.data_elements()),
                     std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_y{std::to_address(storage_y.data_elements()),
                     std::extents<std::size_t, 3, 1>{}};

  matrix<representation, row_indexes, column_indexes> a{span_a};
  column_vector<representation, length, length> x{span_x};
  column_vector<representation, volume, volume, volume> y{span_y};

  matrix_vector_product(a, x, y);

  assert((y.at<0>() == -1. * m3));
  assert((y.at<1>() == -1. * m3));
  assert((y.at<2>() == -1. * m3));
  assert(std::ranges::equal(storage_y.elements(), std::array{-1., -1., -1.}));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
