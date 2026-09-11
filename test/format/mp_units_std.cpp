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

#include <cassert>
#include <concepts>
#include <cstddef>
#include <format>
#include <mdspan>
#include <tuple>

namespace fcarouge::test {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::m2;

using length = quantity<mp_units::isq::length[m]>;

namespace {
using singleton = row_vector<representation, length>;
using row = row_vector<representation, length, length, length>;
using column = column_vector<representation, length, length, length>;
using rectangle = matrix<representation, std::tuple<length, length>,
                         std::tuple<length, length>>;

//! @test Verifies the formatter properties for the mp-units and mdspan
//! composition.
static_assert(std::semiregular<std::formatter<singleton, char>>);
static_assert(std::semiregular<std::formatter<row, char>>);
static_assert(std::semiregular<std::formatter<column, char>>);
static_assert(std::semiregular<std::formatter<rectangle, char>>);
static_assert(std::formattable<singleton, char>);
static_assert(std::formattable<row, char>);
static_assert(std::formattable<column, char>);
static_assert(std::formattable<rectangle, char>);

//! @test Verifies the format algorithm for mp-units quantity typed matrices
//! with the mdspan backend: the rank-0 singleton, the row vector, the column
//! vector, and the general rectangular overloads.
[[maybe_unused]] const auto test{[] -> int {
  representation singleton_storage{};
  std::mdspan singleton_span{&singleton_storage,
                             std::extents<std::size_t, 1, 1>{}};
  singleton s{singleton_span};
  s.at(9. * m);
  assert(std::format("{}", s) == "9 m");

  representation row_storage[]{{}, {}, {}};
  std::mdspan row_span{&row_storage[0], std::extents<std::size_t, 1, 3>{}};
  row r{row_span};
  r.at<0>(1. * m);
  r.at<1>(2. * m);
  r.at<2>(3. * m);
  assert(std::format("{}", r) == "[1 m, 2 m, 3 m]");

  representation column_storage[]{{}, {}, {}};
  std::mdspan column_span{&column_storage[0],
                          std::extents<std::size_t, 3, 1>{}};
  column c{column_span};
  c.at<0>(1. * m);
  c.at<1>(2. * m);
  c.at<2>(3. * m);
  assert(std::format("{}", c) == "[[1 m], [2 m], [3 m]]");

  representation rectangle_storage[]{{}, {}, {}, {}};
  std::mdspan rectangle_span{&rectangle_storage[0],
                             std::extents<std::size_t, 2, 2>{}};
  rectangle e{rectangle_span};
  e.at<0, 0>(1. * m2);
  e.at<0, 1>(2. * m2);
  e.at<1, 0>(3. * m2);
  e.at<1, 1>(4. * m2);
  assert(std::format("{}", e) == "[[1 m², 2 m²], [3 m², 4 m²]]");

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
