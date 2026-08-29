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
#include <tuple>

namespace fcarouge::test {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::m2;

namespace {
//! @test Verifies the transposed algorithm for a rectangular matrix shape
//! with non-trivial types.
[[maybe_unused]] const auto test{[] -> int {
  using length = quantity<mp_units::isq::length[m]>;
  using row_indexes = std::tuple<length, length>;
  using column_indexes = std::tuple<length, length, length>;

  matrix<representation, row_indexes, column_indexes> a;

  a.at<0, 0>(1. * m2);
  a.at<0, 1>(2. * m2);
  a.at<0, 2>(3. * m2);
  a.at<1, 0>(4. * m2);
  a.at<1, 1>(5. * m2);
  a.at<1, 2>(6. * m2);

  matrix<representation, column_indexes, row_indexes> aᵀ{transposed(a)};

  static_assert(std::same_as<decltype(aᵀ)::row_indexes, column_indexes>);
  static_assert(std::same_as<decltype(aᵀ)::column_indexes, row_indexes>);

  assert((aᵀ.at<0, 0>() == 1. * m2));
  assert((aᵀ.at<1, 0>() == 2. * m2));
  assert((aᵀ.at<2, 0>() == 3. * m2));
  assert((aᵀ.at<0, 1>() == 4. * m2));
  assert((aᵀ.at<1, 1>() == 5. * m2));
  assert((aᵀ.at<2, 1>() == 6. * m2));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
