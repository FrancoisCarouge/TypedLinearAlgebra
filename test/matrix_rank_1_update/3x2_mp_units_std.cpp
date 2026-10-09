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
#include <cstddef>
#include <mdspan>
#include <tuple>

namespace fcarouge::test {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::m2;

namespace {
//! @test Verifies the overwriting matrix rank-1 update algorithm, the outer
//! product of a three element and a two element column vector, into a
//! non-square, three-by-two matrix.
[[maybe_unused]] const auto test{[] -> int {
  using length = quantity<mp_units::isq::length[m]>;
  using row_indexes = std::tuple<length, length, length>;
  using column_indexes = std::tuple<length, length>;

  representation storage_lhs[3]{};
  representation storage_rhs[2]{};
  // Stale values the overwriting update must replace, not accumulate onto.
  representation storage_result[6]{9., 9., 9., 9., 9., 9.};

  std::mdspan span_lhs{&storage_lhs[0], std::extents<std::size_t, 3, 1>{}};
  std::mdspan span_rhs{&storage_rhs[0], std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_result{&storage_result[0],
                          std::extents<std::size_t, 3, 2>{}};

  column_vector<representation, length, length, length> lhs{span_lhs};
  column_vector<representation, length, length> rhs{span_rhs};
  matrix<representation, row_indexes, column_indexes> result{span_result};

  lhs.at<0>(1. * m);
  lhs.at<1>(2. * m);
  lhs.at<2>(3. * m);

  rhs.at<0>(4. * m);
  rhs.at<1>(-5. * m);

  matrix_rank_1_update(lhs, rhs, result);

  assert((result.at<0, 0>() == 4. * m2));
  assert((result.at<0, 1>() == -5. * m2));
  assert((result.at<1, 0>() == 8. * m2));
  assert((result.at<1, 1>() == -10. * m2));
  assert((result.at<2, 0>() == 12. * m2));
  assert((result.at<2, 1>() == -15. * m2));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
