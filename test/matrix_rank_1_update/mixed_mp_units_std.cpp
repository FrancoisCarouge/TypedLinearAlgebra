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
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;

namespace {
//! @test Verifies the matrix rank-1 update algorithm types each element of
//! the result as the product of the corresponding position and time elements
//! of the vectors.
[[maybe_unused]] const auto test{[] -> int {
  using length = quantity<mp_units::isq::length[m]>;
  using time = quantity<mp_units::isq::time[s]>;
  using indexes = std::tuple<length, time>;

  representation storage_lhs[2]{};
  representation storage_rhs[2]{};
  representation storage_result[4]{};

  std::mdspan span_lhs{&storage_lhs[0], std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_rhs{&storage_rhs[0], std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_result{&storage_result[0],
                          std::extents<std::size_t, 2, 2>{}};

  column_vector<representation, length, time> lhs{span_lhs};
  column_vector<representation, length, time> rhs{span_rhs};
  matrix<representation, indexes, indexes> result{span_result};

  lhs.at<0>(2. * m);
  lhs.at<1>(3. * s);

  rhs.at<0>(5. * m);
  rhs.at<1>(-7. * s);

  matrix_rank_1_update(lhs, rhs, result);

  assert((result.at<0, 0>() == 10. * m2));
  assert((result.at<0, 1>() == -14. * m * s));
  assert((result.at<1, 0>() == 15. * s * m));
  assert((result.at<1, 1>() == -21. * s2));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
