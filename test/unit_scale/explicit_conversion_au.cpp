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
of the public at large_unit and to the detriment of our heirs and
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

#include <au/prefix.hh>
#include <au/units/meters.hh>

#include <cassert>
#include <functional>
#include <tuple>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the operations compute the expected values on Au quantities
//! in different units, meters and kilometers, once explicitly converted to one
//! unit.
[[maybe_unused]] const auto test{[] -> int {
  using small_unit = au::QuantityD<au::Meters>;
  using large_unit = au::QuantityD<au::Kilo<au::Meters>>;
  using gain =
      matrix<representation, std::tuple<representation, representation>,
             std::tuple<std::identity, std::identity>>;
  using small_column = column_vector<representation, small_unit, small_unit>;
  using large_column = column_vector<representation, large_unit, large_unit>;
  using scalars = std::tuple<representation, representation>;
  using small_matrix =
      matrix<representation, std::tuple<small_unit, small_unit>, scalars>;
  using large_matrix =
      matrix<representation, std::tuple<large_unit, large_unit>, scalars>;
  using mixed_column = column_vector<representation, small_unit, large_unit>;
  using mixed_row = row_vector<representation, small_unit, large_unit>;
  using small_row = row_vector<representation, small_unit, small_unit>;
  using scalar_row = row_vector<representation, representation, representation>;
  using scalar_column =
      column_vector<representation, representation, representation>;

  const small_unit first{au::meters(500.)};
  const large_unit second{au::kilo(au::meters)(3.)};
  const small_column small_operand{first, first};
  const large_column large_operand{second, second};
  const mixed_column mixed{first, second};
  const mixed_row mixed_transposed{first, second};
  const gain identity{{1., 0.}, {0., 1.}};
  const gain ones{{1., 1.}, {1., 1.}};

  // Addition and subtraction.
  const small_column sum{small_operand + small_column{large_operand}};
  const small_column difference{small_operand - small_column{large_operand}};

  assert(sum.at<0>() == first + second);
  assert(sum.at<1>() == first + second);
  assert(difference.at<0>() == first - second);
  assert(difference.at<1>() == first - second);

  // Whichever unit the result is expressed in.
  const large_column large_sum{large_column{small_operand} + large_operand};

  assert(large_sum.at<0>() == first + second);

  // Matrices.
  small_matrix matrix_lhs;
  large_matrix matrix_rhs;

  matrix_lhs.at<0, 0>(first);
  matrix_rhs.at<1, 1>(second);

  const small_matrix matrix_sum{matrix_lhs + small_matrix{matrix_rhs}};

  assert((matrix_sum.at<0, 0>() == first));
  assert((matrix_sum.at<1, 1>() == second));

  // Products.
  const small_column kept{identity * small_column{mixed}};
  const small_column summed{ones * small_column{mixed}};
  const small_unit inner{scalar_row{1., 1.} * small_column{mixed}};
  const small_row kept_row{small_row{mixed_transposed} * identity};

  assert(kept.at<0>() == first);
  assert(kept.at<1>() == second);
  assert(summed.at<0>() == first + second);
  assert(summed.at<1>() == first + second);
  assert(inner == first + second);
  assert(kept_row.at<1>() == second);

  // The dot product.
  const small_unit dotted{dot(small_column{mixed}, scalar_column{1., 1.})};

  assert(dotted == first + second);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
