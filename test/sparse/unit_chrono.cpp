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
#include <chrono>
#include <concepts>
#include <functional>
#include <ratio>
#include <tuple>
#include <type_traits>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies the sparse storage with std::chrono durations: the
//! per-element typing of the indexes, whatever the pattern, a structural zero
//! included.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  using seconds = std::chrono::duration<representation>;
  using minutes = std::chrono::duration<representation, std::ratio<60>>;

  // Heterogeneous periods with a structural zero of its own period.
  using state =
      sparse_column_vector<representation, sparse::make_pattern<2, 1>({{0, 0}}),
                           seconds, minutes>;

  state value{sparse::make_matrix<state>(seconds{2.}, sparse::structural_zero)};

  static_assert(std::same_as<decltype(value.at<1>()), const minutes>);
  assert(value.at<0>() == seconds{2.});
  assert(value.at<1>() == minutes{0.});

  value.at<0>(seconds{3.});

  assert(value.at<0>() == seconds{3.});

  // A conversion of the given period on construction.
  const state converted{
      sparse::make_matrix<state>(minutes{1.}, sparse::structural_zero)};

  assert(converted.at<0>() == seconds{60.});

  // A sum, a difference, a negation, a scaling: each with its element period.
  constexpr auto second_only{sparse::make_pattern<2, 1>({{1, 0}})};

  using second_state =
      sparse_column_vector<representation, second_only, seconds, minutes>;

  const second_state other{
      sparse::make_matrix<second_state>(sparse::structural_zero, minutes{4.})};

  const auto sum{value + other};

  static_assert(decltype(sum)::matrix::structure.count == 2);
  assert(sum.at<0>() == seconds{3.});
  assert(sum.at<1>() == minutes{4.});
  assert((value - other).at<1>() == minutes{-4.});
  assert((-other).at<1>() == minutes{-4.});
  assert((other * 2.).at<1>() == minutes{8.});
  assert((2. * other).at<1>() == minutes{8.});

  // A transposition, the row vector keeps the periods.
  const auto row{transposed(other)};

  assert(row.at<0>() == seconds{0.});
  assert(row.at<1>() == minutes{4.});

  // A dimensionless, diagonal gain applied to durations.
  using identity_indexes = std::tuple<std::identity, std::identity>;

  const sparse_matrix<representation, sparse::diagonal<2>,
                      std::tuple<representation, representation>,
                      identity_indexes>
      gain{sparse::make_matrix<sparse_matrix<
          representation, sparse::diagonal<2>,
          std::tuple<representation, representation>, identity_indexes>>(
          2., sparse::structural_zero, sparse::structural_zero, 3.)};
  const sparse_column_vector<representation, sparse::dense<2, 1>, seconds,
                             seconds>
      durations{seconds{1.}, seconds{4.}};
  const auto image{gain * durations};

  static_assert(
      std::same_as<std::remove_cvref_t<decltype(image.at<1>())>, seconds>);
  assert(image.at<0>() == seconds{2.});
  assert(image.at<1>() == seconds{12.});

  // A uniform matrix of durations, its norms in the element type.
  using uniform_indexes = std::tuple<seconds, seconds>;

  const sparse_matrix<representation, sparse::diagonal<2>, uniform_indexes,
                      identity_indexes>
      uniform{
          sparse::make_matrix<sparse_matrix<representation, sparse::diagonal<2>,
                                            uniform_indexes, identity_indexes>>(
              seconds{-3.}, sparse::structural_zero, sparse::structural_zero,
              seconds{4.})};

  assert(matrix_one_norm(uniform) == seconds{4.});
  assert(matrix_inf_norm(uniform) == seconds{4.});
  assert(matrix_frob_norm(uniform) == seconds{5.});

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
