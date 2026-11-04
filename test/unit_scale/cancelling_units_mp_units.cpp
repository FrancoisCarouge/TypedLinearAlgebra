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

#include <cassert>
#include <utility>

namespace fcarouge::test {
using representation = double;

namespace {
//! @test Verifies a sum of terms in different units of one scale still
//! computes: each term cancels its own unit, so their stored numbers can be
//! summed.
[[maybe_unused]] const auto test{[] -> int {
  using mp_units::si::unit_symbols::km;
  using mp_units::si::unit_symbols::m;

  using small_unit =
      mp_units::quantity<mp_units::isq::length[m], representation>;
  using large_unit =
      mp_units::quantity<mp_units::isq::length[km], representation>;
  using inverse_small = decltype(representation{} / std::declval<small_unit>());
  using inverse_large = decltype(representation{} / std::declval<large_unit>());

  const large_unit first{3. * km};
  const small_unit second{4. * m};
  const large_unit one_large{1. * km};
  const small_unit one_small{1. * m};

  const row_vector<representation, inverse_large, inverse_small> weights{
      1. / one_large, 1. / one_small};
  const column_vector<representation, large_unit, small_unit> values{first,
                                                                     second};

  const auto result{weights * values};

  assert(result.at<>() == first / one_large + second / one_small);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
