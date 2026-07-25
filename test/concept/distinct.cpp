/* Typed Linear Algebra
Version 0.2.0
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

#include <print>

namespace fcarouge::test {
namespace {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using position = quantity<mp_units::isq::length[m]>;
using velocity = quantity<mp_units::isq::velocity[m / s]>;
using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;

//! @test Verifies the distinct typed matrix concept.
[[maybe_unused]] const auto test{[] {
  ///////////////////////
  // distinct_typed_matrix
  {
    //! @todo Permit non-tuple declaration.
    using s = matrix<double, std::tuple<position>, std::tuple<std::identity>>;
    static_assert(distinct_typed_matrix<s>);
    static_assert(uniform_typed_matrix<s>); // that's funny
  }
  {
    using c = matrix<double, std::tuple<position, velocity>,
                     std::tuple<std::identity>>;
    static_assert(distinct_typed_matrix<c>);
    static_assert(not uniform_typed_matrix<c>);
  }
  {
    using c = matrix<double, std::tuple<velocity, velocity>,
                     std::tuple<std::identity>>;
    static_assert(not distinct_typed_matrix<c>);
    static_assert(uniform_typed_matrix<c>);
  }
  {
    using r = matrix<double, std::tuple<std::identity>,
                     std::tuple<position, velocity>>;
    static_assert(distinct_typed_matrix<r>);
    static_assert(not uniform_typed_matrix<r>);
  }
  {
    using r = matrix<double, std::tuple<std::identity>,
                     std::tuple<position, position>>;
    static_assert(not distinct_typed_matrix<r>);
    static_assert(uniform_typed_matrix<r>);
  }
  {
    using m = matrix<double, std::tuple<double, position>,
                     std::tuple<velocity, acceleration>>;
    static_assert(distinct_typed_matrix<m>);
    static_assert(not uniform_typed_matrix<m>);
  }
  {
    using m = matrix<double, std::tuple<position, velocity>,
                     std::tuple<velocity, acceleration>>;
    static_assert(not distinct_typed_matrix<m>);
    static_assert(not uniform_typed_matrix<m>);
  }
  {
    using m = matrix<double, std::tuple<position, position>,
                     std::tuple<position, position>>;
    static_assert(not distinct_typed_matrix<m>);
    static_assert(uniform_typed_matrix<m>);
  }
  return 0;
}()};
} // namespace
} // namespace fcarouge::test
