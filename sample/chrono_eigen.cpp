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

//! @file
//! @brief Time safe linear algebra with std::chrono and Eigen.
//!
//! @details Demonstrate a variety of linear algebra operations over a vector
//! whose elements are `std::chrono` durations of distinct periods. This
//! library composes Eigen as the matrix' linear algebra backend with indexes
//! typed as `std::chrono` durations. This sample explicitly uses double
//! precision floating point representation.

#include "fcarouge/linalg.hpp"

#include <cassert>
#include <chrono>
#include <format>
#include <print>
#include <ratio>

namespace fcarouge::sample {
namespace {
using representation = double;
using seconds = std::chrono::duration<representation>;
using minutes = std::chrono::duration<representation, std::ratio<60>>;
using hours = std::chrono::duration<representation, std::ratio<3600>>;

template <typename... Types>
using column_vector = column_vector<representation, Types...>;

//! @brief Time typed linear algebra samples.
[[maybe_unused]] const auto sample{[] -> int {
  // A heterogeneous vector: each element keeps its own period.
  using durations = column_vector<seconds, minutes, hours>;
  durations x0{seconds{30.}, minutes{2.}, hours{1.}};

  // Printable.
  std::println("x0 = {}", x0);
  assert(std::format("{}", x0) == "[[30s], [2min], [1h]]");

  // Element access keeps the period.
  assert(x0.at<1>() == minutes{2.});
  assert(std::format("{}", x0.at<1>()) == "2min");

  // Multiplication with a scalar factor.
  const durations x1{x0 * 3.};
  assert(x1.at<0>() == seconds{90.});
  assert(x1.at<1>() == minutes{6.});
  assert(x1.at<2>() == hours{3.});

  // Division with a scalar divisor.
  const durations x2{x1 / 2.};
  assert(x2.at<0>() == seconds{45.});

  // Substraction then addition of vectors of the same types.
  const durations x3{x2 - x0};
  assert(x3.at<2>() == hours{0.5});
  const durations x4{x3 + x3};
  assert(x4.at<1>() == minutes{2.});

  // Euclidean L2 norm of a uniform duration vector.
  const column_vector<seconds, seconds, seconds> v{seconds{2.}, seconds{3.},
                                                   seconds{6.}};
  assert(std::format("{}", magnitude(v)) == "7s");

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
