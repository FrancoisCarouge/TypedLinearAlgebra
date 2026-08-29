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
//! @brief Time safe linear algebra with std::chrono and std::linalg.
//!
//! @details Demonstrate the std::linalg-style free functions over a vector
//! whose elements are `std::chrono` durations of distinct periods. This
//! library composes std::mdspan as the matrix' linear algebra backend with
//! indexes typed as `std::chrono` durations. Unlike some third party unit
//! libraries, `std::chrono` durations assign through an rvalue, so the
//! `add()` free function applies here.

#include "fcarouge/linalg.hpp"

#include <cassert>
#include <chrono>
#include <cstddef>
#include <format>
#include <linalg>
#include <mdspan>
#include <print>
#include <ratio>
#include <vector>

namespace fcarouge::sample {
namespace {
using representation = double;
using seconds = std::chrono::duration<representation>;
using minutes = std::chrono::duration<representation, std::ratio<60>>;
using hours = std::chrono::duration<representation, std::ratio<3600>>;

template <typename... Types>
using column_vector = column_vector<representation, Types...>;

template <std::size_t Rows>
using column_extents = std::extents<std::size_t, Rows, 1>;

//! @brief Time typed linear algebra samples with std::mdspan and std::linalg.
[[maybe_unused]] const auto sample{[] -> int {
  using durations = column_vector<seconds, minutes, hours>;

  std::vector v0(std::size_t{3}, representation{});
  std::mdspan s0{v0.data(), column_extents<3>{}};
  durations x0{s0};

  // Element assignment keeps each period.
  x0.at<0>(seconds{30.});
  x0.at<1>(minutes{2.});
  x0.at<2>(hours{1.});

  std::println("x0 = {}", x0);
  assert(std::format("{}", x0) == "[[30s], [2min], [1h]]");
  assert(x0.at<1>() == minutes{2.});

  // In-place multiplication with a scalar factor.
  scale(3., x0);
  assert(x0.at<0>() == seconds{90.});
  assert(x0.at<2>() == hours{3.});

  // Out-of-place addition of two vectors of the same types.
  std::vector v1(std::size_t{3}, representation{});
  std::mdspan s1{v1.data(), column_extents<3>{}};
  durations x1{s1};
  add(x0, x0, x1);
  assert(x1.at<1>() == minutes{12.});

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
