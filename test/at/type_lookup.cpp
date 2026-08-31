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
//! @brief White-box tests for the type-indexed `at<Type>()` lookup helpers.

#include "fcarouge/eigen.hpp"
#include "fcarouge/typed_linear_algebra.hpp"

#include <chrono>
#include <concepts>
#include <functional>
#include <ratio>
#include <tuple>

namespace fcarouge::test {
namespace {
//! @name Fixture types
//! @{

//! @brief A target type two element types below convert to.
struct common_target {};

//! @brief Converts to `common_target`, but not to `beta`.
struct alpha {
  // The implicit conversion is the fixture under test.
  // NOLINTNEXTLINE(google-explicit-constructor,hicpp-explicit-conversions)
  operator common_target() const;
};

//! @brief Converts to `common_target`, but not to `alpha`.
struct beta {
  // The implicit conversion is the fixture under test.
  // NOLINTNEXTLINE(google-explicit-constructor,hicpp-explicit-conversions)
  operator common_target() const;
};

struct base {};
struct derived : base {};

//! @brief Two ordinary durations: neither implicitly converts to the other
//! (each conversion would truncate), yet both share the common type
//! `std::chrono::duration<int, std::ratio<1, 15>>`.
using third = std::chrono::duration<int, std::ratio<1, 3>>;
using fifth = std::chrono::duration<int, std::ratio<1, 5>>;

//! @}

//! @test The element-type tuple is laid out in row-major order and matches the
//! rank-oblivious `element_at`.
[[maybe_unused]] const auto layout{[] -> int {
  using cell = typed_matrix<eigen::matrix<double, 1, 2>,
                            std::tuple<std::identity>, std::tuple<alpha, beta>>;
  using tuple = tla::tuple_typed_matrix<cell>;

  static_assert(std::tuple_size_v<tuple> == 2);
  static_assert(std::same_as<std::tuple_element_t<0, tuple>, alpha>);
  static_assert(std::same_as<std::tuple_element_t<1, tuple>, beta>);

  return 0;
}()};

//! @test An exact element type resolves to its own position.
[[maybe_unused]] const auto exact{[] -> int {
  using tuple = std::tuple<alpha, beta>;

  static_assert(tla::find_first_convertible_index<alpha, tuple>() == 0);
  static_assert(tla::find_first_convertible_index<beta, tuple>() == 1);

  return 0;
}()};

//! @test No convertible element yields the past-the-end sentinel and a zero
//! count.
[[maybe_unused]] const auto missing{[] -> int {
  using tuple = std::tuple<alpha, beta>;

  static_assert(tla::find_first_convertible_index<int, tuple>() == 2);
  static_assert(tla::find_first_convertible_index<base, tuple>() == 2);
  static_assert(tla::count_convertible_indexes<int, tuple>() == 0);
  static_assert(tla::count_convertible_indexes<base, tuple>() == 0);

  return 0;
}()};

//! @test The lookup matches an element that is convertible *to* the request,
//! not one the request is convertible *from*. Direction is asymmetric.
[[maybe_unused]] const auto direction{[] -> int {
  static_assert(
      tla::find_first_convertible_index<base, std::tuple<derived>>() == 0);
  static_assert(
      tla::find_first_convertible_index<derived, std::tuple<base>>() == 1);
  static_assert(tla::count_convertible_indexes<base, std::tuple<derived>>() ==
                1);
  static_assert(tla::count_convertible_indexes<derived, std::tuple<base>>() ==
                0);

  return 0;
}()};

//! @test The residual ambiguity the `distinct_typed_matrix` gate cannot see is
//! caught by the per-request count. `alpha` and `beta` are mutually
//! inconvertible and have no common type, so the matrix is still distinct, yet
//! both convert to `common_target` through a user-defined conversion. The count
//! is two, which is what `at<common_target>()` static-asserts against; the
//! first-match index is only consulted once the count is one.
[[maybe_unused]] const auto residual_ambiguity{[] -> int {
  using cell = typed_matrix<eigen::matrix<double, 1, 2>,
                            std::tuple<std::identity>, std::tuple<alpha, beta>>;
  static_assert(distinct_typed_matrix<cell>);
  static_assert(not tla::have_common_conversion_target<alpha, beta>);

  static_assert(tla::count_convertible_indexes<common_target,
                                               std::tuple<alpha, beta>>() == 2);
  static_assert(tla::count_convertible_indexes<common_target,
                                               std::tuple<beta, alpha>>() == 2);

  return 0;
}()};

//! @test The conversion-target gate catches a pair the plain
//! interconvertibility check would miss: types related only through a shared
//! common type. Genuinely unrelated types stay clear.
[[maybe_unused]] const auto conversion_target{[] -> int {
  static_assert(tla::have_common_conversion_target<int, double>);
  static_assert(tla::have_common_conversion_target<derived, base>);
  static_assert(tla::have_common_conversion_target<alpha, alpha>);

  // Neither duration implicitly converts to the other, yet they share a common
  // type: the plain interconvertibility check misses this, the common-type
  // disjunct catches it.
  static_assert(not tla::are_interconvertible<third, fifth>);
  static_assert(tla::have_common_conversion_target<third, fifth>);

  static_assert(not tla::have_common_conversion_target<alpha, beta>);
  static_assert(not tla::have_common_conversion_target<common_target, int>);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
