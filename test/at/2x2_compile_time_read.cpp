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

#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>

namespace fcarouge::test {
namespace {
//! @brief An underlying matrix counting the reads of its elements at a run-time
//! position, as a storage whose run-time position needs a lookup.
class positioned_storage {
public:
  [[nodiscard]] constexpr auto operator[](std::size_t row, std::size_t column)
      -> double & {
    ++lookup_count;

    return values.at((row * 2) + column);
  }

  [[nodiscard]] constexpr auto operator[](std::size_t row,
                                          std::size_t column) const -> double {
    ++lookup_count;

    return values.at((row * 2) + column);
  }

  template <std::size_t Row, std::size_t Column>
  [[nodiscard]] constexpr auto
  operator[](std::integral_constant<std::size_t, Row> /*row*/,
             std::integral_constant<std::size_t, Column> /*column*/)
      -> double & {
    return std::get<(Row * 2) + Column>(values);
  }

  template <std::size_t Row, std::size_t Column>
  [[nodiscard]] constexpr auto
  operator[](std::integral_constant<std::size_t, Row> /*row*/,
             std::integral_constant<std::size_t, Column> /*column*/) const
      -> double {
    return std::get<(Row * 2) + Column>(values);
  }

  [[nodiscard]] constexpr auto lookups() const -> int { return lookup_count; }

private:
  std::array<double, 4> values{};
  mutable int lookup_count{0};
};

using matrix = typed_matrix<positioned_storage, std::tuple<double, double>,
                            std::tuple<double, double>>;

//! @test Verifies an element written and read at a compile-time position
//! offers the storage its compile-time position, without a run-time lookup,
//! through a modifiable and a constant typed matrix.
static_assert([] -> bool {
  matrix value;

  value.at<0, 1>(2.);
  value.at<1, 0>(3.);

  const matrix &constant{value};

  return value.at<0, 1>() == 2. && constant.at<1, 0>() == 3. &&
         constant.at<1, 1>() == 0. && value.data().lookups() == 0;
}());
} // namespace
} // namespace fcarouge::test
