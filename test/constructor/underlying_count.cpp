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

#include <cstddef>
#include <tuple>
#include <utility>

namespace fcarouge::test {
namespace {
//! @brief An underlying matrix recording how many times its value was copied
//! and moved.
class counted_storage {
public:
  constexpr counted_storage() = default;

  constexpr counted_storage(const counted_storage &other)
      : copy_count{other.copy_count + 1}, move_count{other.move_count},
        value{other.value} {}

  constexpr counted_storage(counted_storage &&other) noexcept
      : copy_count{other.copy_count}, move_count{other.move_count + 1},
        value{other.value} {}

  constexpr auto operator=(const counted_storage &other)
      -> counted_storage & = default;

  constexpr auto operator=(counted_storage &&other) noexcept
      -> counted_storage & = default;

  constexpr ~counted_storage() = default;

  [[nodiscard]] constexpr auto operator[](std::size_t /*row*/,
                                          std::size_t /*column*/) -> double & {
    return value;
  }

  [[nodiscard]] constexpr auto operator[](std::size_t /*row*/,
                                          std::size_t /*column*/) const -> const
      double & {
    return value;
  }

  [[nodiscard]] constexpr auto copies() const -> int { return copy_count; }

  [[nodiscard]] constexpr auto moves() const -> int { return move_count; }

private:
  int copy_count{0};
  int move_count{0};
  double value{0.};
};

using matrix = typed_matrix<counted_storage, std::tuple<double>,
                            std::tuple<double, double>>;

//! @test Verifies the construction from an underlying matrix lvalue copies it
//! once, without a move.
static_assert([] -> bool {
  counted_storage storage;
  const matrix result{storage};

  return result.data().copies() == 1 && result.data().moves() == 0;
}());

//! @test Verifies the construction from a constant underlying matrix lvalue
//! copies it once, without a move.
static_assert([] -> bool {
  const counted_storage storage;
  const matrix result{storage};

  return result.data().copies() == 1 && result.data().moves() == 0;
}());

//! @test Verifies the construction from an underlying matrix xvalue moves it
//! once, without a copy.
static_assert([] -> bool {
  counted_storage storage;
  const matrix result{std::move(storage)};

  return result.data().copies() == 0 && result.data().moves() == 1;
}());

//! @test Verifies the construction from an underlying matrix prvalue moves it
//! once, without a copy.
static_assert([] -> bool {
  const matrix result{counted_storage{}};

  return result.data().copies() == 0 && result.data().moves() == 1;
}());
} // namespace
} // namespace fcarouge::test
