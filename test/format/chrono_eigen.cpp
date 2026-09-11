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

#include "fcarouge/linalg.hpp"

#include <cassert>
#include <chrono>
#include <concepts>
#include <format>

namespace fcarouge::test {
using representation = double;
using seconds = std::chrono::duration<representation>;

namespace {
using singleton = row_vector<representation, seconds>;
using row = row_vector<representation, seconds, seconds, seconds>;
using column = column_vector<representation, seconds, seconds, seconds>;

//! @test Verifies the formatter properties for the std::chrono and Eigen
//! composition.
static_assert(std::semiregular<std::formatter<singleton, char>>);
static_assert(std::semiregular<std::formatter<row, char>>);
static_assert(std::semiregular<std::formatter<column, char>>);
static_assert(std::formattable<singleton, char>);
static_assert(std::formattable<row, char>);
static_assert(std::formattable<column, char>);

//! @test Verifies the format algorithm for std::chrono duration typed matrices
//! with the Eigen backend: the rank-0 singleton, the row vector, and the column
//! vector overloads.
[[maybe_unused]] const auto test{[] -> int {
  const singleton s{seconds{9.}};
  assert(std::format("{}", s) == "9s");

  const row r{seconds{1.}, seconds{2.}, seconds{3.}};
  assert(std::format("{}", r) == "[1s, 2s, 3s]");

  const column c{seconds{1.}, seconds{2.}, seconds{3.}};
  assert(std::format("{}", c) == "[[1s], [2s], [3s]]");

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
