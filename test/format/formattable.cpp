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
#include <concepts>
#include <format>
#include <iterator>
#include <vector>

namespace fcarouge::test {
namespace {
// One type per format() overload: rank 0, column, row, mxn.
using rank0 = matrix<int, 1, 1>;
using column = matrix<int, 3, 1>;
using row = matrix<int, 1, 3>;
using rank2 = matrix<int, 3, 3>;

//! @test Verifies the formatter is semiregular, per the BasicFormatter
//! named requirement.
static_assert(std::semiregular<std::formatter<rank0, char>>);
static_assert(std::semiregular<std::formatter<column, char>>);
static_assert(std::semiregular<std::formatter<row, char>>);
static_assert(std::semiregular<std::formatter<rank2, char>>);
static_assert(std::semiregular<std::formatter<rank0, wchar_t>>);
static_assert(std::semiregular<std::formatter<column, wchar_t>>);
static_assert(std::semiregular<std::formatter<row, wchar_t>>);
static_assert(std::semiregular<std::formatter<rank2, wchar_t>>);

//! @test Verifies std::formattable, the BasicFormatter named requirement.
static_assert(std::formattable<rank0, char>);
static_assert(std::formattable<column, char>);
static_assert(std::formattable<row, char>);
static_assert(std::formattable<rank2, char>);
static_assert(std::formattable<rank0, wchar_t>);
static_assert(std::formattable<column, wchar_t>);
static_assert(std::formattable<row, wchar_t>);
static_assert(std::formattable<rank2, wchar_t>);

//! @test Verifies the format method is a FormatContext template.
[[maybe_unused]] const auto test{[] -> int {
  const rank2 m{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  assert(std::format("{}", m) == "[[1, 2, 3], [4, 5, 6], [7, 8, 9]]");

  const row r{1, 2, 3};
  std::vector<char> vector_out;
  std::format_to(std::back_inserter(vector_out), "{}", r);
  assert((vector_out ==
          std::vector<char>{'[', '1', ',', ' ', '2', ',', ' ', '3', ']'}));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
