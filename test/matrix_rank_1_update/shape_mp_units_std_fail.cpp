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
#include <mdspan>
#include <tuple>

namespace fcarouge::test {
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

using mp_units::si::unit_symbols::m;

namespace {
//! @test Verifies the updating matrix rank-1 update algorithm rejects an `E`
//! matrix whose shape differs from the result matrix.
[[maybe_unused]] const auto test{[] -> int {
  using length = quantity<mp_units::isq::length[m]>;
  using indexes = std::tuple<length, length>;

  representation storage_lhs[2]{};
  representation storage_rhs[2]{};
  representation storage_addend[6]{};
  representation storage_result[4]{};

  std::mdspan span_lhs{&storage_lhs[0], std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_rhs{&storage_rhs[0], std::extents<std::size_t, 2, 1>{}};
  std::mdspan span_addend{&storage_addend[0],
                          std::extents<std::size_t, 2, 3>{}};
  std::mdspan span_result{&storage_result[0],
                          std::extents<std::size_t, 2, 2>{}};

  column_vector<representation, length, length> lhs{span_lhs};
  column_vector<representation, length, length> rhs{span_rhs};

  // Intended:
  // matrix<representation, indexes, indexes> addend{...};

  matrix<representation, indexes, std::tuple<length, length, length>> addend{
      span_addend};

  matrix<representation, indexes, indexes> result{span_result};

  matrix_rank_1_update(lhs, rhs, addend, result);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
