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

#include <concepts>
#include <cstddef>
#include <functional>
#include <mdspan>
#include <tuple>
#include <type_traits>

namespace fcarouge::test {
namespace {
//! @test Verifies a typed matrix over index tuples whose own definition
//! evaluates their tuple size, as libc++'s `std::tuple` does, compiles.
//!
//! @details The typed matrix is the first to name the index tuple type, so
//! selecting its `std::tuple_size` specialization is what completes it.
//! Completing it evaluates that same tuple size again, which is only a
//! constant expression if the selection did not need it complete already.
template <typename... Types> struct self_sized {
  using size =
      std::integral_constant<std::size_t, std::tuple_size_v<self_sized>>;
};
} // namespace
} // namespace fcarouge::test

template <typename... Types>
struct std::tuple_size<fcarouge::test::self_sized<Types...>>
    : std::integral_constant<std::size_t, sizeof...(Types)> {};

template <std::size_t Index, typename... Types>
struct std::tuple_element<Index, fcarouge::test::self_sized<Types...>>
    : std::tuple_element<Index, std::tuple<Types...>> {};

namespace fcarouge::test {
namespace {
using column =
    typed_matrix<std::mdspan<double, std::extents<std::size_t, 2, 1>>,
                 self_sized<double, double>, std::tuple<std::identity>>;

static_assert(column::rows == 2);
static_assert(column::columns == 1);
static_assert(rank_typed_matrix<column, 1>);
static_assert(std::tuple_size_v<column> == 2);
static_assert(std::same_as<std::tuple_element_t<1, column>, double>);
static_assert(not same_as_typed_matrix<self_sized<int, int, int>>);
static_assert(std::tuple_size_v<self_sized<int, int, int>> == 3);
} // namespace
} // namespace fcarouge::test
