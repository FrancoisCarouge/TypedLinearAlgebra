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

namespace fcarouge::test {
//! @brief A value counting the arithmetic operations it results from, without
//! any global state, so the count is a compile-time verification.
struct counted {
  double value{0.};
  std::size_t multiplications{0};
  std::size_t additions{0};

  [[nodiscard]] friend constexpr auto operator*(const counted &lhs,
                                                const counted &rhs) -> counted {
    return {.value = lhs.value * rhs.value,
            .multiplications = lhs.multiplications + rhs.multiplications + 1,
            .additions = lhs.additions + rhs.additions};
  }

  [[nodiscard]] friend constexpr auto operator+(const counted &lhs,
                                                const counted &rhs) -> counted {
    return {.value = lhs.value + rhs.value,
            .multiplications = lhs.multiplications + rhs.multiplications,
            .additions = lhs.additions + rhs.additions + 1};
  }

  friend constexpr auto operator+=(counted &lhs, const counted &rhs)
      -> counted & {
    lhs = lhs + rhs;
    return lhs;
  }
};

namespace {
template <auto Pattern> using storage = sparse::matrix<counted, Pattern>;

//! @brief The sums of the counts of a matrix's stored elements.
struct total {
  template <auto Pattern>
  [[nodiscard]] static consteval auto
  multiplications(const storage<Pattern> &value) -> std::size_t {
    std::size_t result{0};

    for (const counted &element : value.values) {
      result += element.multiplications;
    }

    return result;
  }

  template <auto Pattern>
  [[nodiscard]] static consteval auto additions(const storage<Pattern> &value)
      -> std::size_t {
    std::size_t result{0};

    for (const counted &element : value.values) {
      result += element.additions;
    }

    return result;
  }
};

using sparse::dense;
using sparse::diagonal;
using sparse::tridiagonal;
using sparse::zero;

//! @test Verifies a product computes exactly its elementary products, far fewer
//! than the dense product's for the same shape.
static_assert(total::multiplications(storage<dense<3, 3>>{} *
                                     storage<dense<3, 3>>{}) == 27);
static_assert(total::multiplications(storage<tridiagonal<3>>{} *
                                     storage<tridiagonal<3>>{}) == 17);
static_assert(total::multiplications(storage<diagonal<3>>{} *
                                     storage<dense<3, 3>>{}) == 9);
static_assert(total::multiplications(storage<diagonal<3>>{} *
                                     storage<diagonal<3>>{}) == 3);
static_assert(total::multiplications(storage<zero<3, 3>>{} *
                                     storage<dense<3, 3>>{}) == 0);

//! @test Verifies the count scales with the stored elements, not the shape: a
//! diagonal product is linear, where a dense product is cubic in the size.
static_assert(total::multiplications(storage<diagonal<32>>{} *
                                     storage<diagonal<32>>{}) == 32);
static_assert(total::multiplications(storage<dense<8, 8>>{} *
                                     storage<dense<8, 8>>{}) == 512);

//! @test Verifies a sum only adds the elements stored by both terms: an element
//! stored by one term is copied.
static_assert(total::additions(storage<diagonal<3>>{} +
                               storage<diagonal<3>>{}) == 3);
static_assert(total::additions(storage<diagonal<3>>{} +
                               storage<tridiagonal<3>>{}) == 3);
static_assert(total::additions(storage<tridiagonal<3>>{} +
                               storage<diagonal<3>>{}) == 3);
static_assert(total::additions(storage<tridiagonal<3>>{} +
                               storage<zero<3, 3>>{}) == 0);
static_assert(total::additions(storage<zero<3, 3>>{} +
                               storage<tridiagonal<3>>{}) == 0);
static_assert(total::additions(storage<dense<3, 3>>{} +
                               storage<dense<3, 3>>{}) == 9);

inline constexpr auto first_row{
    sparse::make_pattern<3, 3>({{0, 0}, {0, 1}, {0, 2}})};
inline constexpr auto last_row{
    sparse::make_pattern<3, 3>({{2, 0}, {2, 1}, {2, 2}})};

static_assert(total::additions(storage<first_row>{} + storage<last_row>{}) ==
              0);

//! @test Verifies a typed matrix adds no arithmetic to its storage's.
inline constexpr auto typed_product{
    sparse_matrix<counted, tridiagonal<3>,
                  std::tuple<counted, counted, counted>,
                  std::tuple<counted, counted, counted>>{} *
    sparse_matrix<counted, diagonal<3>, std::tuple<counted, counted, counted>,
                  std::tuple<counted, counted, counted>>{}};

static_assert(total::multiplications(typed_product.data()) == 7);
} // namespace
} // namespace fcarouge::test
