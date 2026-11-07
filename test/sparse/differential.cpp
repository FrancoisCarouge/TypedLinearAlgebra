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
#include <cassert>
#include <cstddef>

namespace fcarouge::test {
namespace {
//! @brief The verification of the operations against the textbook dense
//! computations, of explicit loops over arrays.
struct differential {
  template <std::size_t Rows, std::size_t Columns>
  using reference = std::array<std::array<double, Columns>, Rows>;

  //! @brief A matrix of the pattern filled with distinct, signed, integral
  //! values of exact products and sums.
  template <auto Pattern>
  [[nodiscard]] static constexpr auto fill(double seed)
      -> sparse::matrix<double, Pattern> {
    sparse::matrix<double, Pattern> result;

    for (std::size_t i{0}; i < Pattern.count; ++i) {
      const double value{seed + static_cast<double>(i)};

      result.values.at(i) = i % 2 == 0 ? value : -value;
    }

    return result;
  }

  //! @brief The matrix as an explicit dense array, read at each position.
  template <auto Pattern>
  [[nodiscard]] static constexpr auto
  to_reference(const sparse::matrix<double, Pattern> &value)
      -> reference<Pattern.rows, Pattern.columns> {
    reference<Pattern.rows, Pattern.columns> result{};

    for (std::size_t row{0}; row < Pattern.rows; ++row) {
      for (std::size_t column{0}; column < Pattern.columns; ++column) {
        result.at(row).at(column) = value[row, column];
      }
    }

    return result;
  }

  //! @brief The sum, the difference, the transposition, and the negation.
  template <auto Lhs, auto Rhs>
    requires(Lhs.rows == Rhs.rows and Lhs.columns == Rhs.columns)
  static constexpr void verify_additive() {
    const auto lhs{fill<Lhs>(3.)};
    const auto rhs{fill<Rhs>(11.)};
    const auto lhs_values{to_reference(lhs)};
    const auto rhs_values{to_reference(rhs)};

    const auto sum{to_reference(lhs + rhs)};
    const auto difference{to_reference(lhs - rhs)};
    const auto negated{to_reference(-lhs)};
    const auto mirrored{to_reference(lhs.transpose())};

    for (std::size_t row{0}; row < Lhs.rows; ++row) {
      for (std::size_t column{0}; column < Lhs.columns; ++column) {
        [[maybe_unused]] const double left{lhs_values.at(row).at(column)};
        [[maybe_unused]] const double right{rhs_values.at(row).at(column)};

        assert(sum.at(row).at(column) == left + right);
        assert(difference.at(row).at(column) == left - right);
        assert(negated.at(row).at(column) == -left);
        assert(mirrored.at(column).at(row) == left);
      }
    }
  }

  //! @brief The product, and its pattern exactly storing the reached positions.
  template <auto Lhs, auto Rhs>
    requires(Lhs.columns == Rhs.rows)
  static constexpr void verify_product() {
    const auto lhs{fill<Lhs>(3.)};
    const auto rhs{fill<Rhs>(11.)};
    const auto lhs_values{to_reference(lhs)};
    const auto rhs_values{to_reference(rhs)};
    const auto product{lhs * rhs};
    const auto values{to_reference(product)};

    for (std::size_t row{0}; row < Lhs.rows; ++row) {
      for (std::size_t column{0}; column < Rhs.columns; ++column) {
        double expected{0.};
        bool reached{false};

        for (std::size_t inner{0}; inner < Lhs.columns; ++inner) {
          expected +=
              lhs_values.at(row).at(inner) * rhs_values.at(inner).at(column);
          reached =
              reached || (Lhs.stores(row, inner) && Rhs.stores(inner, column));
        }

        assert(values.at(row).at(column) == expected);
        assert(decltype(product)::structure.stores(row, column) == reached);
      }
    }
  }

  //! @brief A typed matrix computes the same as its storage, and as a typed
  //! matrix without structural zeros.
  template <auto Lhs, auto Rhs> static constexpr void verify_typed() {
    using dense_matrix =
        sparse_matrix<double, sparse::dense<Lhs.rows, Lhs.columns>>;

    const sparse_matrix<double, Lhs> lhs{fill<Lhs>(3.)};
    const sparse_matrix<double, Rhs> rhs{fill<Rhs>(11.)};

    assert((lhs + rhs).data() == lhs.data() + rhs.data());
    assert((lhs - rhs).data() == lhs.data() - rhs.data());
    assert((lhs * rhs).data() == lhs.data() * rhs.data());
    assert(transposed(lhs).data() == lhs.data().transpose());
    assert((-lhs).data() == -lhs.data());
    assert((lhs * 2.).data() == lhs.data() * 2.);
    assert((lhs + rhs) == (dense_matrix{lhs} + rhs));
    assert((lhs * rhs) == (dense_matrix{lhs} * rhs));
  }

  template <auto Lhs, auto... Patterns> static constexpr void verify_with() {
    (verify_additive<Lhs, Patterns>(), ...);
    (verify_product<Lhs, Patterns>(), ...);
    (verify_typed<Lhs, Patterns>(), ...);
  }

  template <auto... Patterns> static constexpr void verify_all_pairs() {
    (verify_with<Patterns, Patterns...>(), ...);
  }
};

inline constexpr auto is_arrow{[](std::size_t row, std::size_t column) -> bool {
  return row == 0 || column == 0 || row == column;
}};
inline constexpr auto arrow{sparse::generate<3, 3, is_arrow>()};
inline constexpr auto anti_diagonal{
    sparse::make_pattern<3, 3>({{0, 2}, {1, 1}, {2, 0}})};
inline constexpr auto corner{sparse::make_pattern<3, 3>({{2, 2}})};

//! @test Verifies every pair of a subset of the patterns at compile time,
//! within the compilers' constant evaluation limits.
static_assert([] -> bool {
  differential::verify_all_pairs<sparse::tridiagonal<3>, sparse::zero<3, 3>,
                                 arrow, corner>();

  return true;
}());

//! @test Verifies every pair of a family of patterns, and rectangular factors,
//! at run time.
[[maybe_unused]] const auto test{[] -> int {
  // Not a constant expression: every operation below is `constexpr`, and this
  // test would otherwise be evaluated at compile time, and never at run time.
  [[maybe_unused]] const volatile int run_time{0};

  assert(run_time == 0);

  differential::verify_all_pairs<
      sparse::diagonal<3>, sparse::tridiagonal<3>, sparse::lower_triangular<3>,
      sparse::upper_triangular<3>, sparse::dense<3, 3>, sparse::zero<3, 3>,
      arrow, anti_diagonal, corner>();

  constexpr auto wide{
      sparse::make_pattern<2, 4>({{0, 0}, {0, 3}, {1, 1}, {1, 2}})};
  constexpr auto tall{
      sparse::make_pattern<4, 3>({{0, 2}, {1, 0}, {2, 1}, {3, 0}, {3, 2}})};

  differential::verify_product<wide, tall>();
  differential::verify_product<wide, sparse::dense<4, 3>>();
  differential::verify_product<sparse::dense<2, 4>, tall>();
  differential::verify_product<wide, sparse::zero<4, 3>>();
  differential::verify_additive<wide, sparse::dense<2, 4>>();
  differential::verify_additive<wide, sparse::zero<2, 4>>();

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
