/* Typed Linear Algebra
Version 0.1.0
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
#include <type_traits>

namespace fcarouge::test {
namespace {
//! @test Verifies the nested type of the matrix.
[[maybe_unused]] auto test{[] {
  matrix<double, 3, 3> z;

  static_assert(
      std::same_as<
          decltype(z),
          fcarouge::typed_matrix<
              fcarouge::typed_matrix<
                  fcarouge::typed_matrix<Eigen::Matrix<double, 3, 3, 0, 3, 3>,
                                         std::tuple<double, double, double>,
                                         std::tuple<double, double, double>>,
                  std::tuple<double, double, double>,
                  std::tuple<double, double, double>>,
              std::tuple<double, double, double>,
              std::tuple<double, double, double>>>);

  static_assert(std::same_as<
                decltype(z.data()),
                fcarouge::typed_matrix<
                    fcarouge::typed_matrix<Eigen::Matrix<double, 3, 3, 0, 3, 3>,
                                           std::tuple<double, double, double>,
                                           std::tuple<double, double, double>>,
                    std::tuple<double, double, double>,
                    std::tuple<double, double, double>> &>);

  static_assert(std::same_as<
                decltype(z)::matrix,
                fcarouge::typed_matrix<
                    fcarouge::typed_matrix<Eigen::Matrix<double, 3, 3, 0, 3, 3>,
                                           std::tuple<double, double, double>,
                                           std::tuple<double, double, double>>,
                    std::tuple<double, double, double>,
                    std::tuple<double, double, double>>>);

  static_assert(std::same_as<decltype(z)::underlying, double>);

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
