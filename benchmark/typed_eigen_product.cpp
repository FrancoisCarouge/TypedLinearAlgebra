/* Typed Linear Algebra
Version 0.2.0
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

#include <nanobench.h>

#include <cstddef>
#include <format>
#include <fstream>
#include <random>
#include <string>

namespace fcarouge::benchmark {
namespace {
template <auto Size>
const std::string csv{std::format(
    "{{{{#result}}}}| {{{{title}}}} | {:5d}x{:<5d} | {{{{median(elapsed)}}}} | "
    "{{{{medianAbsolutePercentError(elapsed)}}}} |{{{{/result}}}}\n",
    Size, Size)};

//! @benchmark Typed Eigen square matrix-matrix product.
template <auto Size> void bench() {
  matrix<double, Size, Size> a;
  matrix<double, Size, Size> b;
  std::random_device device;
  std::mt19937 generator{device()};
  std::uniform_real_distribution<> distribution{0., 1.};

  for (std::size_t i{0}; i < Size; ++i) {
    for (std::size_t j{0}; j < Size; ++j) {
      a(i, j) = distribution(generator);
      b(i, j) = distribution(generator);
    }
  }

  std::ofstream results{"results.txt", std::ios::app};
  ankerl::nanobench::Bench()
      .output(nullptr)
      .title("typed matrix from Eigen::Matrix")
      .run([&]() {
        matrix<double, Size, Size> r{a * b};
        ankerl::nanobench::doNotOptimizeAway(r);
      })
      .render(csv<Size>.c_str(), results);
}
} // namespace
} // namespace fcarouge::benchmark

int main() { fcarouge::benchmark::bench<${SIZE}>(); }
