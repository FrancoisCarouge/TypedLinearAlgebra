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

#include "fcarouge/realtime.hpp"

#include <Eigen/Eigen>
#include <nanobench.h>

#include <format>
#include <fstream>
#include <random>
#include <string>

namespace fcarouge::benchmark {
namespace {
//! @brief The benchmark measures, prints, and tears down its counters past its
//! entry point: it opts out of the real-time verification for its lifetime.
const not_realtime opt_out;

//! @brief Nanoseconds in a second, nanobench measuring elapsed time in seconds.
constexpr double nanoseconds_per_second{1e9};

//! @brief Percent in a ratio, nanobench reporting its median absolute
//! percentage error as a ratio.
constexpr double percent_per_ratio{100.};

//! @benchmark Eigen square matrix-matrix product.
template <auto Size> void bench() {
  Eigen::Matrix<double, Size, Size> a;
  Eigen::Matrix<double, Size, Size> b;
  std::random_device device;
  std::mt19937 generator{device()};
  std::uniform_real_distribution<> distribution{0., 1.};

  for (Eigen::Index i{0}; i < Size; ++i) {
    for (Eigen::Index j{0}; j < Size; ++j) {
      a(i, j) = distribution(generator);
      b(i, j) = distribution(generator);
    }
  }

  ankerl::nanobench::Bench runner;
  runner.output(nullptr).title("Eigen::Matrix").run([&]() {
    Eigen::Matrix<double, Size, Size> r{a * b};
    ankerl::nanobench::doNotOptimizeAway(r);
  });

  const auto &result{runner.results().front()};
  std::ofstream results{"results.txt", std::ios::app};
  results << std::format(
      "| {} | {}x{} | {:.2f} | {:.2f} |\n", runner.title(), Size, Size,
      result.median(ankerl::nanobench::Result::Measure::elapsed) *
          nanoseconds_per_second,
      result.medianAbsolutePercentError(
          ankerl::nanobench::Result::Measure::elapsed) *
          percent_per_ratio);
}
} // namespace
} // namespace fcarouge::benchmark

int main() { fcarouge::benchmark::bench<${SIZE}>(); }
