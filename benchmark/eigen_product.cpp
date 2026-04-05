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

#include <fstream>
#include <iostream>
#include <print>
#include <random>

#include <nanobench.h>

#include <fcarouge/eigen.hpp>

namespace fcarouge::benchmark {
namespace {
void bench() {
  constexpr std::size_t N{16};

  eigen::matrix<double, N, N> a;
  eigen::matrix<double, N, N> b;
  
//   std::vector<double> storage_a(N * N);
//   std::vector<double> storage_b(N * N);
//   std::vector<double> storage_c(N * N);

  std::random_device device;
  std::mt19937 generator{device()};
  std::uniform_real_distribution<> distribution{0., 1.};

  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
    a(i, j) = distribution(generator);
    b[i, j] = distribution(generator);
    }
  }

//   using matrix = std::mdspan<double, std::extents<std::size_t, N, N>>;

//   matrix a{storage_a.data()};
//   matrix b{storage_b.data()};
//   matrix c{storage_c.data()};

//   std::ofstream results{"benchmark/matrix_product_mdspan.csv", std::ios::app};

  ankerl::nanobench::Bench()
      // .output(nullptr)
      .name("eigen_product")
      .run([&]() { 
        eigen::matrix<double, N, N> r = a * b;
        ankerl::nanobench::doNotOptimizeAway(r);
      })
      // .render(csv.c_str(), results)
      ;
}
} // namespace
} // namespace fcarouge::benchmark

int main() { fcarouge::benchmark::bench(); }
