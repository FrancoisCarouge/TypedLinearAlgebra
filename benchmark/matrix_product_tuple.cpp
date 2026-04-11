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
#include <linalg>
#include <mdspan>
#include <print>
#include <random>
#include <tuple>
#include <utility>

#include <nanobench.h>

namespace fcarouge::benchmark {
namespace {

const std::string csv{
    R"(# Function, Elapsed Time (s)
{{{{#result}}}}{{{{name}}}}, {{{{minimum(elapsed)}}}}{{{{/result}}}}
)"};

template <typename T, std::size_t... Is>
auto make_n_tuple_impl(std::index_sequence<Is...>) {
  // Expansion of (Is, T) will effectively result in T being repeated N times
  return std::tuple<std::decay_t<decltype((void)Is, T())>...>{};
}

template <typename T, std::size_t N>
using n_tuple = decltype(make_n_tuple_impl<T>(std::make_index_sequence<N>{}));

namespace mp {

namespace detail {

// Deduce the return type from I = 0
template <class F>
using result_t =
    decltype(std::declval<F>()(std::integral_constant<std::size_t, 0>{}));

template <std::size_t I, std::size_t N> struct with_index_impl {
  template <class F> static constexpr result_t<F> call(std::size_t i, F &&f) {
    if (i == I) {
      return std::forward<F>(f)(std::integral_constant<std::size_t, I>{});
    } else {
      return with_index_impl<I + 1, N>::call(i, std::forward<F>(f));
    }
  }
};

// Base case
template <std::size_t N> struct with_index_impl<N, N> {
  template <class F> static constexpr result_t<F> call(std::size_t, F &&) {
    throw std::out_of_range("mp_with_index: index out of range");
  }
};

} // namespace detail

template <std::size_t N, class F>
constexpr auto mp_with_index(std::size_t i, F &&f) -> detail::result_t<F> {
  return detail::with_index_impl<0, N>::call(i, std::forward<F>(f));
}

} // namespace mp

template <typename T> struct accessor {
  using element_type = double;
  using reference = double &;
  using data_handle_type = T *;
  using offset_policy = accessor;

  constexpr accessor() = default;

  constexpr reference access(data_handle_type p, std::size_t i) const noexcept {
    constexpr std::size_t N = std::tuple_size_v<std::remove_reference_t<T>>;

    return mp::mp_with_index<N>(
        i, [&](auto I) -> double & { return std::get<I>(*p); });
  }

  constexpr data_handle_type offset(data_handle_type p,
                                    std::size_t i) const noexcept {
    return p;
  }
};

//! @brief Benchmark for matrix product with mdspan.
//!
//! @details This benchmark demonstrates the performance of matrix
//! multiplication under `std::linalg::matrix_product` using `std::mdspan`. It
//! initializes two matrices of size 32x32 with random values and measures the
//! time taken to compute their product with the nanobench library. The results
//! are printed to the console.
void bench() {
  constexpr std::size_t N{1};

  using storage = n_tuple<double, N * N>;

  storage storage_a;
  storage storage_b;
  storage storage_c;

  std::random_device device;
  std::mt19937 generator{device()};
  std::uniform_real_distribution<> distribution{0., 1.};

  std::apply([&](auto &...value) { ((value = distribution(generator)), ...); },
             storage_a);

  std::apply([&](auto &...value) { ((value = distribution(generator)), ...); },
             storage_b);

  using matrix = std::mdspan<double, std::extents<std::size_t, N, N>,
                             Kokkos::layout_right, accessor<storage>>;

  matrix a{&storage_a};
  matrix b{&storage_b};
  matrix c{&storage_c};

  std::ofstream results{"benchmark/matrix_product_mdspan.csv", std::ios::app};

  ankerl::nanobench::Bench()
      // .output(nullptr)
      .name("matrix_product_tuple")
      .run([&]() { std::linalg::matrix_product(a, b, c); })
      // .render(csv.c_str(), results)
      ;
}
} // namespace
} // namespace fcarouge::benchmark

int main() { fcarouge::benchmark::bench(); }
